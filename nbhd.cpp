#include <algorithm>
#include <vector>
#include <iostream>
#include <signal.h>

extern "C" {
#include "premier.h"
#include "sort.h"
#include "kvec.h"
#include "kmer.h"
#include "numeric.h"
#include "trans.h"
#include "em.h"
#include "bitarray.h"
#include "read.h"
}

#define PMR_SPARSE_NBHD_SIZE 16

static void hmm_build_1st_nbhd_fast(data *d, read_t *read, int tmin, std::vector<kmer_t*> &k_nbhd);

static int kmer_suffix_comp_qsort(const void *a, const void *b)
{
	register kmer_t *ka = *((kmer_t **) a);
	register kmer_t *kb = *((kmer_t **) b);
	register kbits_t sfx_a = kbits_id_of(ka), sfx_b = kbits_id_of(kb);

	return (sfx_a > sfx_b) ? 1 : -1;
}

static int kmer_suffix_comp(kmer_t* ka, kmer_t* kb)
{
	register kbits_t sfx_a = kbits_id_of(ka), sfx_b = kbits_id_of(kb);

	return (sfx_a > sfx_b);
}

static void hmm_build_1st_nbhd_fast(data *d, read_t *read, int tmin, std::vector<kmer_t*> &k_nbhd)
{
	static char BASES[5] = {'A', 'C', 'T', 'G', 'N'};

	int actg_trans_count[4];
	int cum_actg_trans_count[4];
	int total_hmm_trans = 0;

	const int kmer_len = d->opt->kmer_length;
	const int num_kmers = read->length - kmer_len + 1;
	const kbits_t klen_flag = (1UL << (kmer_len << 1)) - 1;
	const int shift = (kmer_len - 1) << 1;
	const int qmax = d->opt->max_qual_score + 1;
	const int qmax2 = qmax << 1;
	const int bwidth = d->opt->qscore_bin_width;
	const double _uniq_tp = 0.0;

	int tmax = read->length - kmer_len;
	if (tmin >= tmax) return; 

	// starting position on reverse complement
	tmin = read->length - 1 - (tmin + kmer_len - 1);

	mempool_t *mp = mempool_create(MEMPOOL_DEFAULT_SIZE);

	dict *kdict = d->kmer_dict;

	/* reverse complement sequence */
	char *rseq = (char *) calloc(read->length + 1, sizeof(char));
	for (int p = 0; p < read->length; ++p) {
		rseq[p] = BASES[
			base_to_rc_bits(read->sequence[read->length - p - 1])];
	}

	char *conv_q = convert_quality_scores(read->qscore, read->length, 
			d->opt->qual_score_offset, bwidth, mp);

	state_nbhd_t *T_nbhds = (state_nbhd_t *) mempool_nalloc(mp,
			sizeof(*T_nbhds) * num_kmers, 16);

	/* flag for N (not determined) base, first kmer only */
	kbits_t obs_n_flag = kmer_n_base_flag(rseq, kmer_len);
	kbits_t next_obs_n_flag = kmer_n_base_flag(rseq + 1, kmer_len);

	/* convert string literal to numeric representation */
	kbits_t observed_kmers[BITS_TO_U64(read->length << 1)];
	read_seq_to_numeric(rseq, observed_kmers, read->length, kmer_len);

	/* set up first kmer neighborhood */
	kbits_t obs_kid = bitarray_get(observed_kmers, tmin << 1, kmer_len << 1);
	ksubstr_t obs_base = kmer_effective_base(rseq[tmin+kmer_len-1]);
	ksubstr_t obs_next_base = kmer_effective_base(rseq[tmin+kmer_len]);

	// set the first neighborhood to be the observed kmer itself
	T_nbhds[tmin].size = hmm_load_preconstructed_nbhd(
			d, read->id, obs_kid, 
			obs_n_flag, conv_q, &T_nbhds[tmin], kmer_len, 0, kdict, mp, 
			&_uniq_tp);

	/* initial state distribution */
	for (int t = tmin; t < tmax; t++) {
		int tnext_kp = t + kmer_len;

		obs_next_base = kmer_effective_base(rseq[t + kmer_len]);
		obs_kid = bitarray_get(observed_kmers, t << 1, kmer_len << 1);
		kbits_t obs_next_kid = bitarray_get(observed_kmers, (t + 1) << 1, 
				kmer_len << 1);

		kmer_t *obs_pk = (kmer_t *) &(dict_find(kdict,
					kbits_cast_to_ptr(obs_kid))->value);

		kbits_t _new_nf = ((obs_next_base >> 2) << 1) | 
			(obs_next_base >> 2);
		next_obs_n_flag = (next_obs_n_flag >> 2) | (_new_nf << shift);

		state_nbhd_t *curr_nbhd = &T_nbhds[t],
					 *next_nbhd = &T_nbhds[t+1];

		const int dmax = d->opt->max_hamming_dist;

		int t_nbhd_size = curr_nbhd->size;

		kbits_t *t_states_sorted_sfx = T_nbhds[t].states_sorted_sfx;
		double *t_states_alphas = T_nbhds[t].alphas;
		double **t_kmer_trans_p = T_nbhds[t].kmer_trans_prob;

		/* ----- II. compute one iteration of forward algorithm ----- */

		/* --- II a. set up kmer neighborhood at next position --- */
		memset(actg_trans_count, 0, sizeof(int) << 2);

		/* one linear scan to determine # of distinct suffix, 
		 * compute downstream actg nucleotide counts etc. */
		int n_distinct_suffixes = hmm_count_distinct_suffixes(curr_nbhd, 
				actg_trans_count, dmax, obs_kid, obs_n_flag,
				obs_next_base, d->transition_p);

		curr_nbhd->n_uniq_sfx = n_distinct_suffixes;

		/* number of kmers in next position */
		int tnext_nbhd_size = actg_trans_count[0];
		/* cumulative transition counts */
		cum_actg_trans_count[0] = actg_trans_count[0];
		for (int i = 1; i < 4; i++) {
			tnext_nbhd_size += actg_trans_count[i];
			cum_actg_trans_count[i] = cum_actg_trans_count[i-1] +
				actg_trans_count[i];
		}

		if (tnext_nbhd_size == 0) {
			mempool_destroy(mp);
			return;
		}

		/* allocate memory given tnext_nbhd_size */
		char *pmem = (char *) mempool_nalloc(mp,
				nbhd_alloc_size(tnext_nbhd_size), 16);

		next_nbhd->size = tnext_nbhd_size;
		hmm_setup_nbhd_ptrs(next_nbhd, tnext_nbhd_size, pmem);

		int kmer_idx = 0, index_pfx = 0;
		register ksubstr_t mismatch_flag = ~(1UL << obs_next_base) |
			~((obs_next_base >> 2) - 1UL);
		for (int i = 0; i < n_distinct_suffixes; i++) {
			/*
			int n_common_sfx_kmers = bitarray_get_pwr2(
					curr_nbhd->ba_distinct_sfx, i, 1) + 1;
			*/
			kbits_t _repr_kid = t_states_sorted_sfx[kmer_idx];
			//kbits_t mut_flag = ((_repr_kid ^ obs_kid) | obs_n_flag) & 3UL;
			kbits_t sfx_hd = bitarray_get_pwr2(curr_nbhd->ba_hamming_dist, 
					i, 3);
			kbits_t common_sfx = kbits_suffix(_repr_kid);

			/* v 4 bits                   v 4 bits (higher)
			 *
			 * A  -                       - A
			 * C  |  -- k-1 substring --  | C
			 * T  |                       | T
			 * G  -                       - G
			 */
			kbits_t _tf = bitarray_get_pwr2(curr_nbhd->ba_pfx_sfx_flags,
					i, 3);
			kbits_t ups_trans_packed = trans_flag2packed(_tf & 15);
			/* which kmers/nucleotides current suffix transition into */
			kbits_t dns_trans_packed = trans_flag2packed(_tf >> 4);

			int n_common_sfx_kmers = ups_trans_packed & 7; 

			int j = dns_trans_packed & 7;
			for (dns_trans_packed >>= 3; j > 0; --j, dns_trans_packed >>= 2) {
				kbits_t dns_base = dns_trans_packed & 3;
				kbits_t tnext_kid = common_sfx | (dns_base << shift);

				/* compute the position sorted by suffix */
				int index_sfx = cum_actg_trans_count[dns_base] - 
					actg_trans_count[dns_base];
				actg_trans_count[dns_base]--;

				kmer_t *tnext_kmer = (kmer_t *) &(dict_find(kdict, 
						kbits_cast_to_ptr(tnext_kid))->value);
				next_nbhd->kmer_ptrs[index_sfx] = tnext_kmer;

				/* set up transition flag */
				bitarray_set_pwr2(next_nbhd->ba_kmer_trans_flag,
						index_sfx, tnext_kmer->trans_flag, 2);
				/* set up Hamming distance */
				bitarray_set_pwr2(next_nbhd->ba_hamming_dist, 
						index_sfx, 
						sfx_hd + ((mismatch_flag >> dns_base) & 1), 3);

				next_nbhd->suffixes_order[index_pfx] = index_sfx;
				next_nbhd->states_sorted_sfx[index_sfx] = tnext_kid;
				index_pfx++;
			}

			kmer_idx += n_common_sfx_kmers;
		}

		obs_base = obs_next_base;
		obs_n_flag = next_obs_n_flag;
	}

	state_nbhd_t *tmax_nbhd = &T_nbhds[tmax];

	// find the reverse complements of all kmers in tmax_nbhd
	for (int i = 0; i < tmax_nbhd->size; ++i) {
		kbits_t rev_kid = kmer_reverse_complement(
				tmax_nbhd->states_sorted_sfx[i], kmer_len);
		kmer_t *rev_pk = (kmer_t *) &(dict_find(d->kmer_dict,
					kbits_cast_to_ptr(rev_kid))->value);
		tmax_nbhd->kmer_ptrs[i] = rev_pk;
	}

	if (tmax_nbhd->size > 0) {
		std::copy(tmax_nbhd->kmer_ptrs,
				tmax_nbhd->kmer_ptrs + tmax_nbhd->size,
				std::back_inserter(k_nbhd));
		std::sort(k_nbhd.begin(), k_nbhd.end(), kmer_suffix_comp);
	}

destroy:
	mempool_destroy(mp);
}

void hmm_build_1st_nbhd_slow(data *d, read_t *read, void *fdata)
{
	const int kmer_len = d->opt->kmer_length;
	const double thresh = d->opt->penalty_eta / log1p(1/d->opt->penalty_gamma);

	const int dmax = 4;
	const int max_tries = kmer_len * (kmer_len - 1) / 2 * 9;

	int *rsel = (int *) fdata;
	if (rsel != NULL && rsel[read->id] == 0) return;
	mempool_t *_mp = rsel == NULL ? d->mp : d->tmpmp;

	// # of mutated kmers examined/explored at each certain d level.
	int d_cands_explored[dmax] = {0};
	int d_max_tries[dmax+1] = {0};
	d_max_tries[2] = max_tries;
	for (int d = 3; d <= dmax; ++d) {
		d_max_tries[d] = d_max_tries[d-1] * 3;
	}

	int *qorder = sort_order_char(read->qscore, kmer_len, NULL);

	kbits_t obs_kid = kmer_seq_to_id(read->sequence, kmer_len);
	kmer_t *pk = (kmer_t *) &(dict_find(d->kmer_dict, 
				kbits_cast_to_ptr(obs_kid))->value);

	/* # of times a kmer and its reverse complement are observed */
	double k_cexp = EXP(d->avg_initial_p[kmer_initial_index(pk,
			d->uniq_kmer_init_off)]) * d->total_kmer_exp_counts;

	kvec_t(kmer_t *) k_nbhd;
	kv_init(k_nbhd);

	/* observed kmer has zero Hamming distance */
	mut_kmer_t obs_mk = {obs_kid, 0};

	/* observed kmer is always in the first neighborhood */
	kv_push(kmer_t *, k_nbhd, pk);

	std::vector<mut_kmer_t> k_mut_cands;

	/* expand the observed kmer if :
	 *   1) it has # occurence below threshold, or
	 *   2) it has no eligible transitions (removed in initialization) */
	if (k_cexp < thresh || pk->trans_flag == 0) {
		k_mut_cands.reserve(d_max_tries[dmax]);
		k_mut_cands.push_back(obs_mk);
	}

	int hamming_d = 0; /* hamming distance */
	int counter = 0;
	int mutcand_idx = 0;
	while (mutcand_idx < k_mut_cands.size()) {
		mut_kmer_t mk = k_mut_cands[mutcand_idx]; //k_mut_cands.front();
		int _d = _mm_popcnt_u64(mk.mut_flag);

		++mutcand_idx;
		//k_mut_cands.pop_front();

		++d_cands_explored[_d];
		if ((_d <= 2 && _d > hamming_d) ||
				(_d > 2 && d_cands_explored[_d] >= d_max_tries[_d])) {
		//if (_d > hamming_d) {
			/* Hamming distance increased, check if we already have viable
			 * candidates in the neighborhood */
			if (kv_size(k_nbhd) > 1) break;
			hamming_d = _d;
		}

		for (int j = 0; j < kmer_len; j++) {
			int dfs_flag  = 0;

			/* prioritize the loci with lower quality scores */
			int i = qorder[j];
			int i2 = i << 1;

			/* to avoid duplicated mutations, we will only mutate position j,
			 * if no bits higher than j are set to 1. */
			if (mk.mut_flag & ((1UL << (j+1)) - 1UL)) continue;

			for (kbits_t base = 0; base < 4; ++base) {
				if (base == ((mk.id >> i2) & 3)) continue;
				if (base == ((obs_kid >> i2) & 3)) continue;

				kbits_t mut_kid = (mk.id & ~(3UL << i2)) | (base << i2);
				kbits_t mut_flag = mk.mut_flag | (1UL << j);

				dictentry *de_mut = dict_find(d->kmer_dict, 
						kbits_cast_to_ptr(mut_kid));

				if (de_mut == NULL && (_d+1) < dmax) {
					/* continue searching */
					mut_kmer_t new_mk = {mut_kid, mut_flag};

					// breadth-first search for smaller d, due to the
					// relatively low complexity; for larger d, switch to
					// depth-first search to be more opportunistic (and greedy)
					if (_d <= 2) 
						k_mut_cands.push_back(new_mk);
					else {
						--mutcand_idx;
						k_mut_cands[mutcand_idx] = new_mk;
						//k_mut_cands.push_front(new_mk);
						dfs_flag = 1;
					}
				}
				else if (de_mut != NULL) {
					kmer_t *mut_pk = (kmer_t *) &(de_mut->value);

					double mut_k_cexp = EXP(d->avg_initial_p[
							kmer_initial_index(mut_pk, d->uniq_kmer_init_off)]) * 
						d->total_kmer_exp_counts;

					if (mut_k_cexp > thresh && mut_pk->trans_flag != 0) {
						kv_push(kmer_t *, k_nbhd, mut_pk);
					}
					else if ((_d + 1) < dmax) {
						mut_kmer_t new_mk = {mut_kid, mut_flag};
						if (_d <= 2)
							k_mut_cands.push_back(new_mk);
						else {
							--mutcand_idx;
							k_mut_cands[mutcand_idx] = new_mk;
							//k_mut_cands.push_front(new_mk);
							dfs_flag = 1;
						}
					}
				}
			}

			if (dfs_flag) break;
		}
	}

	/* sort the nbhd by suffix */
	register int size_nbhd = kv_size(k_nbhd);
	kmer_t **new_nbhd = NULL;

#pragma omp critical
{
	new_nbhd = (kmer_t **) mempool_alloc(_mp, 
			size_nbhd * sizeof(kmer_t *));
	memcpy(new_nbhd, &(kv_A(k_nbhd, 0)), size_nbhd * sizeof(kmer_t *));
}

	qsort(new_nbhd, size_nbhd, sizeof(kmer_t *), kmer_suffix_comp_qsort);

	d->preconstructed_nbhds[read->id] = (kmer_nbhd_t) {size_nbhd, new_nbhd};

	kv_destroy(k_nbhd);
	free(qorder);
}

void hmm_build_1st_nbhd(data *d, read_t *read, void *fdata)
{
	const int kmer_len = d->opt->kmer_length;
	const double thresh = d->opt->penalty_eta / log1p(1/d->opt->penalty_gamma);
	const int qoff = d->opt->qual_score_offset;

	const double valid_avg_qscore = 15.0;

	int tmax = read->length - kmer_len;

	kbits_t obs_kid = kmer_seq_to_id(read->sequence, kmer_len);
	kmer_t *pk = (kmer_t *) &(dict_find(d->kmer_dict, 
				kbits_cast_to_ptr(obs_kid))->value);

	int *rsel = (int *) fdata;
	if (rsel != NULL && rsel[read->id] == 0) return;
	mempool_t *_mp = rsel == NULL ? d->mp : d->tmpmp;

	/* # of times a kmer and its reverse complement are observed */
	double k_cexp = EXP(d->avg_initial_p[kmer_initial_index(pk,
			d->uniq_kmer_init_off)]) * d->total_kmer_exp_counts;

	std::vector<kmer_t *> k_nbhd;

	double sum_qscore = 0.0;
	for (int t = 0; t < kmer_len; ++t) {
		sum_qscore += read->qscore[t] - qoff; 
	}

	if (k_cexp < thresh || pk->trans_flag == 0 ||
			(sum_qscore / kmer_len) < valid_avg_qscore) {

		int tvalid = 1;
		do {
			bool valid_kmer_located = false;
			for (; tvalid < tmax; ++tvalid) {
				sum_qscore += (read->qscore[tvalid + kmer_len - 1] - 
						read->qscore[tvalid - 1]);
				kbits_t _rkid = kmer_reverse_complement(
						kmer_seq_to_id(read->sequence + tvalid, kmer_len), 
						kmer_len);
				kmer_t *_rpks = (kmer_t *) &(dict_find(d->kmer_dict, 
							kbits_cast_to_ptr(_rkid))->value);

				double avg_k_cexp = EXP(d->avg_initial_p[kmer_initial_index(_rpks,
							d->uniq_kmer_init_off)]) * d->total_kmer_exp_counts;

				if (avg_k_cexp > thresh && _rpks->trans_flag > 0 &&
						(sum_qscore / kmer_len) >= valid_avg_qscore) {
					valid_kmer_located = true;
					break;
				}
			}

			if (valid_kmer_located) {
				// try a faster method to build the neighborhood.
				hmm_build_1st_nbhd_fast(d, read, tvalid, k_nbhd);
				if (k_nbhd.size() == 0) { // || (k_nbhd.size() == 1 &&
							//k_nbhd[0] == pk)) {
					// use the slower method as a fallback
					//
					//if (pk->trans_flag) k_nbhd.push_back(pk);
					//else {
						return hmm_build_1st_nbhd_slow(d, read, fdata);
					//}
				}
			}
			else {
				// failed to find an "anchor" valid kmer, try the slow method 
				return hmm_build_1st_nbhd_slow(d, read, fdata);
			}
		} while(0);
	}
	else {
		// use the observed kmer 
		k_nbhd.push_back(pk);
	}

	kmer_t **new_nbhd = NULL;
#pragma omp critical
{
	new_nbhd = (kmer_t **) mempool_alloc(_mp, 
			k_nbhd.size() * sizeof(kmer_t *));
	std::copy(k_nbhd.begin(), k_nbhd.end(), new_nbhd);
}

	d->preconstructed_nbhds[read->id] = (kmer_nbhd_t) {(int) k_nbhd.size(), new_nbhd};
}

void hmm_determine_dmax(data *d, read_t *read, void *fdata)
{
	const int kmer_len = d->opt->kmer_length;
	int tmax = read->length - kmer_len + 1;
	double thresh = d->opt->penalty_eta / log1p(1/d->opt->penalty_gamma);
	dict *kdict = d->kmer_dict;

	int *rsel = (int *) fdata;
	if (rsel != NULL && rsel[read->id] == 0) return;

	int is_prev_kcov_low = 0;
	int low_kcov_span_len = 0;
	std::vector<int> low_kcov_regions;
	low_kcov_regions.reserve(read->length >> 1);

	/* convert string literal to numeric representation */
	kbits_t observed_kmers[BITS_TO_U64(read->length << 1)];
	read_seq_to_numeric(read->sequence, observed_kmers, read->length, kmer_len);

	double avg_k_cexp = 0.0;
	double min_k_cexp = 1e6;

	for (int t = 0; t < tmax; t++) {
		kbits_t obs_kid = bitarray_get(observed_kmers, t << 1, kmer_len << 1);
		kmer_t *obs_pk = (kmer_t *) &(dict_find(kdict,
					kbits_cast_to_ptr(obs_kid))->value);

		double obs_k_cexp = EXP(d->avg_initial_p[kmer_initial_index(obs_pk,
			d->uniq_kmer_init_off)]) * d->total_kmer_exp_counts;
		avg_k_cexp += obs_k_cexp;

		if (obs_k_cexp < min_k_cexp) min_k_cexp = obs_k_cexp;

		int low_kcov = (obs_k_cexp < thresh);
		if (low_kcov) {
			++low_kcov_span_len;
		}
		else {
			if (low_kcov_span_len > 0) {
				low_kcov_regions.push_back(low_kcov_span_len);
				low_kcov_span_len = 0;
			}
		}

		is_prev_kcov_low = low_kcov;
	}

	avg_k_cexp /= tmax;
	if (avg_k_cexp > d->opt->kcov_q90 && min_k_cexp > thresh) {
		d->read_dmax[read->id] = d->opt->max_hamming_dist >> 2;
	}
	else if (min_k_cexp > thresh) {
		d->read_dmax[read->id] = d->opt->max_hamming_dist >> 1;
	}
	else {
		d->read_dmax[read->id] = d->opt->max_hamming_dist;
	}
	d->read_dmax[read->id] = d->opt->max_hamming_dist;
}

/**
 * hmm_adaptive_dmax:
 *
 * Choose $d_{i, t}$, the maximum Hamming distance for read $i$ at position
 * $t$.
 *
 * This procedure precedes the EM algorithm, and is crucial to reduce the
 * overall computational complexity of the HMM.
 *
 * It is important to note that the running time of PREMIER is dominated by a
 * tiny fraction of reads which have large neighborhoods. 
 * For the majority of the reads, the neighborhoods are sparse, due to the
 * kmer-uniqueness induced by the l0-like penalty. 
 * For reads with sparse neighborhoods, the complexity of the HMM is
 * linear to the (average) neighborhood size; whereas reads with
 * dense neighborhoods still have linear complexity, albeit with larger scaling
 * constant.
 *
 */
void hmm_adaptive_dmax(data *d, read_t *read, void *fdata) 
{
	static char BASES[5] = {'A', 'C', 'T', 'G', 'N'};

	int actg_trans_count[4];
	int cum_actg_trans_count[4];
	int total_hmm_trans = 0;

	const int skip_rebuild = (fdata != NULL);

	const int kmer_len = d->opt->kmer_length;
	const int num_kmers = read->length - kmer_len + 1;
	const kbits_t klen_flag = (1UL << (kmer_len << 1)) - 1;
	const int shift = (kmer_len - 1) << 1;
	const int qmax = d->opt->max_qual_score + 1;
	const int qmax2 = qmax << 1;
	const int bwidth = d->opt->qscore_bin_width;
	const double _uniq_tp = 0.0;

	int tmin = 0;
	int tmax = read->length - kmer_len;

	mempool_t *mp = mempool_create(MEMPOOL_DEFAULT_SIZE);
	dict *kdict = d->kmer_dict;

	char *rseq = read->sequence; 
	char *conv_q = convert_quality_scores(read->qscore, read->length, 
			d->opt->qual_score_offset, bwidth, mp);

	state_nbhd_t *T_nbhds = (state_nbhd_t *) mempool_nalloc(mp,
			sizeof(*T_nbhds) * num_kmers, 16);

	/* flag for N (not determined) base, first kmer only */
	kbits_t obs_n_flag = kmer_n_base_flag(rseq, kmer_len);
	kbits_t next_obs_n_flag = kmer_n_base_flag(rseq + 1, kmer_len);

	/* convert string literal to numeric representation */
	kbits_t observed_kmers[BITS_TO_U64(read->length << 1)];
	read_seq_to_numeric(rseq, observed_kmers, read->length, kmer_len);

	/* set up first kmer neighborhood */
	kbits_t obs_kid = bitarray_get(observed_kmers, tmin << 1, kmer_len << 1);
	ksubstr_t obs_base = kmer_effective_base(rseq[tmin+kmer_len-1]);
	ksubstr_t obs_next_base = kmer_effective_base(rseq[tmin+kmer_len]);

	// set the first neighborhood to be the observed kmer itself
	T_nbhds[tmin].size = hmm_load_preconstructed_nbhd(
			d, read->id, obs_kid, 
			obs_n_flag, conv_q, &T_nbhds[tmin], kmer_len, 0, kdict, mp, 
			&_uniq_tp);

	// look-ahead window
	int window_size = d->opt->hd_window_size; 

	// use the window size // half of the kmer length as a starting point
	uint8_t dmax = kmer_len / 2;
	//uint8_t dmax = window_size; 
	
	// position-dependent d_{i, w}, where w is the window index
	uint8_t *d_iw = &d->read_dmax[read->id * d->n_hd_windows];
	std::fill_n(d_iw, d->n_hd_windows, dmax);

	const int nbhd_max_size[2] = {256, 32768};

	int widx = 0;
	for (int ts = tmin; ts < tmax; ts += window_size, ++widx) {
		int tw = ts + window_size < tmax ? ts + window_size : tmax;

		// if all quality scores within current window has quality score 2.
		bool q2_flag = std::all_of(&conv_q[ts+1], &conv_q[tw+kmer_len],
				[](char q) { return q == 2; });
		//bool low_cov_flag = std::all_of();

		int window_nbhd_size;
		bool rebuild = false, dmax_adjusted = false;
		bool final_k = (ts + kmer_len + 1) >= read->length;
		std::vector<uint8_t> mth_nsizes;
		mth_nsizes.reserve(window_size);
		do {
			window_nbhd_size = 0;
			// build neighborhood in the next window
			for (int t = ts; t < tw; ++t) {
				int tnext_kp = t + kmer_len;

				obs_next_base = kmer_effective_base(rseq[t + kmer_len]);
				obs_kid = bitarray_get(observed_kmers, t << 1, kmer_len << 1);
				kbits_t obs_next_kid = bitarray_get(observed_kmers, 
						(t + 1) << 1, kmer_len << 1);

				kmer_t *obs_pk = (kmer_t *) &(dict_find(kdict,
							kbits_cast_to_ptr(obs_kid))->value);

				kbits_t _new_nf = ((obs_next_base >> 2) << 1) | 
					(obs_next_base >> 2);
				next_obs_n_flag = (next_obs_n_flag >> 2) | 
					(_new_nf << shift);

				state_nbhd_t *curr_nbhd = &T_nbhds[t],
							 *next_nbhd = &T_nbhds[t+1];

				int t_nbhd_size = curr_nbhd->size;

				kbits_t *t_states_sorted_sfx = T_nbhds[t].states_sorted_sfx;
				double *t_states_alphas = T_nbhds[t].alphas;
				double **t_kmer_trans_p = T_nbhds[t].kmer_trans_prob;

				memset(actg_trans_count, 0, sizeof(int) << 2);

				/* one linear scan to determine # of distinct suffix, 
				 * compute downstream actg nucleotide counts etc. */
				int n_distinct_suffixes = hmm_count_distinct_suffixes(curr_nbhd, 
						actg_trans_count, dmax, obs_kid, obs_n_flag,
						obs_next_base, d->transition_p);

				curr_nbhd->n_uniq_sfx = n_distinct_suffixes;

				/* number of kmers in next position */
				int tnext_nbhd_size = actg_trans_count[0];
				/* cumulative transition counts */
				cum_actg_trans_count[0] = actg_trans_count[0];
				for (int i = 1; i < 4; i++) {
					tnext_nbhd_size += actg_trans_count[i];
					cum_actg_trans_count[i] = cum_actg_trans_count[i-1] +
						actg_trans_count[i];
				}

				window_nbhd_size += tnext_nbhd_size;
				if (tnext_nbhd_size == 0) {
					uint64_t t0_ham_dist = T_nbhds[0].ba_hamming_dist[0];
					mempool_destroy(mp);

					/* given that the current read runs out of pathway, it may
					 * signify that the first neighborhood is incorrectly
					 * constructed. 
					 * if the first neighborhood contains kmers with mismatches
					 * to the observed kmer, we will rebuild the first 
					 * neighborhood with the slow method. */
					if (t0_ham_dist > 0 && !skip_rebuild) {
						int old_nbhd_size = d->preconstructed_nbhds[read->id].n;
						kmer_t **old_nbhd = d->preconstructed_nbhds[read->id].nbhd;
#if DEBUG
						fprintf(stderr, 
								"[DEBUG] Rebuilding first neighborhood for read [%d] %s\n",
								read->id, read->identifier);

						fprintf(stderr, "[DEBUG] Old nbhd: ");
						for (int i = 0; i < d->preconstructed_nbhds[read->id].n; ++i) {
							fprintf(stderr, "%lu ", kbits_id_of(old_nbhd[i]));
						}
						fprintf(stderr, "\n");
#endif

						hmm_build_1st_nbhd_slow(d, read, NULL);

						// if the two neighborhoods are different, restart to 
						// determine dmax
						if (d->preconstructed_nbhds[read->id].n != old_nbhd_size) {
							return hmm_adaptive_dmax(d, read, fdata);
						}
						else {
							if (!std::equal(old_nbhd, old_nbhd + old_nbhd_size,
										d->preconstructed_nbhds[read->id].nbhd))
								return hmm_adaptive_dmax(d, read, fdata);
						}
#if DEBUG
						fprintf(stderr, "[DEBUG] New nbhd: ");
						for (int i = 0; i < d->preconstructed_nbhds[read->id].n; ++i) {
							fprintf(stderr, "%lu ", kbits_id_of(new_nbhd[i]));
						}
						fprintf(stderr, "\n");
#endif

					}

					return;
				}

				/* allocate memory given tnext_nbhd_size */
				char *pmem = (char *) mempool_nalloc(mp,
						nbhd_alloc_size(tnext_nbhd_size), 16);

				next_nbhd->size = tnext_nbhd_size;
				hmm_setup_nbhd_ptrs(next_nbhd, tnext_nbhd_size, pmem);

				int kmer_idx = 0, index_pfx = 0;
				register ksubstr_t mismatch_flag = ~(1UL << obs_next_base) |
					~((obs_next_base >> 2) - 1UL);

				bool incr_dmax = false;
				for (int i = 0; i < n_distinct_suffixes; i++) {
					kbits_t _repr_kid = t_states_sorted_sfx[kmer_idx];
					kbits_t sfx_hd = bitarray_get_pwr2(curr_nbhd->ba_hamming_dist, 
							i, 3);
					kbits_t common_sfx = kbits_suffix(_repr_kid);

					/* v 4 bits                   v 4 bits (higher)
					 *
					 * A  -                       - A
					 * C  |  -- k-1 substring --  | C
					 * T  |                       | T
					 * G  -                       - G
					 */
					kbits_t _tf = bitarray_get_pwr2(curr_nbhd->ba_pfx_sfx_flags,
							i, 3);
					kbits_t ups_trans_packed = trans_flag2packed(_tf & 15);
					/* which kmers/nucleotides current suffix transition into */
					kbits_t dns_trans_packed = trans_flag2packed(_tf >> 4);

					int n_common_sfx_kmers = ups_trans_packed & 7; 

					int j = dns_trans_packed & 7;
					for (dns_trans_packed >>= 3; j > 0; --j, dns_trans_packed >>= 2) {
						kbits_t dns_base = dns_trans_packed & 3;
						kbits_t tnext_kid = common_sfx | (dns_base << shift);

						/* compute the position sorted by suffix */
						int index_sfx = cum_actg_trans_count[dns_base] - 
							actg_trans_count[dns_base];
						actg_trans_count[dns_base]--;

						kmer_t *tnext_kmer = (kmer_t *) &(dict_find(kdict, 
								kbits_cast_to_ptr(tnext_kid))->value);
						next_nbhd->kmer_ptrs[index_sfx] = tnext_kmer;

						/* set up transition flag */
						bitarray_set_pwr2(next_nbhd->ba_kmer_trans_flag,
								index_sfx, tnext_kmer->trans_flag, 2);

						/* set up Hamming distance */
						kbits_t _hd = sfx_hd + ((mismatch_flag >> dns_base) & 1);
						bitarray_set_pwr2(next_nbhd->ba_hamming_dist, 
								index_sfx, _hd, 3);

						// if a certain kmer's Hamming distance reaches dmax,
						// and it has a unique transition, and the neighborhood
						// is still "sparse", increase dmax
						/*
						if (_hd == dmax && tnext_kmer->n_trans == 1 && 
								curr_nbhd->size <= PMR_SPARSE_NBHD_SIZE) {
							kbits_t _b = (trans_flag2packed(tnext_kmer->trans_flag) >> 3) & 3;
							int _np = t + kmer_len + 1;
							// check if an error will be emitted
							if (_np < read->length && 
									_b != kmer_effective_base(rseq[_np])) {
								incr_dmax = true;
							}
						}
						*/

						next_nbhd->suffixes_order[index_pfx] = index_sfx;
						next_nbhd->states_sorted_sfx[index_sfx] = tnext_kid;
						index_pfx++;
					}

					kmer_idx += n_common_sfx_kmers;
				}

				// if (incr_dmax) ++dmax;

				int _max_size = nbhd_max_size[q2_flag];
				// check if the neighborhood sizes are above the upper limit,
				// with two exceptions: the neighborhood size is allowed to
				// grow if the entire window has low quality score (Q=2).
				if (!rebuild && next_nbhd->size > _max_size) {

					size_t _n = BITS_TO_U64(next_nbhd->size << 3) << 3;
					uint8_t *hdist = new uint8_t[_n];
					memcpy(hdist, next_nbhd->ba_hamming_dist, _n);

					std::nth_element(hdist, hdist + _max_size - 1, 
							hdist + next_nbhd->size);

					mth_nsizes.push_back(hdist[_max_size - 1]);

					delete[] hdist;
				}

				obs_base = obs_next_base;
				obs_n_flag = next_obs_n_flag;
			}

			if (rebuild) break;

			if (mth_nsizes.size() > 0) {
				auto min_hd = std::min_element(mth_nsizes.begin(), mth_nsizes.end());
				// lower dmax so that the neighborhood sizes would not exceed
				// the upper limit, m, and then rebuild the neighborhoods.
				dmax = *min_hd;
				rebuild = true;
				dmax_adjusted = true;
			}
		} while (rebuild);

		if (dmax_adjusted) {
			uint8_t _d = dmax > 1 ? dmax - 1 : 1;
			d_iw[widx] = _d;

			// retrospectively apply a smaller dmax (for earlier windows)
			// so that later windows can use slightly larger dmax, coinciding
			// with the fact that error rate ticks up at 3' end.
			if (_d > 1) --_d;
			for (int w = widx - 1; w > 0; --w) {
				if (d_iw[w] > _d) d_iw[w] = _d; 
			}
		}

		// being opportunistic, and increase dmax by 1
		if (dmax < (kmer_len >> 1)) ++dmax;
	}


destroy:
	mempool_destroy(mp);
}
