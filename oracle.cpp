/* Implements a k-local oracle for error-correction */

#include <unordered_map>
#include <vector>
#include <utility>
#include <algorithm>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <omp.h>
#include <numeric>

#include <signal.h>

extern "C" {
#include "premier.h"
#include "kmer.h"
#include "numeric.h"
#include "trans.h"
#include "em.h"
#include "mstep.h"
#include "bitarray.h"
#include "read.h"
#include "iodata.h"
#include "atomic.h"
}

#include "nbhd.h"

using namespace std;
using GroundTruth = unordered_map<string, vector<pair<int, int>>>;

inline uint64_t kbits_mutation_flag(uint64_t p, uint64_t q) {
	uint64_t xor_pq = p ^ q;
	return xor_pq | ((xor_pq & KBITS_LO_MASK) << 1) | ((xor_pq & KBITS_HI_MASK) >> 1);
}

struct OraclePayload {
	GroundTruth &ground_truth;
	int *n_errors_corrected;
	int *n_errors_missed;
	bool tally_incidence;

	OraclePayload(GroundTruth &gt, int *tp, int *fn, bool tally = false) : 
		ground_truth(gt), n_errors_corrected(tp), n_errors_missed(fn), tally_incidence(tally) {};
};

int load_ground_truths(const char *f, GroundTruth &gnd_truth)
{
	int n_total_errors = 0;

	ifstream efs(f);
	string entry;
	
	while (getline(efs, entry)) {
		// skip comments/annotations
		if (entry[0] == '#') continue;

		size_t sp = entry.find(' ');
		// read identifier
		size_t sp2 = entry.find(' ', sp + 1);

		string error_tag = entry.substr(0, sp);
		string error_pos = entry.substr(sp + 1, sp2 - sp - 1);

		auto hit = gnd_truth.find(error_tag);
		if (hit == gnd_truth.end()) {
			vector<pair<int, int>> _v;
			_v.reserve(8);
			_v.push_back(make_pair(stoi(error_pos), (int) entry[sp2 + 1]));
			gnd_truth.emplace(make_pair(error_tag, _v));
			++n_total_errors;
		}
		else {
			hit->second.push_back(make_pair(stoi(error_pos), (int) entry[sp2 + 1]));
			++n_total_errors;
		}
	}

	return n_total_errors;
}

void mark_genomic_kmers(data *d, read_t *read, void *fdata)
{
	const int kmer_len = d->opt->kmer_length;
	const int shift = (kmer_len - 1) << 1;
	const int tmin = 0;
	const int tmax = read->length - kmer_len;

	dict *kdict = d->kmer_dict;
	const double _uniq_tp = 0.0;
	OraclePayload *p_payload = static_cast<OraclePayload *>(fdata);
	GroundTruth &gnd_truth = p_payload->ground_truth;
	const bool tally = p_payload->tally_incidence;

	auto hit = gnd_truth.find(string(read->identifier));
	char *corr_rseq = read->sequence;	/* why assignment? */
	if (hit != gnd_truth.end()) {
		auto &_v = hit->second;
		corr_rseq = new char[read->length + 1];
		corr_rseq[read->length] = '\0';
		memcpy(corr_rseq, read->sequence, read->length);

		for (auto it = _v.begin(); it != _v.end(); ++it) {
			corr_rseq[it->first - 1] = it->second;
		}
	}

	// mark all kmers as genomic kmer and:
	// - count transitions (tally) or
	// - set transition flag (!tally)
	// between k-spectrum kmers
	kbits_t obs_kid = kmer_seq_to_id(corr_rseq, kmer_len);
	dictentry *de = dict_find(kdict, kbits_cast_to_ptr(obs_kid));
	kmer_t *pk = NULL;
	if (de != NULL)
		pk = (kmer_t *) &de->value;
	kbits_t obs_next_base; // = kmer_effective_base(corr_rseq[tmin+kmer_len]);
	double *pk_exp_trans = (de == NULL) ? NULL : 
		(d->exp_trans_p + (pk->unique_trans ?  pk->index : d->nonuniq_kmers[pk->index].poff));

	for (int t = tmin; t <= tmax; ++t) {
		if (de != NULL)
			pk->has_multi_copies = 1;	// use this to mark true kmer

		if (t < tmax) {
			obs_next_base = kmer_effective_base(corr_rseq[t+kmer_len]);

			obs_kid = (obs_kid >> 2) | (obs_next_base << shift);
			dictentry *de_new = dict_find(kdict, kbits_cast_to_ptr(obs_kid));
			if (de != NULL && de_new != NULL) {
				// update transition flag/count
				if (tally && ((pk->trans_flag & (1 << obs_next_base)) > 0)) {
					int bidx = _mm_popcnt_u64(pk->trans_flag & ((1 << obs_next_base) - 1));
					atomic_add_dbl(&pk_exp_trans[bidx], 1.0);
				}
				else if (!tally)
					pk->alt_trans_flag |= (1 << obs_next_base);
			}
			/* count missing transitions */
			else if (tally && de != NULL && de_new == NULL)	/* temporary */
				pk->alt_trans_flag |= (1 << obs_next_base);

			de = de_new;
			if (de != NULL) {
				pk = (kmer_t *) &de_new->value; 
				pk_exp_trans = d->exp_trans_p + (pk->unique_trans ? 
					pk->index : d->nonuniq_kmers[pk->index].poff);
			}
		}

	}

	if (hit != gnd_truth.end()) delete [] corr_rseq;
}

void oracle_corrector(data *d, read_t *read, void *fdata)
{
	int actg_trans_count[4];
	int cum_actg_trans_count[4];

	const int kmer_len = d->opt->kmer_length;
	const int dmax = kmer_len >> 1;
	const int num_kmers = read->length - kmer_len + 1;
	const int shift = (kmer_len - 1) << 1;
	const kbits_t klen_flag = (1UL << (kmer_len << 1)) - 1;
	const int bwidth = d->opt->qscore_bin_width;

	const int tmin = 0;
	int tmax = read->length - kmer_len;

	dict *kdict = d->kmer_dict;
	double _uniq_tp = 0.0;

	OraclePayload *p_payload = static_cast<OraclePayload *>(fdata);
	GroundTruth &gnd_truth = p_payload->ground_truth;
	int *n_corred = &p_payload->n_errors_corrected[omp_get_thread_num()];
	int *n_err = &p_payload->n_errors_missed[omp_get_thread_num()];
	auto hit = gnd_truth.find(string(read->identifier));

	kbits_t observed_kmers[BITS_TO_U64(read->length << 1)];
	kbits_t error_flags[BITS_TO_U64(read->length << 1)];

	read_seq_to_numeric(read->sequence, observed_kmers, read->length, kmer_len);
	fill_n(error_flags, BITS_TO_U64(read->length << 1),  0UL);

	/* read has corrections */
	kbits_t obs_kid = bitarray_get(observed_kmers, tmin << 1, kmer_len << 1);
	kbits_t corrected_kid = obs_kid;
	int hd_1st_kmer = 0;
	char *corr_rseq;
	if (hit != gnd_truth.end()) {
		corr_rseq = new char[read->length + 1];
		corr_rseq[read->length] = '\0';
		memcpy(corr_rseq, read->sequence, read->length);
		for (auto it = hit->second.begin(); it != hit->second.end(); ++it) {
			int p = it->first - 1;
			corr_rseq[p] = it->second;
			/* in first kmer */
			if (p < kmer_len) {
				corrected_kid &= ~(3UL << (p << 1));
				corrected_kid |= static_cast<uint64_t>(base_to_bits(it->second)) << (p << 1); 
				++hd_1st_kmer;
			}

			bitarray_set_pwr2(error_flags, p, 3, 1);
		}
	}

	/* read has no errors or in first kmer and the correct first kmer is not
	 * in k-spectrum: give up correction on this read */
	if (hit == gnd_truth.end() ||
			(corrected_kid != obs_kid &&
			dict_find(kdict, kbits_cast_to_ptr(corrected_kid)) == NULL)) {
		FILE *stream = d->errors_fp ? d->errors_fp : stdout;
		fprintf(stream, "@%s\n%.*s\n+\n%.*s\n",
				read->identifier, read->length, read->sequence,
				read->length, read->qscore);
		*n_err += hit == gnd_truth.end() ? 0 : hit->second.size();
		/*
		if (hit != gnd_truth.end() && hit->second.size())
			INFO_TEE(stderr, "%s 0 %d %d\n", read->identifier, 
				hit == gnd_truth.end() ? 0 : hit->second.size(), 
				hit == gnd_truth.end() ? 0 : hit->second.size());
		*/
		return;
	}

	mempool_t *mp = mempool_create(MEMPOOL_DEFAULT_SIZE);

	char *rseq = read->sequence; 
	char *conv_q = convert_quality_scores(read->qscore, read->length, 
			d->opt->qual_score_offset, bwidth, mp);

	state_nbhd_t *T_nbhds = (state_nbhd_t *) mempool_nalloc(mp,
			sizeof(*T_nbhds) * num_kmers, 16);

	/* flag for N (not determined) base, first kmer only */
	kbits_t obs_n_flag = kmer_n_base_flag(rseq, kmer_len);
	kbits_t next_obs_n_flag = kmer_n_base_flag(rseq + 1, kmer_len);
	kbits_t obs_base = kmer_effective_base(rseq[tmin+kmer_len-1]);
	kbits_t obs_next_base = kmer_effective_base(rseq[tmin+kmer_len]);

	double tnext_max_delta = NAN;
	int tnext_argmax_delta = -1;

	state_nbhd_t *n1 = &T_nbhds[tmin];
	char *pmem = (char *) mempool_nalloc(mp, nbhd_alloc_size(1), 16);
	hmm_setup_nbhd_ptrs(n1, 1, pmem);

	kmer_t *corr_pk = kptr_locate(d, corrected_kid);
	n1->size = 1;
	n1->states_sorted_sfx[0] = corrected_kid;
	n1->kmer_ptrs[0] = corr_pk; 
	n1->kmer_trans_prob[0] = (corr_pk->unique_trans ? 
			&_uniq_tp : d->transition_p + 
			(d->nonuniq_kmers[corr_pk->index].poff));

	bitarray_set_pwr2(n1->ba_hamming_dist, 0, hd_1st_kmer, 3);
	bitarray_set_pwr2(n1->ba_kmer_trans_flag, 0, corr_pk->trans_flag, 2);

	int t = tmin;
	for (; t < tmax; t++) {
		kbits_t tnext_err_flag = bitarray_get(error_flags, (t+1) << 1, kmer_len << 1);
		kbits_t obs_next_kid = bitarray_get(observed_kmers, 
				(t + 1) << 1, kmer_len << 1);
		obs_kid = bitarray_get(observed_kmers, t << 1, kmer_len << 1);
		obs_next_base = kmer_effective_base(rseq[t + kmer_len]);

		kmer_t *obs_pk = (kmer_t *) &(dict_find(kdict,
					kbits_cast_to_ptr(obs_kid))->value);

		state_nbhd_t *curr_nbhd = &T_nbhds[t], *next_nbhd = &T_nbhds[t+1];
		int t_nbhd_size = curr_nbhd->size;

		fill_n(actg_trans_count, 4, 0);

		/* one linear scan to determine # of distinct suffix, 
		 * compute downstream actg nucleotide counts etc. */
		int n_distinct_suffixes = hmm_count_distinct_suffixes(curr_nbhd, 
				actg_trans_count, dmax, obs_kid, obs_n_flag,
				obs_next_base, d->transition_p);

		curr_nbhd->n_uniq_sfx = n_distinct_suffixes;

		/* actg_trans_count: number of suffices that can transition to A, C, G, T */
		partial_sum(actg_trans_count, actg_trans_count + 4,
				cum_actg_trans_count);
		int tnext_nbhd_size = cum_actg_trans_count[3];

		/* no valid transitions */
		if (tnext_nbhd_size == 0) {
			break;
		}

		/* allocate memory given tnext_nbhd_size */
		char *pmem = (char *) mempool_nalloc(mp,
				nbhd_alloc_size(tnext_nbhd_size), 16);

		next_nbhd->size = tnext_nbhd_size;
		hmm_setup_nbhd_ptrs(next_nbhd, tnext_nbhd_size, pmem);

		tnext_max_delta = NAN;
		tnext_argmax_delta = -1;

		int kmer_idx = 0, index_pfx = 0;
		register ksubstr_t mismatch_flag = ~(1UL << obs_next_base) |
			~((obs_next_base >> 2) - 1UL);
		for (int i = 0; i < n_distinct_suffixes; i++) {
			kbits_t _repr_kid = curr_nbhd->states_sorted_sfx[kmer_idx];
			kbits_t sfx_hd = bitarray_get_pwr2(curr_nbhd->ba_hamming_dist, 
					i, 3);
			kbits_t common_sfx = kbits_suffix(_repr_kid);

			/* nucleotides upstream and downstream of current (k-1)suffix */
			kbits_t _tf = bitarray_get_pwr2(curr_nbhd->ba_pfx_sfx_flags, i, 3);
			kbits_t ups_trans_packed = trans_flag2packed(_tf & 15);
			kbits_t dns_trans_packed = trans_flag2packed(_tf >> 4);

			/* number upstream kmers with this suffix */
			int n_common_sfx_kmers = ups_trans_packed & 7; 

			/* number transitions out of this suffix */
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

				/* next kmer is legit if it differs from observed only at errors
				 * (perfect oracle) and is true kmer */
				kbits_t mut_flag = kbits_mutation_flag(obs_next_kid, tnext_kid);
				bool is_legit_kmer = (mut_flag == tnext_err_flag) && 
					tnext_kmer->has_multi_copies;

				double delta;
				int argmax_delta = hmm_compute_delta(ups_trans_packed, 
						kmer_idx, curr_nbhd->kmer_ptrs, curr_nbhd->alphas,
						curr_nbhd->kmer_trans_prob, 0.0, dns_base, &delta);
				update_max_quantity(delta, index_sfx,
						tnext_max_delta, tnext_argmax_delta);

				/* set up transition flag, if the kmer in the neighborhood has
				 * mutations other than the ground truth errors, remove it from
				 * the neighborhood as it does not lead to a correct pathway
				 * for error correction */
				bitarray_set_pwr2(next_nbhd->ba_kmer_trans_flag,
						index_sfx, is_legit_kmer ? tnext_kmer->trans_flag : 0, 2);

				/* set up Hamming distance */
				kbits_t _hd = sfx_hd + ((mismatch_flag >> dns_base) & 1);
				bitarray_set_pwr2(next_nbhd->ba_hamming_dist, 
						index_sfx, _hd, 3);

				next_nbhd->suffixes_order[index_pfx] = index_sfx;
				next_nbhd->states_sorted_sfx[index_sfx] = tnext_kid;
				//next_nbhd->betas[index_sfx] = is_legit_kmer;
				next_nbhd->betas[index_sfx] = (double) argmax_delta;
				next_nbhd->alphas[index_sfx] = is_legit_kmer ? delta : NAN;

				next_nbhd->kmer_trans_prob[index_sfx] = (
						tnext_kmer->unique_trans ?  &_uniq_tp : 
						d->transition_p + (d->nonuniq_kmers[tnext_kmer->index].poff));

				index_pfx++;
			}

			kmer_idx += n_common_sfx_kmers;
		}
	}

	if (t == tmax && tnext_argmax_delta >= 0 && !IS_ZERO(tnext_max_delta)) {
		state_nbhd_t *nbhd = &T_nbhds[tmax];
		int n_corr = kmer_hamming_dist(
				nbhd->states_sorted_sfx[tnext_argmax_delta],
				(obs_kid >> 2) | (obs_next_base << shift));
		int n_fpn = kmer_hamming_dist(	/* original error count */
				(obs_kid >> 2) | (obs_next_base << shift),
				kmer_seq_to_id(&corr_rseq[tmax], kmer_len));

		char *corrected_rseq = new char[read->length + 1];
		memcpy(corrected_rseq, rseq, read->length);
		read->sequence = corrected_rseq;

		hmm_correct_errors(d, read, tmax,
				nbhd->states_sorted_sfx[tnext_argmax_delta],
				(obs_kid >> 2) | (obs_next_base << shift),
				obs_n_flag,
				kmer_len);

		int optimal_st_idx = (int) nbhd->betas[tnext_argmax_delta];

		for (int t = tmax - 1; t >= tmin; t--) {
			nbhd = &T_nbhds[t];

			hmm_correct_one_error(d, read, t,
					kbits_first_base(nbhd->states_sorted_sfx[optimal_st_idx]),
					kmer_effective_base(rseq[t]));

			n_corr += (kbits_first_base(nbhd->states_sorted_sfx[optimal_st_idx])
				!= kmer_effective_base(rseq[t]));
			n_fpn += (kmer_effective_base(rseq[t])
				!= kmer_effective_base(corr_rseq[t]));

			optimal_st_idx = (int) nbhd->betas[optimal_st_idx];
		}

		// if (n_corr < hit->second.size()) raise(SIGINT);
		*n_corred += n_corr;
		*n_err += n_fpn;

		FILE *stream = d->errors_fp ? d->errors_fp : stdout;
		fprintf(stream, "@%s\n%.*s\n+\n%.*s\n",
				read->identifier, read->length, corrected_rseq,
				read->length, read->qscore);

//		if (n_err) INFO_TEE(stderr, "%s %d %d %d\n", read->identifier, n_corr, n_err, hit->second.size());

		delete[] corrected_rseq;
	}
	/* ran out of pathway: output uncorrected read */
	else {
		FILE *stream = d->errors_fp ? d->errors_fp : stdout;
		fprintf(stream, "@%s\n%.*s\n+\n%.*s\n",
				read->identifier, read->length, read->sequence,
				read->length, read->qscore);

		*n_err += hit->second.size();
		//INFO_TEE(stderr, "%s 0 %d %d\n", read->identifier, hit->second.size(), hit->second.size());
		//raise(SIGINT);

		goto destroy;

	}

	// perform error-correction.
#if 0
	// --- "unique candidacy" rule ---

	int n_corr = 0;
	auto &_v = hit->second;
	char *corr_rseq = new char[read->length + 1];
	memcpy(corr_rseq, rseq, read->length);

	for (auto it = _v.begin(); it != _v.end(); ++it) {
		int pos_err = it->first - 1;

		// can always correct the first kmer
		if (pos_err < kmer_len) {
			corr_rseq[pos_err] = it->second;	
			++n_corr;
			//++*n_corred;
			continue;
		}

		int pmin = pos_err - kmer_len + 1;
		int pmax = pos_err >= tmax ? tmax : pos_err;
		bool has_uniq = false;
		for (int p = pmin; p <= pmax; ++p) {
			state_nbhd_t &nbhd = T_nbhds[p];
			// all possible error-corrections
			int n_possible_corrs = (int) accumulate(nbhd.betas, 
					nbhd.betas + nbhd.size, 0);

			if (p > t) {
				has_uniq = false;
				break;
			}

			// unambiguous, unique correction
			if (n_possible_corrs == 1) {
				has_uniq = true;
			}
		}

		// heuristic: last kmer? 
		// if (pos_err > tmax) has_uniq = true;

		// so long as one kmer covering this error can uniquely resolve the
		// mismatch, the Oracle is capable of recovering the error.
		if (has_uniq) {
			corr_rseq[pos_err] = it->second;	
			++n_corr;
		}
	}

	*n_corred += n_corr;
	
	FILE *stream = d->errors_fp ? d->errors_fp : stdout;
	fprintf(stream, "@%s\n%.*s\n+\n%.*s\n",
			read->identifier, read->length, corr_rseq,
			read->length, read->qscore);

#endif

destroy:
	mempool_destroy(mp);
	return;
}

void oracle_update_kspectrum(data *d, OraclePayload &payload)
{
	double jnt_exp_trans[4];
	const int kmer_len = d->opt->kmer_length;
	const int shift = (kmer_len - 1) << 1;

	// reset alternative transition flags.
	INFO_TEE(d->log_fp, "Resetting alternative transition flags...\n");
	int64_t n_kmers = d->kmer_dict->sizemask + 1;
	for_each(d->kmer_dict->table, d->kmer_dict->table + n_kmers,
			[](dictentry &de) {
				if (de.value != NULL) {	
					kmer_t *pk = (kmer_t *) &(de.value);
					pk->alt_trans_flag = 0;
				}
			});

#if PMR_ORACLE_SUPPLEMENT_TRANS
	// identify "true genomic" kmers using the ground truth errors (set has_multi_copies);
	// reidentify transition flags using genomic transitions only (set alt_trans_flag).
	INFO_TEE(d->log_fp, "Identifying true genomic kmers and transition flags using ground truth...\n");
	io_iterate_reads_fastq(d, d->opt->n_threads, d->opt->is_tty, 
			mark_genomic_kmers, NULL, static_cast<void *>(&payload));

	int n_true_trans_removed = 0;

	INFO_TEE(d->log_fp, "Updating reverse complement transition flags...\n");
	// handles reciprocal transitions: A -> B ==> rc(B) -> rc(A).
	// simultaneously, replicate the "true genomic" kmer flag to the opposite
	// strand.
	for (int64_t i = 0; i < n_kmers; ++i) {
		dictentry *de = &d->kmer_dict->table[i];
		if (de->value == NULL) continue;

		kmer_t *pk = (kmer_t *) &(de->value);
		kbits_t kid = kbits_id_of(pk);

		if (pk->alt_trans_flag > 0) {
			kbits_t tpacked = trans_flag2packed(pk->alt_trans_flag);
			int n_genom_trans = tpacked & 7, b = 0;

			for (tpacked >>= 3; b < n_genom_trans; ++b, tpacked >>= 2) {
				kbits_t dns_base = tpacked & 3;
				kmer_t *dns_pk = kptr_locate(d, kmer_reverse_complement(
							(kid >> 2) | (dns_base << shift), kmer_len));

				// reciprocal flag
				dns_pk->alt_trans_flag |= (1U << complementary_base(kid & 3));
			}
		}

		kbits_t rev_kid = kmer_reverse_complement(kid, kmer_len);
		kmer_t *rev_pk = kptr_locate(d, rev_kid); 

		if (pk->has_multi_copies) 	// used to indicate true genome kmer
			rev_pk->has_multi_copies = 1;
	}

	// rearrange memory layout, as some unique kmer may become non-unique
	// kmers.
	INFO_TEE(d->log_fp, "Rearranging transition parameter memory layout...\n");
	int n_nonuniq_trans = 1;
	int n_nonuniq_kmers = 0;
	for (int64_t i = 0; i < n_kmers; ++i) {
		dictentry *de = &d->kmer_dict->table[i];
		if (de->value == NULL) continue;

		kmer_t *pk = (kmer_t *) &(de->value);
		int n = trans_flag2packed(pk->alt_trans_flag) & 7;
		if (n > 1) {
			n_nonuniq_trans += n;
			++n_nonuniq_kmers;
		}
	};

	free(d->nonuniq_kmers);
	d->nonuniq_kmers = (kmer_ext_t *)calloc(
		n_nonuniq_kmers, sizeof(*d->nonuniq_kmers));

	// reset kmer indices, offset by 1 as per data.c
	int uidx = n_nonuniq_trans, poff = 1, kidx = 0;
	for (int64_t i = 0; i < n_kmers; ++i) {
		dictentry *de = &d->kmer_dict->table[i];
		if (de->value == NULL) continue;

		kmer_t *pk = (kmer_t *) &(de->value);
		int n = trans_flag2packed(pk->alt_trans_flag) & 7;
		if (n == 1) {
			pk->unique_trans = 1;
			pk->index = uidx++;
		}
		else if (n > 1){
			pk->index = kidx++;
			d->nonuniq_kmers[pk->index].poff = poff;
			poff += n;
		}
		else {
			pk->unique_trans = 1;
			pk->index = 0;
		}

		for (int b = 0; b < 4; ++b) {
			uint32_t mask = (1U << b);
			if (((pk->alt_trans_flag & mask) > 0) &&
					((pk->trans_flag & mask) == 0)) {
				++n_true_trans_removed;
			}
		}

		pk->n_trans = n;
		pk->trans_flag = pk->alt_trans_flag;
	}

	INFO_TEE(d->log_fp, "True transitions missing from the transition flags: %d\n",
			n_true_trans_removed);

	free(d->exp_trans_p);
	free(d->transition_p);

	d->exp_trans_p = (double *) calloc(uidx, sizeof(double));
	d->transition_p = (double *) calloc(poff, sizeof(double));
	d->n_total_trans = uidx;
#else
	/* has_multi_copies flag used in oracle corrector */
/*
	for (int64_t i = 0; i < n_kmers; ++i) {
		dictentry *de = &d->kmer_dict->table[i];
		if (de->value == NULL) continue;

		kmer_t *pk = (kmer_t *) &(de->value);
		kbits_t kid = kbits_id_of(pk);
		kbits_t rev_kid = kmer_reverse_complement(kid, kmer_len);
		kmer_t *rev_pk = kptr_locate(d, rev_kid); 
		if (pk->has_multi_copies) 	// used to indicate true genome kmer
			rev_pk->has_multi_copies = 1;
	}
*/
	/* reset expected transitions to zero (preparing for tallying transition
	 * incidence.) */ 
	memset(d->exp_trans_p, 0, d->n_total_trans * sizeof(*d->exp_trans_p));

	/* only necessary if alt_trans_flag used above during search for missing transitions */
/*
	for_each(d->kmer_dict->table, d->kmer_dict->table + n_kmers,
			[](dictentry &de) {
				if (de.value != NULL) {	
					kmer_t *pk = (kmer_t *) &(de.value);
					pk->alt_trans_flag = 0;
				}
			});
*/

#endif

	INFO_TEE(d->log_fp, 
			"Tallying kmer transitions on the true de Bruijn graph.\n");

	payload.tally_incidence = true;

	/* count number of transitions observed in corrected data,
	 * and count missing transitions using alt_trans_flag to avoid double-counting */
	io_iterate_reads_fastq(d, d->opt->n_threads, d->opt->is_tty, 
			mark_genomic_kmers, NULL, static_cast<void *>(&payload));

#if !PMR_ORACLE_SUPPLEMENT_TRANS

	/* Many transitions will be dropped because they are not observed in
	 * the ground truth.  This code does the dropping & resetting. */

	INFO_TEE(d->log_fp, "Rearranging transition parameter memory layout...\n");
	/* count number of non-unique kmers */
	int n_nonuniq_trans = 1;
	int n_nonuniq_kmers = 0;
	int n_total_missing_transitions = 0;
	for (int64_t i = 0; i < n_kmers; ++i) {
		dictentry *de = &d->kmer_dict->table[i];
		if (de->value == NULL) continue;

		kmer_t *pk = (kmer_t *) &(de->value);
		double *exp_trans_p = d->exp_trans_p + (pk->unique_trans ? 
				pk->index : d->nonuniq_kmers[pk->index].poff);
		// combine coverage from both strands.
		mstep_reciprocal_exp_trans(d, pk, &jnt_exp_trans[0]);
		kbits_t _tfp = trans_flag2packed(pk->trans_flag & 15);
		pk->trans_flag = 0;
		int ntrans = 0, b;
		double covg = 0.0;
		for (_tfp >>= 3, b = 0; b < pk->n_trans; ++b, _tfp >>= 2) {
			kbits_t base = _tfp & 3;
			jnt_exp_trans[b] += exp_trans_p[b];
			covg += jnt_exp_trans[b];
			if (jnt_exp_trans[b] > 0) {
				++ntrans;
				pk->trans_flag |= (1UL << base);
			}
			if (pk->alt_trans_flag & (1UL << base))
				++n_total_missing_transitions;
		}

		if (ntrans > 1) {
			n_nonuniq_trans += ntrans;
			++n_nonuniq_kmers;
		}
		pk->n_trans = ntrans;
//		pk->trans_flag = pk->alt_trans_flag;
	}

	INFO_TEE(d->log_fp, "True transitions missing from the transition flags: %d\n",
			n_total_missing_transitions);

	/* rearrange memory */
	free(d->nonuniq_kmers);
	d->nonuniq_kmers = (kmer_ext_t *)calloc(
		n_nonuniq_kmers, sizeof(*d->nonuniq_kmers));

	int uidx = n_nonuniq_trans, poff = 1, kidx = 0;
	for (int64_t i = 0; i < n_kmers; ++i) {
		dictentry *de = &d->kmer_dict->table[i];
		if (de->value == NULL) continue;

		kmer_t *pk = (kmer_t *) &(de->value);

		if (pk->n_trans == 1) {
			pk->unique_trans = 1;
			pk->index = uidx++;
		}
		else if (pk->n_trans > 1){
			pk->index = kidx++;
			d->nonuniq_kmers[pk->index].poff = poff;
			poff += pk->n_trans;
		}
		else {
			pk->unique_trans = 1;
			pk->index = 0;
		}
	};

	free(d->exp_trans_p);
	free(d->transition_p);

	d->exp_trans_p = (double *) calloc(uidx, sizeof(double));
	d->transition_p = (double *) calloc(poff, sizeof(double));
	d->n_total_trans = uidx;

	/* wasteful, since we have already done this calculation in another memory space */
	INFO_TEE(d->log_fp, 
			"Retallying kmer transitions on the true de Bruijn graph.\n");
	io_iterate_reads_fastq(d, d->opt->n_threads, d->opt->is_tty, 
			mark_genomic_kmers, NULL, static_cast<void *>(&payload));

#endif

	INFO_TEE(d->log_fp, "Computing transition probabilities...\n");
	for (int64_t i = 0; i < n_kmers; ++i) {
		dictentry *de = &d->kmer_dict->table[i];
		if (de->value == NULL) continue;

		kmer_t *pk = (kmer_t *) &(de->value);
		// if (pk->n_trans <= 1) continue;

		double *exp_trans_p = d->exp_trans_p + (pk->unique_trans ? 
				pk->index : d->nonuniq_kmers[pk->index].poff);
		double *cur_tp = d->transition_p + (pk->unique_trans ?
				0 : d->nonuniq_kmers[pk->index].poff);

		// combine coverage from both strands.
		mstep_reciprocal_exp_trans(d, pk, &jnt_exp_trans[0]);
		double sum_exp_trans = 0.0;
//fprintf(stderr, "%s (%d):", kmer_id_to_seq(kbits_id_of(pk), kmer_len, NULL), pk->trans_flag);
		for (int b = 0; b < pk->n_trans; ++b) {
			jnt_exp_trans[b] += exp_trans_p[b];
			sum_exp_trans += jnt_exp_trans[b];
//fprintf(stderr, " %f", jnt_exp_trans[b]);
		}
//fprintf(stderr, "\n");

		if (pk->n_trans > 1) {
			// MLEs of the transition probabilities using transition counts from
			// the TRUE de Bruijn graph.
			for (int b = 0; b < pk->n_trans; ++b) {
				cur_tp[b] = LOG(jnt_exp_trans[b] / sum_exp_trans);
			}
		}
		else if (pk->n_trans == 1 && sum_exp_trans == 0.0) {
			/* Unique transition: since it DOES not have a transition
			probability parameter, we cannot set such probability to zero. 
			Instead, set the transition flag to zero. */
			pk->n_trans = 0;
			pk->trans_flag = 0;
			pk->index = 0;
		}
		// else: do nothing
	}
}

void _neighborhood_mapper(data *d, read_t * read, void *fdata)
{
	OraclePayload *p_payload = static_cast<OraclePayload *>(fdata);
	GroundTruth &gnd_truth = p_payload->ground_truth;

	const int kmer_len = d->opt->kmer_length;

	auto hit = gnd_truth.find(string(read->identifier));
	char *corr_rseq = read->sequence;
	if (hit != gnd_truth.end()) {
		auto &_v = hit->second;
		corr_rseq = new char[read->length + 1];
		corr_rseq[read->length] = '\0';
		memcpy(corr_rseq, read->sequence, read->length);

		/* wasted iterations on errors >k */
		for (auto it = _v.begin(); it != _v.end(); ++it) {
			corr_rseq[it->first - 1] = it->second;
		}
	}

	kbits_t obs_kid = kmer_seq_to_id(corr_rseq, kmer_len);
	dictentry *de = dict_find(d->kmer_dict, kbits_cast_to_ptr(obs_kid));
	if (de == NULL) {
		d->preconstructed_nbhds[read->id] = (kmer_nbhd_t) {0, NULL};
	}
	else {
		kmer_t **new_nbhd = (kmer_t **) mempool_alloc(d->mp, 1 * sizeof(kmer_t *));
		new_nbhd[0] = (kmer_t *) &(de->value);
		d->preconstructed_nbhds[read->id] = (kmer_nbhd_t) {1, new_nbhd};
	}

	if (hit != gnd_truth.end()) delete[] corr_rseq;
}

void oracle_build_first_neighborhoods(data *d, OraclePayload &payload)
{
	io_iterate_reads_fastq(d, 1, d->opt->is_tty, 
			_neighborhood_mapper, NULL, static_cast<void *>(&payload));
}

void oracle_init_em(data *d)
{
	GroundTruth gnd_truth_errors;
	int n_total_errors = load_ground_truths(d->opt->ground_truth_file, 
			gnd_truth_errors);

	OraclePayload payload(gnd_truth_errors, NULL, NULL);

	oracle_build_first_neighborhoods(d, payload);

	if (d->opt->adaptive_hamming_dist) {
		INFO_TEE(d->log_fp, "Determining read-specific max. Hamming distances...\n");
		io_iterate_reads_fastq(d, d->opt->n_threads, d->opt->is_tty, hmm_adaptive_dmax, NULL, (void *) 1);
	}
	else {
		/* set all dmax to be the default */
		memset(d->read_dmax, d->opt->max_hamming_dist, d->n_reads * d->n_hd_windows * sizeof(uint8_t));
	}

	oracle_update_kspectrum(d, payload);
}

void run_oracle(data *d)
{
	GroundTruth gnd_truth_errors;
	int *n_err_corred = new int[d->opt->n_threads];
	fill_n(n_err_corred, d->opt->n_threads, 0);
	int *n_err = new int[d->opt->n_threads];
	fill_n(n_err, d->opt->n_threads, 0);

	OraclePayload payload(gnd_truth_errors, n_err_corred, n_err);

	int n_total_errors = load_ground_truths(d->opt->ground_truth_file, 
			gnd_truth_errors);

	INFO_TEE(d->log_fp, "Total number of errors after load_ground_truths(): %d\n", n_total_errors);

	/* if PMR_ORACLE_SUPPLEMENT_TRANS, supplement transitions with true transitions
	 * (but not kmers) from ground truth, and recompute transition probabilities 
	 * from ground truth */
	oracle_update_kspectrum(d, payload);

	INFO_TEE(d->log_fp, "Total number of errors: %d\n", 
		accumulate(n_err_corred, n_err_corred + d->opt->n_threads, 0));

	fill_n(n_err_corred, d->opt->n_threads, 0);

	/* load d->preconstructed_nbhds with true kmer if in k-spectrum and
	 * NULL if not in k-spectrum */
	oracle_build_first_neighborhoods(d, payload);
	
	INFO_TEE(d->log_fp, "Performing error correction with Oracle.\n");
	io_iterate_reads_fastq(d, d->opt->n_threads,
			d->opt->is_tty, oracle_corrector, NULL,
			static_cast<void *>(&payload));

	int n_total_corred = accumulate(n_err_corred, 
			n_err_corred + d->opt->n_threads, 0);
	int n_total_fn = accumulate(n_err, n_err + d->opt->n_threads, 0);

	INFO_TEE(d->log_fp, 
			"Errors recovered by Oracle: %d / %d (gain: %4.3g, %d)\n",
			n_total_corred, n_total_errors,
			(double) n_total_corred / n_total_errors, n_total_fn);

	delete[] n_err_corred;
	delete[] n_err;
}
