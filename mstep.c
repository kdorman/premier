#include "mstep.h"
#include "trans.h"
#include "kmer.h"
#include "bitarray.h"
#include "em.h"
#include "sort.h"
#include "atomic.h"
#include "iodata.h"

#include <omp.h>
#include <assert.h>
#include <signal.h>
#include <stdlib.h>

#define MATHLIB_STANDALONE 1
#include <Rmath.h>

#define EPSILON 1e-6


//static int n_exhaustive = 0;
//static int n_monte_carlo = 0;

// #define PMR_CHECK_ASCENT 1

/*
static double mstep_estimate_trans_p(data *d, kmer_t *pk, 
		double thresh, mstep_data_t *md);
*/

static double mstep_compute_reverse_trans(data *d, kmer_t *rev_pk, double *trans_p,
		double mu, int kmer_len, tflag_t *tf, double *exp_counts, double *k_cexp_cnt);
static double mstep_extend_monte_carlo(data *d, kmer_t *pk);

static inline int random_transition(double *tp, int n, unsigned int *seed)
{
	double sum_tp = 0.0;
	double r = (double) rand_r(seed) / RAND_MAX;	/* TODO: rand_r() is in POSIX.1-2001 and marked obsolete in POSIX.1-2008, thus it triggers a warning under C99 */

	for (int b = 0; b < n; ++b) {
		sum_tp += EXP(tp[b]);
		if (r < sum_tp) return b;
	}

	return -1;
	//return n - 1;
}

static inline int random_state(double *tp, int n, double w, unsigned int *seed)
{
	//double sum_tp = 0.0;
	double r = (double) rand_r(seed) / RAND_MAX * w;

	for (int b = 0; b < n - 1; ++b) {
		//sum_tp += tp[b];
		//if (r < sum_tp) return b;
		if (r < tp[b]) return b;
	}

	//return -1;
	return n - 1;
}

double mstep_compute_reverse_trans(data *d, kmer_t *rev_pk, double *trans_p,
		double mu, int kmer_len, tflag_t *tf, double *exp_counts, double *k_cexp_cnt)
{
	int shift = (kmer_len - 1) * 2;
	kbits_t rev_kid = kbits_id_of(rev_pk);
	kbits_t k_trans_packed = trans_flag2packed(rev_pk->trans_flag);
	kbits_t n_trans = k_trans_packed & 7;
	int i = 0;

	/*
	if (n_trans == 1) {
		trans_p[0] = 1.0;
		tf->trans_flag = rev_pk->trans_flag;
		tf->n_nonzero_trans = 1;

		return 1.0;
	}
	*/

	double sum_tp = 0.0;
	double kcexp_rev = 0.0; 
	for (k_trans_packed >>= 3; i < n_trans; i++, k_trans_packed >>= 2) {
		kbits_t dns_base = k_trans_packed & 3;
		kbits_t fwd_dns_base = complementary_base(rev_kid & 3);

		kbits_t rc_fwd_kid = kbits_suffix(rev_kid) | (dns_base << shift);
		kbits_t fwd_kid = kmer_reverse_complement(rc_fwd_kid, kmer_len);

		kmer_t *fwd_pk = (kmer_t *) &(dict_find(d->kmer_dict, 
				kbits_cast_to_ptr(fwd_kid))->value);

		int par_off = kmer_initial_index(fwd_pk, d->uniq_kmer_init_off);
		double mu_dns = d->avg_initial_p[par_off];
		double _fwd_cnt = exp_counts[par_off];

		if ((fwd_pk->trans_flag & (1 << fwd_dns_base)) == 0) {
			trans_p[i] = 0.0;
		}
		else {
			double rtp = PRODUCT(mu_dns, -mu); 
			if (fwd_pk->n_trans > 1) {
				int fwd_nbidx = base_to_transidx(fwd_pk->trans_flag,
						(uint32_t) fwd_dns_base);
				double _tp = *(d->transition_p +
						d->nonuniq_kmers[fwd_pk->index].poff + fwd_nbidx);
				rtp += _tp;

				_fwd_cnt *= EXP(_tp);
			}

			kcexp_rev += _fwd_cnt;
			trans_p[i] = EXP(rtp);
			sum_tp += trans_p[i];
			register int _nz = (trans_p[i] > 0);
			tf->trans_flag |= (_nz << dns_base);
			tf->n_nonzero_trans += _nz;
		}
	}

	*k_cexp_cnt = kcexp_rev;
	return sum_tp;
}

uint32_t mstep_reciprocal_exp_trans(data *d, kmer_t *pk, double *recip_exp_trans)
{
	int kmer_len = d->opt->kmer_length;
	int shift = (kmer_len - 1) << 1;
	kbits_t k_trans_packed = trans_flag2packed(pk->trans_flag);
	kbits_t n_trans = k_trans_packed & 7;
	int i = 0;

	uint32_t ss = 0;
	kbits_t kid = kbits_id_of(pk);

	for (k_trans_packed >>= 3; i < n_trans; i++, k_trans_packed >>= 2) {
		kbits_t dns_base = k_trans_packed & 3;
		kbits_t rc_dns_base = complementary_base(kid & 3);

		kbits_t dns_kid = kbits_suffix(kid) | (dns_base << shift);
		kbits_t rc_dns_kid = kmer_reverse_complement(dns_kid, kmer_len);
		kmer_t *rc_pk = (kmer_t *) &(dict_find(d->kmer_dict, 
				kbits_cast_to_ptr(rc_dns_kid))->value);

		ss |= (rc_pk->single_stranded << i);

		if ((rc_pk->trans_flag & (1 << rc_dns_base)) == 0) {
			/* no such transition on reverse complement kmer */
			recip_exp_trans[i] = 0.0;
		}
		else {

			int rc_nbidx = _mm_popcnt_u64(rc_pk->trans_flag & 
					((1 << rc_dns_base) - 1));

			double *rc_exp_trans_p = d->exp_trans_p + (rc_pk->unique_trans ? 
					rc_pk->index : d->nonuniq_kmers[rc_pk->index].poff);

			recip_exp_trans[i] = rc_exp_trans_p[rc_nbidx];
		}
	}

	return ss;
}

double mstep_extend_kmer(data *d, kmer_t *pk, int recur)
{
	//static double _uniq_tp = 0.0;
	const int kmer_len = d->opt->kmer_length;
	const int shift = (kmer_len - 1) << 1;
	const int _uoff = d->uniq_kmer_init_off;

	const int R = 5;  // maximum recursion level.

	if (pk->alt_trans_flag == 0) {
		pk->dirty = 1;
		if (pk->n_trans > 0) {
			d->prob_kmer_extension[kmer_initial_index2(pk, _uoff)] = NAN;
		}
		return NAN;
	}

	/*
	int tid = omp_get_thread_num();
	printf("Thread %d, kmer index: %d, recur: %d\n",
			tid, (dictentry *) ((uintptr_t) pk - 0x8) - d->kmer_dict->table, recur);
	*/

	if (recur > R) {
		//int ext_len = 0;
		//double prop_extend = mstep_extend_monte_carlo(d, pk, &ext_len);
		double prop_extend;
		if (pk->dirty == 0) {
			prop_extend = mstep_extend_monte_carlo(d, pk);
			// if (prop_extend == 0.0) prop_extend = -ext_len;
			d->prob_kmer_extension[kmer_initial_index2(pk, _uoff)] = prop_extend;
			/*
			atomic_cas_dbl(&d->prob_kmer_extension[
					kmer_initial_index2(pk, d->uniq_kmer_init_off)],
					prop_extend);
			*/
			pk->dirty = 1;
		}
		else {
			prop_extend = d->prob_kmer_extension[
					kmer_initial_index2(pk, d->uniq_kmer_init_off)];
		}

		return prop_extend;
	}

	kbits_t kid = kbits_id_of(pk);

	double prop_extend = 0.0;
	kbits_t tpacked = trans_flag2packed(pk->trans_flag) >> 3;
	/*
	double *tp = pk->unique_trans ? &_uniq_tp : 
		d->transition_p + d->nonuniq_kmers[pk->index].poff;
	*/
	double *tp = d->exp_trans_p + (pk->unique_trans ? pk->index : 
		d->nonuniq_kmers[pk->index].poff);

	// this extends slightly beyond the maximum length (+R at most)
	
	//int extends = 1;
	for (int b = 0; b < pk->n_trans; ++b, tpacked >>= 2) {
		kbits_t base = tpacked & 3;
		kmer_t *pkn = kptr_locate(d, (kid >> 2) | (base << shift));
		
		double w = pkn->dirty ? 
			d->prob_kmer_extension[kmer_initial_index2(pkn, _uoff)] :
			mstep_extend_kmer(d, pkn, recur+1);

		//extends &= (w == 0.0);
		prop_extend += EXP(PRODUCT(w, tp[b]));
	}

	/*
	if (extends) {
		prop_extend = 0.0;
	}
	else {
		prop_extend = mstep_extend_monte_carlo(d, pk);
	}
	*/

	prop_extend = LOG(prop_extend);

	/*
	atomic_cas_dbl(&d->prob_kmer_extension[
			kmer_initial_index2(pk, d->uniq_kmer_init_off)],
			prop_extend);
	*/

	d->prob_kmer_extension[kmer_initial_index2(pk, _uoff)] = prop_extend;
	pk->dirty = 1;

	return prop_extend;
}

double mstep_extend_monte_carlo(data *d, kmer_t *pk)
{
	int actg_trans_count[4];
	int cum_actg_trans_count[4];

	const int kmer_len = d->opt->kmer_length;
	const int shift = (kmer_len - 1) << 1;
	//static double _uniq_tp = 0.0, _zero_tp = NAN;

	const int max_nbhd_size = 1024;
	const int max_len = kmer_len; //100;
	const int N = 1000;

	mempool_t *mp = mempool_create(MEMPOOL_DEFAULT_SIZE);

	state_nbhd_t *T_nbhds = (state_nbhd_t *) mempool_nalloc(mp,
			sizeof(*T_nbhds) * (max_len + 1), 16);

	state_nbhd_t *first_nbhd = &T_nbhds[0];
	first_nbhd->size = 1;
	char *pmem = (char *) mempool_nalloc(mp, nbhd_alloc_size(1), 16);
	hmm_setup_nbhd_ptrs(first_nbhd, 1, pmem);

	first_nbhd->states_sorted_sfx[0] = kbits_id_of(pk);
	first_nbhd->kmer_ptrs[0] = pk; 
	first_nbhd->alphas[0] = 0.0;
	/*
	if (pk->unique_trans) {
		first_nbhd->kmer_trans_prob[0] = (pk->alt_trans_flag == 0) ?
			 &_zero_tp: &_uniq_tp;
	}
	else {
		first_nbhd->kmer_trans_prob[0] = d->transition_p + d->nonuniq_kmers[pk->index].poff;
	}
	*/
	first_nbhd->kmer_trans_prob[0] = d->exp_trans_p + 
		(pk->unique_trans ? pk->index : d->nonuniq_kmers[pk->index].poff);

	bitarray_set_pwr2(first_nbhd->ba_hamming_dist, 0, 0x7f, 3);
	bitarray_set_pwr2(first_nbhd->ba_kmer_trans_flag, 0, pk->trans_flag, 2);

	double tnext_max_alpha = NAN;
	int tnext_argmax_alpha = 0;

	/* initial state distribution */
	int t = 0;
	for (; t < max_len; ++t) {
		state_nbhd_t *curr_nbhd = &T_nbhds[t],
					 *next_nbhd = &T_nbhds[t+1];

		tnext_max_alpha = NAN;
		tnext_argmax_alpha = 0;

		kbits_t *t_states_sorted_sfx = T_nbhds[t].states_sorted_sfx;

		memset(actg_trans_count, 0, sizeof(int) << 2);

		int n_distinct_suffixes = 0;
		int tnext_nbhd_size = 0;

		/* one linear scan to determine # of distinct suffix, 
		 * compute downstream actg nucleotide counts etc. */
		n_distinct_suffixes = hmm_count_distinct_suffixes(curr_nbhd, 
				actg_trans_count, UINT64_MAX, 0UL, 0UL, 0UL, d->transition_p);

		curr_nbhd->n_uniq_sfx = n_distinct_suffixes;

		/* number of kmers in next position */
		tnext_nbhd_size = actg_trans_count[0];
		/* cumulative transition counts */
		cum_actg_trans_count[0] = actg_trans_count[0];
		for (int i = 1; i < 4; i++) {
			tnext_nbhd_size += actg_trans_count[i];
			cum_actg_trans_count[i] = cum_actg_trans_count[i-1] +
				actg_trans_count[i];
		}

		/* allocate memory given tnext_nbhd_size */
		char *pmem = (char *) mempool_nalloc(mp,
				nbhd_alloc_size(tnext_nbhd_size), 16);

		next_nbhd->size = tnext_nbhd_size;
		hmm_setup_nbhd_ptrs(next_nbhd, tnext_nbhd_size, pmem);

		if (tnext_nbhd_size == 0) {
			// use curr_nbhd as the contexts
			break;
		}

		int kmer_idx = 0, index_pfx = 0;
		for (int i = 0; i < n_distinct_suffixes; i++) {
			kbits_t _repr_kid = t_states_sorted_sfx[kmer_idx];
			kbits_t common_sfx = kbits_suffix(_repr_kid);
			kbits_t _tf = bitarray_get_pwr2(curr_nbhd->ba_pfx_sfx_flags, i, 3);
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

				kmer_t *tnext_kmer = (kmer_t *) &(dict_find(d->kmer_dict, 
						kbits_cast_to_ptr(tnext_kid))->value);
				next_nbhd->kmer_ptrs[index_sfx] = tnext_kmer;

				double alpha = hmm_compute_alpha(
						ups_trans_packed, 
						kmer_idx,
						curr_nbhd->kmer_ptrs,
						curr_nbhd->alphas,
						curr_nbhd->kmer_trans_prob,
						0.0, dns_base);

				update_max_quantity(alpha, index_sfx,
						tnext_max_alpha, tnext_argmax_alpha);

				/* set up transition flag */
				if (!IS_ZERO(alpha)) {
					bitarray_set_pwr2(next_nbhd->ba_kmer_trans_flag,
							index_sfx, tnext_kmer->trans_flag, 2);
				}
				bitarray_set_pwr2(next_nbhd->ba_hamming_dist, index_sfx, 
						0x7f, 3);

				next_nbhd->suffixes_order[index_pfx] = index_sfx;
				next_nbhd->states_sorted_sfx[index_sfx] = tnext_kid;
				next_nbhd->alphas[index_sfx] = alpha;
				/*
				if (tnext_kmer->unique_trans) {
					next_nbhd->kmer_trans_prob[index_sfx] = 
						(tnext_kmer->alt_trans_flag == 0) ? 
						&_zero_tp : &_uniq_tp;
				}
				else {
					next_nbhd->kmer_trans_prob[index_sfx] = d->transition_p + d->nonuniq_kmers[tnext_kmer->index].poff;
				}
				*/
				next_nbhd->kmer_trans_prob[index_sfx] = 
					d->exp_trans_p + (tnext_kmer->unique_trans ? 
							tnext_kmer->index : d->nonuniq_kmers[tnext_kmer->index].poff);

				index_pfx++;
			}

			kmer_idx += n_common_sfx_kmers;
		}

		// sparsify the neighborhood.
		for (int i = 0; i < tnext_nbhd_size; ++i) {
			// 20 order of magnitude.
			if (tnext_max_alpha - next_nbhd->alphas[i] > PMR_NBHD_THIN_THRESHOLD) {
				bitarray_set_pwr2(next_nbhd->ba_kmer_trans_flag,
						i, 0, 2);
			}
		}

		if (tnext_nbhd_size > max_nbhd_size && t >= 25) break;
	}

	if (t == max_len) {
		//++n_exhaustive;
		double prop_extend = log_sum(T_nbhds[max_len].size, T_nbhds[max_len].alphas, 
				tnext_argmax_alpha, tnext_max_alpha);
		mempool_destroy(mp);

		return prop_extend;
	}
	else if (T_nbhds[t+1].size > 0) {
		//++n_monte_carlo;

		// perform Monte carlo ahead.
		kbits_t kid, base;
		double prop_extend = 0.0;

		int n = T_nbhds[t+1].size;
		double *ll = (double *) mempool_nalloc(mp, 
				sizeof(*ll) * n , 16);
		memcpy(ll, T_nbhds[t+1].alphas, sizeof(*ll) * n);
		for (int i = 0; i < n; ++i) {
			register double _ll = EXP(ll[i]);
			ll[i] = _ll;
		}
		int64_t *lord = sort_order(ll, n, NULL);

		double sum_ll = 0.0;
		for (int i = 0; i < n; ++i) {
			sum_ll += ll[i];
			ll[i] = sum_ll;
		}

		kmer_t **kptrs = T_nbhds[t+1].kmer_ptrs;
		double *exp_trans_p = d->exp_trans_p;
		kmer_ext_t *nonuniq_k = d->nonuniq_kmers;
		unsigned int tid = omp_get_thread_num();

		for (int i = 0; i < N; ++i) {
			int terminated = 0;
			int llidx = random_state(ll, n, sum_ll, &tid);

			int idx = lord[llidx];
			kmer_t *pkn = kptrs[idx];
			double chain_ext_prob = T_nbhds[t+1].alphas[idx]; 

			kid = kbits_id_of(pkn);

			for (int tp = t + 1; tp <= max_len; ++tp) {
				if (pkn->n_trans == 0 || pkn->alt_trans_flag == 0) {
					terminated = 1;
					break;
				}
				else {
					double *trans_p = exp_trans_p + (
							pkn->unique_trans ? pkn->index : 
							nonuniq_k[pkn->index].poff);
					int bidx = random_transition(trans_p, pkn->n_trans, &tid);
					if (bidx < 0 || IS_ZERO(trans_p[bidx])) {
						terminated = 1;
						break;
					}
					base = (trans_flag2packed(pkn->trans_flag) >> 
						(3 + 2 * bidx)) & 3;
				}

				kid = (kid >> 2) | (base << shift);
				pkn = kptr_locate(d, kid);
			}

			if (terminated == 0) //prop_extend += 1.0 / N; // EXP(chain_ext_prob) / N;
				prop_extend += EXP(chain_ext_prob) / N;
		}

		free(lord);

		mempool_destroy(mp);
		return LOG(prop_extend);
	}
	else {
		// the kmer cannot extend further than the target length, implying that
		// *ALL* kmers along the pathways cannot extend long enough as well.
		/*
		for (int tp = 0; tp <= t; ++tp) {
			state_nbhd_t *nbhd = &T_nbhds[tp];
			for (int i = 0; i < nbhd->size; ++i) {
				kmer_t *pk = nbhd->kmer_ptrs[i];
				if (pk->dirty == 0 && pk->n_trans > 0) {
					d->prob_kmer_extension[kmer_initial_index2(pk, d->uniq_kmer_init_off)] = NAN;
					pk->dirty = 1;
				}
			}
		}
		*/

		mempool_destroy(mp);

		return NAN;
	}
}

void mstep_extend_all_kmers(data *d)
{
	int64_t n_kmers = d->kmer_dict->sizemask + 1;

	d->prob_kmer_extension = realloc(d->prob_kmer_extension,
			(d->n_nonzero_kmers + 1) * sizeof(*d->prob_kmer_extension));
	memset(d->prob_kmer_extension, 0, 
			(d->n_nonzero_kmers + 1) * sizeof(*d->prob_kmer_extension));
	d->prob_kmer_extension[0] = NAN;

	double ts = omp_get_wtime();

	int g_idx = 0, g_prev_idx = 0;

	omp_set_num_threads(d->opt->n_threads);
	#pragma omp parallel 
	{
		int th_idx = 0;
		int tid = omp_get_thread_num();

		dictentry *kdict_tbl = d->kmer_dict->table;

		#pragma omp for schedule(static) 
		for (int64_t i = 0; i < n_kmers; ++i) {
			// report progress
			if (th_idx >= 10000) {
				#pragma omp atomic update
				g_idx += th_idx;

				th_idx = 0;
			}
			else 
				++th_idx;

			if (tid == 0 && (g_idx - g_prev_idx >= 50000)) {
				g_prev_idx = g_idx;

				double duration = omp_get_wtime() - ts;
				io_draw_progress_bar(g_idx, n_kmers, duration, 0, 0, "Kmers");
			}

			dictentry *de = &kdict_tbl[i];
			if (de->value != NULL) {
				kmer_t *pk = (kmer_t *) &(de->value);
				if (pk->n_trans > 0 && pk->dirty == 0) 
					mstep_extend_kmer(d, pk, 0);
			}
		}
	}

	double duration = omp_get_wtime() - ts;
	io_draw_progress_bar(0, n_kmers, duration, 0, 1, "Kmers");


	// reset all dirty flags
	int n_kmers_terminal = 0;
	for (int i = 0; i < n_kmers; ++i) {
		dictentry *de = &d->kmer_dict->table[i];
		if (de->value == NULL) continue;

		kmer_t *pk = (kmer_t *) &(de->value);

		pk->dirty = 0;
		if (pk->n_trans == 0) {
			++n_kmers_terminal;
			continue;
		}

		double *prob_extend = &d->prob_kmer_extension[
				kmer_initial_index2(pk, d->uniq_kmer_init_off)];
		if (IS_ZERO(*prob_extend)) {
			*prob_extend = LOG(1e-10);
			++n_kmers_terminal;
		}
	}

	printf("Total kmers: %zu, kmers that can't extend: %d\n",
			d->n_nonzero_kmers, n_kmers_terminal);
}

void mstep_strand_specific_errors(data *d)
{
	//double thresh = d->opt->penalty_eta / log1p(1.0/d->opt->penalty_gamma);

	//int kmer_len = d->opt->kmer_length;
	//int shift = (kmer_len - 1) * 2;

	double jnt_exp_trans_p[4] = {0.0, 0.0, 0.0, 0.0};
	double recip_exp_trans[4] = {0.0, 0.0, 0.0, 0.0};
	int64_t n_kmers = d->kmer_dict->sizemask + 1;
	const double pval_thresh = 1e-16;

	int n_trans_zeroed_out = 0;
	for (int64_t i = 0; i < n_kmers; ++i) {
		dictentry *de = &d->kmer_dict->table[i];
		if (de->value == NULL) continue;

		kmer_t *pk = (kmer_t *) &(de->value);
		int n = pk->n_trans;
		//kbits_t kid = kbits_id_of(pk);

		// if (pk->label == 1) continue;

		double *exp_trans_p = d->exp_trans_p + (pk->unique_trans ? 
				pk->index : d->nonuniq_kmers[pk->index].poff);
		double *cur_tp = d->transition_p + (pk->unique_trans ?
				0 : d->nonuniq_kmers[pk->index].poff);
		//uint32_t ssflag = mstep_reciprocal_exp_trans(d, pk, &recip_exp_trans[0]);

		double max_trans_j = 0.0;
		for (int j = 0; j < n; ++j) {
			jnt_exp_trans_p[j] = exp_trans_p[j] + recip_exp_trans[j];

			if (jnt_exp_trans_p[j] > max_trans_j)
				max_trans_j = jnt_exp_trans_p[j];
		}

		// ------- EXPERIMENTAL TEST CODE -----
		if (n > 1) {
			int tp_adjusted = 0;

			for (int j = 0; j < n; ++j) {
				double x = MIN(exp_trans_p[j], recip_exp_trans[j]);
				double pval = pnorm(x, 0.5 * jnt_exp_trans_p[j],
						sqrt(0.25 * jnt_exp_trans_p[j]), 1, 0);

				if ((pval < pval_thresh) && 
						(jnt_exp_trans_p[j] != max_trans_j) && x < 50) { //&&
						//(exp_trans_p[j] < thresh || recip_exp_trans[j] < thresh)) {
					tp_adjusted = 1;
					cur_tp[j] = LOG(1e-60);
					++n_trans_zeroed_out;

#if 0
					// adjust the transition probabilities for the other label as well.
					kbits_t base = (trans_flag2packed(pk->trans_flag) >> (3 + j * 2)) & 3;
					kmer_t *pkr = kptr_locate(d, 
							kmer_reverse_complement((kid >> 2) | (base << shift), kmer_len));
					int tidx = base_to_transidx(pk->trans_flag, complementary_base(kid & 3));

					if (pkr->n_trans <= 1) {
						pkr->trans_flag = 0;
						pkr->index = 0;
						pkr->n_trans = 0;
					}
					else {
						double *tpr = d->transition_p + d->nonuniq_kmers[pkr->index].poff;
						tpr[tidx] = LOG(1e-60);
						renormalize_trans_p(pkr->n_trans, tpr);
					}
#endif

				}
			}

			if (tp_adjusted) 
				renormalize_trans_p(n, cur_tp);
		}
	}

	INFO_TEE(d->log_fp, "Transitions considered as strand-specific errors: %d\n",
			n_trans_zeroed_out);
}

/* Conditional maximization of transition parameters */
void cmstep_penalized_transition(data *d, int pfunc, int flip_label)
{
	double jnt_exp_trans_p[4] = {0.0, 0.0, 0.0, 0.0};
	double recip_exp_trans[4] = {0.0, 0.0, 0.0, 0.0};
	double penalized_tp[4] = {0.0, 0.0, 0.0, 0.0};
	double _etas[4];
	double thresh = d->opt->penalty_eta / log1p(1.0/d->opt->penalty_gamma);

	size_t n_uniq_kmers = 0;
	size_t n_nonuniq_kmers = 0;
	size_t n_total_trans = 0;
	size_t n_kmers = d->kmer_dict->used;
	int kmer_len = d->opt->kmer_length;
	//int shift = (kmer_len - 1) * 2;

	double sum_k_cexp_cnt = 0.0;
	double loglik_penalty = 0.0;

	mstep_data_t md = {
		.gamma = d->opt->penalty_gamma,
		.eta = _etas, 
		.xi = NULL,
		.signs = 0,
		.n_xi = 0,
	};

	int lbl_free_param = 0, lbl_dep_param = 1;
	if (flip_label) {
		lbl_free_param = 1;
		lbl_dep_param = 0;
	}

	d->n_trans_removed_below_rho = 0;

	/* The algorithm proceeds as follows:
	 *  1) For label 0 (forward strand) kmers, the transition probabilities
	 *	   (free parameters) will be estimated using MPLE.
	 *  2) When all estimation is finished for forward kmers, enter another
	 *     loop which processes reverse kmers only, computing their transition
	 *     probabilities and the implied transition flags during this process.
	 *  3) Update initial state distribution parameters.
	 */

	int n_fwd = 0;

	double *_init_p = calloc(d->n_nonzero_kmers + 1,
			sizeof(*_init_p));

	//dict *kdict = d->kmer_dict;

	int n_trans_zeroed_out = 0;

	dictentry *de = d->kmer_dict->table;
	INFO_TEE(d->log_fp, "Estimating transition parameters for \"forward\" strand kmers...\n");
	for (int i = 0; i < n_kmers; de++) {
		if (de->key == NULL) continue;
		i++;

		kmer_t *pk = (kmer_t *) &de->value;
		int n = pk->n_trans;
		//kbits_t kid = kbits_id_of(pk);

		/* reverse complement kmer */
		/*
		kbits_t rc_kid = kmer_reverse_complement(kid, kmer_len);
		kmer_t *rc_pk = (kmer_t *) &(dict_find(kdict, 
					kbits_cast_to_ptr(rc_kid))->value);
		*/

		if (pk->label == lbl_dep_param) continue;
		
		n_fwd++;

		double *exp_trans_p = d->exp_trans_p + (pk->unique_trans ? 
				pk->index : d->nonuniq_kmers[pk->index].poff);
		double *cur_tp = d->transition_p + (pk->unique_trans ?
				0 : d->nonuniq_kmers[pk->index].poff);

		double strand_cexp_cnt = 0.0;
		double in_exp_trans = 0.0;
		double sum_exp_trans = 0.0;

		double max_trans_j = 0;

		// initialize all penalty parameters to the default
		for (int j = 0; j < n; ++j) {
			_etas[j] = d->opt->penalty_eta;
		}

		//uint32_t ssflag = mstep_reciprocal_exp_trans(d, pk, &recip_exp_trans[0]);

		for (int j = 0; j < n; ++j) {
			strand_cexp_cnt += exp_trans_p[j];
			jnt_exp_trans_p[j] = exp_trans_p[j] + recip_exp_trans[j];

			sum_exp_trans += jnt_exp_trans_p[j];
			in_exp_trans += recip_exp_trans[j];

			if (jnt_exp_trans_p[j] > max_trans_j)
				max_trans_j = jnt_exp_trans_p[j];

			// EXPERIMENTAL: FIXME!!!
			// if (ssflag & (1 << j)) jnt_exp_trans_p[j] /= 2.0;
		}

		double k_cexp_cnt = 0.0;
		tflag_t new_tf = {pk->trans_flag, n};

		// ------- EXPERIMENTAL TEST CODE -----
		if (n > 1) {
			for (int j = 0; j < n; ++j) {
				double x = MIN(exp_trans_p[j], recip_exp_trans[j]);
				double pval = pnorm(x, 0.5 * jnt_exp_trans_p[j],
						sqrt(0.25 * jnt_exp_trans_p[j]), 1, 0);
				if ((pval < 1e-10) && 
						(jnt_exp_trans_p[j] != max_trans_j) && x < 50) { //&&
						//(exp_trans_p[j] < thresh || recip_exp_trans[j] < thresh)) {
					//raise(SIGINT);
					++n_trans_zeroed_out;
				}
					// raise(SIGINT);
			}
		}
		// ------- EXPERIMENTAL TEST CODE -----

		/* pass in current estimates for trans. prob. to compute the
		 * penalty term. */
		memcpy(&penalized_tp[0], cur_tp, n * sizeof(double));

		/* estimating initial state distribution parameters */
		/*
		double rc_recip_exp_trans[4];
		mstep_reciprocal_exp_trans(d, rc_pk, &rc_recip_exp_trans[0]);
		for (int j = 0; j < rc_pk->n_trans; ++j) {
			in_exp_trans += rc_recip_exp_trans[j];
		}
		*/

		/*
		md.eta = ssflag ? d->opt->penalty_eta * 2 : d->opt->penalty_eta;
		thresh = md.eta / log1p(1.0/d->opt->penalty_gamma);
		*/

		if (pfunc > 0) {
			/* use penalized estimation */
			loglik_penalty += cmstep_estimate_trans_p(d, pk, &penalized_tp[0], 
					&jnt_exp_trans_p[0], &k_cexp_cnt, pk->trans_flag, thresh, 
					&md, &new_tf); 
		}
		else {
			double sum_exp_trans_p = 0.0;
			for (int j = 0; j < n; ++j) {
				//sum_exp_trans_p += exp_trans_p[j];	
				sum_exp_trans_p += jnt_exp_trans_p[j];	
			}

			k_cexp_cnt = sum_exp_trans_p;

			for (int j = 0; j < n; ++j) {
				//penalized_tp[j] = LOG(exp_trans_p[j] / sum_exp_trans_p);
				penalized_tp[j] = LOG(jnt_exp_trans_p[j] / sum_exp_trans_p);
			}
		}


		/* overwrite the expected count */
		if (n > 0) {
			//exp_trans_p[0] = k_cexp_cnt; //strand_cexp_cnt;
			sum_k_cexp_cnt += k_cexp_cnt; //strand_cexp_cnt;
			memcpy(exp_trans_p, jnt_exp_trans_p, n * sizeof(double));

			_init_p[kmer_initial_index(pk, d->uniq_kmer_init_off)] = k_cexp_cnt;
		}


		/* update kmer attributes */
		pk->prev_n_trans = n;
		pk->dirty = 1;
		pk->trans_flag = new_tf.trans_flag;
		pk->unique_trans = (new_tf.n_nonzero_trans <= 1);
		pk->n_trans = new_tf.n_nonzero_trans;

		if (new_tf.n_nonzero_trans == 1) {
			n_uniq_kmers++;
		}
		else if (new_tf.n_nonzero_trans > 1) {
			n_nonuniq_kmers++;
			
			/* update transition probability in place */
			double *kmer_trans_p = d->transition_p + 
				d->nonuniq_kmers[pk->index].poff;
			for (int j = 0; j < n; j++) {
				if (!IS_ZERO(penalized_tp[j])) {
					*kmer_trans_p = penalized_tp[j];
					++kmer_trans_p;
				}
			}

			n_total_trans += new_tf.n_nonzero_trans;
		}
	}

	INFO_TEE(d->log_fp, "Total kmers: %d, forward kmers: %d\n", n_kmers, n_fwd);

	INFO_TEE(d->log_fp, "Transitions considered as strand-specific errors: %d\n",
			n_trans_zeroed_out);

	de = d->kmer_dict->table;
	INFO_TEE(d->log_fp, "Computing transition parameters for \"reverse\" strand kmers...\n");
	int n_zero_mu = 0, n_mu = 0;
	for (int i = 0; i < n_kmers; de++) {
		if (de->key == NULL) continue;

		i++;

		kmer_t *pk = (kmer_t *) &de->value;
		if (pk->label == lbl_free_param) continue;

		int n = pk->n_trans;
		//kbits_t kid = kbits_id_of(pk);

		++n_mu;
		double mu = d->avg_initial_p[kmer_initial_index(pk,
				d->uniq_kmer_init_off)];

		if (isnan(mu)) ++n_zero_mu;

		double *exp_trans_p = d->exp_trans_p + (pk->unique_trans ? 
				pk->index : d->nonuniq_kmers[pk->index].poff);
		/*double *cur_tp = d->transition_p + (pk->unique_trans ?
				0 : d->nonuniq_kmers[pk->index].poff);*/

		double k_cexp_cnt = 0.0;
		double sum_exp_trans_p = 0.0;
		double rev_cexp_cnt;
		//double _penalty = 0.0;

		//uint32_t ssflag = mstep_reciprocal_exp_trans(d, pk, &recip_exp_trans[0]);
		for (int j = 0; j < n; ++j) {
			// _penalty += thresh * log1p(EXP(cur_tp[j]) / md.gamma);
			jnt_exp_trans_p[j] = exp_trans_p[j] + recip_exp_trans[j];
			sum_exp_trans_p += jnt_exp_trans_p[j];
			k_cexp_cnt += recip_exp_trans[j]; // use expected counts from the label 1 kmer.
			// k_cexp_cnt += jnt_exp_trans_p[j];
		}

		tflag_t new_tf = {0, 0};
		double sum_tp = mstep_compute_reverse_trans(d, pk, &penalized_tp[0], 
				mu, kmer_len, &new_tf, _init_p, &rev_cexp_cnt);

		/* renormalize */
		for (int j = 0; j < n; ++j) {
			penalized_tp[j] = LOG(penalized_tp[j] / sum_tp);
		}

		if (n > 0) {
			// exp_trans_p[0] = k_cexp_cnt;
		
			// for unique label 2 transitions, use alternative formula to
			// compute the initial state probability to ensure equality.
			if (n == 1) {
				/*
				uint64_t base = (trans_flag2packed(pk->trans_flag) >> 3) & 3;
				kbits_t _id = kmer_reverse_complement((kid >> 2) | (base << shift), kmer_len);
				kmer_t *_pk = (kmer_t *) &(dict_find(d->kmer_dict, kbits_cast_to_ptr(_id))->value);

				k_cexp_cnt = _init_p[kmer_initial_index(_pk, d->uniq_kmer_init_off)];
				if (_pk->n_trans > 1) {
					int _tidx = base_to_transidx(_pk->trans_flag, 
							(uint32_t) complementary_base(kid & 3));
					k_cexp_cnt *= EXP(*(d->transition_p +
						d->nonuniq_kmers[_pk->index].poff + _tidx));
				}
				*/
				k_cexp_cnt = rev_cexp_cnt;
			}

			sum_k_cexp_cnt += k_cexp_cnt;
			_init_p[kmer_initial_index(pk, d->uniq_kmer_init_off)] = k_cexp_cnt;
		}
		//else {
		//	_penalty = d->opt->penalty_eta; 
		//}

		pk->prev_n_trans = n;
		pk->dirty = 1;
		pk->trans_flag = new_tf.trans_flag;
		pk->unique_trans = (new_tf.n_nonzero_trans <= 1);
		pk->n_trans = new_tf.n_nonzero_trans;

		if (new_tf.n_nonzero_trans == 1) {
			n_uniq_kmers++;
		}
		else if (new_tf.n_nonzero_trans > 1) {
			n_nonuniq_kmers++;
			
			/* update transition probability in place */
			double *kmer_trans_p = d->transition_p + 
				d->nonuniq_kmers[pk->index].poff;
			for (int j = 0; j < n; j++) {
				if (!IS_ZERO(penalized_tp[j])) {
					*kmer_trans_p = penalized_tp[j];
					++kmer_trans_p;
				}
			}

			n_total_trans += new_tf.n_nonzero_trans;
		}
	}

	d->n_total_trans = n_total_trans + n_uniq_kmers;
	d->n_nonuniq_trans = n_total_trans;

	/* allocate memory for new sets of parameters */
	kmer_ext_t *_nonuniq_kmers = calloc(
			n_nonuniq_kmers, sizeof(*d->nonuniq_kmers));

	double *_avg_init_p = calloc(n_uniq_kmers + n_nonuniq_kmers + 1,
			sizeof(*_avg_init_p));
#if PMR_EXTEND_READS
	double *_prob_kmer_ext = calloc(n_uniq_kmers + n_nonuniq_kmers + 1,
			sizeof(*_prob_kmer_ext));
#endif

	d->n_nonzero_kmers = n_uniq_kmers + n_nonuniq_kmers;

	double *_transition_p = calloc(n_total_trans + 1, sizeof(*_transition_p));
	_transition_p[0] = NAN;

	const int _uoff = n_nonuniq_kmers - n_total_trans;

	//INFO("Adjusting indexing..\n");
	/* adjust the index accordingly for unique kmers */
	int n_valid_kmers = 0;
	size_t _kidx = 0, _poff = 1, _uidx = n_total_trans;
	de = d->kmer_dict->table;
	for (int i = 0; i < n_kmers; de++) {
		if (de->value == NULL) continue;

		kmer_t *pk = (kmer_t *) &de->value;
		kbits_t kid = kbits_id_of(pk);

		++i;
		if (pk->dirty == 0) continue;
		//if (pk->label == 1) continue;

		int _par_off = pk->prev_n_trans <= 1 ? 
				pk->index : d->nonuniq_kmers[pk->index].poff;

		/*
		double k_initp = pk->prev_n_trans == 0 ? NAN :
			_penalized_init_p[kmer_initial_index(pk, d->uniq_kmer_init_off)];
		*/
		double k_cexp_cnt = pk->prev_n_trans == 0 ? 0.0 : 
			_init_p[kmer_initial_index(pk, d->uniq_kmer_init_off)];

#if PMR_EXTEND_READS
		double k_prob_ext = d->prob_kmer_extension[
			kmer_initial_index(pk, d->uniq_kmer_init_off)];
#endif

		if (k_cexp_cnt > thresh) ++n_valid_kmers;
		//double k_cexp_cnt =  d->exp_trans_p[_par_off];

		/* twin kmer */
		kbits_t rc_kid = kmer_reverse_complement(kid, kmer_len);

		kmer_t *rc_pk = (kmer_t *) &(dict_find(d->kmer_dict, 
					kbits_cast_to_ptr(rc_kid))->value);

		int _rc_par_off = rc_pk->prev_n_trans <= 1 ?
					rc_pk->index : d->nonuniq_kmers[rc_pk->index].poff;
		/*
		double rc_initp = (rc_pk->prev_n_trans == 0) ? NAN : 
			_penalized_init_p[kmer_initial_index(rc_pk, d->uniq_kmer_init_off)];
		*/

		double rc_cexp_cnt = (rc_pk->prev_n_trans == 0) ? 0.0 : 
			_init_p[kmer_initial_index(rc_pk, d->uniq_kmer_init_off)]; //d->exp_trans_p[_rc_par_off];

		if (rc_cexp_cnt > thresh) ++n_valid_kmers;
		//double rc_cexp_cnt =  d->exp_trans_p[_rc_par_off];

#if PMR_EXTEND_READS
		double rc_prob_ext = d->prob_kmer_extension[
			kmer_initial_index(rc_pk, d->uniq_kmer_init_off)];
#endif

		if (pk->n_trans == 1) {
			/* update index */
			pk->index = _uidx++;
		}
		else if (pk->n_trans > 1) {
			memcpy(_transition_p + _poff, 
					d->transition_p + _par_off,
					sizeof(*_transition_p) * pk->n_trans);

			/* update index */
			pk->index = _kidx;
			_nonuniq_kmers[_kidx++].poff = _poff;
			_poff += pk->n_trans;
		}
		else {
			pk->index = 0;
		}

		/* check if the kmer is a palindrome */
		if (kid != rc_kid) {
			if (rc_pk->n_trans == 1) {
				rc_pk->index = _uidx++;
			}
			else if (rc_pk->n_trans > 1) {
				memcpy(_transition_p + _poff, 
						d->transition_p + _rc_par_off,
						sizeof(*_transition_p) * rc_pk->n_trans);

				/* update index */
				rc_pk->index = _kidx;
				_nonuniq_kmers[_kidx++].poff = _poff;
				_poff += rc_pk->n_trans;
			}
			else {
				rc_pk->index = 0;
			}

			pk->dirty = 0;
			rc_pk->dirty = 0;

			double avg_kcexp_cnt = (k_cexp_cnt + rc_cexp_cnt) / 2.0;
			/* EXPERIMENTAL: FIXME!!! */
			if (avg_kcexp_cnt < thresh &&
					(pk->trans_flag == 0 || rc_pk->trans_flag == 0)) {
				// if a kmer is below the threshold, and appears to be the
				// first kmer, push its initial state prob. to zero, 
				// mimicking the penalty on the initial state distribution.
				avg_kcexp_cnt = d->opt->penalty_gamma;
			}
			/* END OF EXPERIMENTAL CODE: FIXME!! */

			/* compute the average initial probability */
			double mean_initial_p = LOG(avg_kcexp_cnt / sum_k_cexp_cnt);

			int _idx = kmer_initial_index(pk, _uoff);
			int _ridx = kmer_initial_index(rc_pk, _uoff);
			if (_idx > 0) {
				_avg_init_p[_idx] = mean_initial_p;
#if PMR_EXTEND_READS
				_prob_kmer_ext[_idx] = k_prob_ext;
#endif
			}
			if (_ridx > 0) {
				_avg_init_p[_ridx] = mean_initial_p;
#if PMR_EXTEND_READS
				_prob_kmer_ext[_ridx] = rc_prob_ext;
#endif
			}
		}
		else {
			/* index updated */
			pk->dirty = 0;

			double mean_initial_p = LOG(k_cexp_cnt / sum_k_cexp_cnt);
			int _idx = kmer_initial_index(pk, _uoff);
			if (_idx > 0) {
				_avg_init_p[_idx] = mean_initial_p;
#if PMR_EXTEND_READS
				_prob_kmer_ext[_idx] = k_prob_ext;
#endif
			}
		}
	}

	d->uniq_kmer_init_off = _uoff;
	d->prop_nonuniq_kmers = (double) n_nonuniq_kmers / n_valid_kmers * 100.0;
	d->n_nonuniq_kmers = n_nonuniq_kmers;

	INFO_TEE(d->log_fp, "Non-unique kmers: %d / (%.2g%% | V: %.2g%%)\n", n_nonuniq_kmers, 
			(double) n_nonuniq_kmers / n_kmers * 100.0, d->prop_nonuniq_kmers);

	INFO_TEE(d->log_fp, "# transitions removed (when all expected transitions are below threshold): %d\n",
			d->n_trans_removed_below_rho);

	_avg_init_p[0] = NAN;
#if PMR_EXTEND_READS
	_prob_kmer_ext[0] = NAN;
#endif
	
	//INFO("Reallocating memory blocks...\n");
	/* reallocate the memory blocks */
	d->exp_trans_p = realloc(d->exp_trans_p, 
			d->n_total_trans * sizeof(*d->exp_trans_p));

	d->total_kmer_exp_counts = sum_k_cexp_cnt;

	free(d->transition_p);
	d->transition_p = _transition_p;
	free(d->avg_initial_p);
	d->avg_initial_p = _avg_init_p;
#if PMR_EXTEND_READS
	free(d->prob_kmer_extension);
	d->prob_kmer_extension = _prob_kmer_ext;
#endif
	free(d->nonuniq_kmers);
	d->nonuniq_kmers = _nonuniq_kmers;
	free(_init_p);

	// d->loglik -= loglik_penalty;
}

/* 
 * Utility function, tries to use bisection algorithm to find the Lagrange
 * multiplier lambda.
 */
double cmstep_try_bisection(mstep_data_t *d2, double lamb_lb,
		int32_t pos_signs, int *signs, 
		msp_lcstr_func plfunc, msp_lcstr_derv plderv)
{
	mstep_data_t _d1;
	mstep_data_t *d1 = &_d1;
	memcpy(d1, d2, sizeof(_d1));
	d1->signs = pos_signs; 

	double lo_delta = 1e-12;
	double lambda = NAN;

	/* (\sum_j p_{ij}) - 1 at lower bound */
	double lin_cstr_lo = plfunc(lamb_lb - lamb_lb * lo_delta, d2);

	int iter = 0;
	while (isnan(lin_cstr_lo) && (iter++) < 60) {
		lo_delta = pow(lo_delta, 0.5);
		lin_cstr_lo = plfunc(lamb_lb - lamb_lb * lo_delta, d2);
	}
	
	if (isnan(lin_cstr_lo)) return NAN;

	double bound_sign = lin_cstr_lo * plfunc(-1e-10, d2);

	if (!isnan(bound_sign) && bound_sign < 0) {
		lambda = mstep_bisection(plfunc,
				lamb_lb - lamb_lb * lo_delta, -1e-10, 
				1e-10, d2);
		*signs = d2->signs;
	}
	else if (plfunc(lamb_lb - lamb_lb * lo_delta, d1) *
			plfunc(-1e-10, d1) < 0) {
		lambda = mstep_bisection(plfunc,
				lamb_lb - lamb_lb * lo_delta, -1e-10, 
				1e-10, d1);
		*signs = d1->signs;
	}

	return lambda;
}

double cmstep_estimate_trans_p(data *d, kmer_t *pk,
		double *tp, double *exp_trans_p, 
		double *initial_p, uint64_t trans_flag, double thresh, 
		mstep_data_t *md, tflag_t *tf) 
{
	static int cannot_solve = 0;

	kbits_t k_trans_packed = trans_flag2packed(trans_flag);
	kbits_t n_trans = k_trans_packed & 7;

	double _k_penalty = 0.0, max_xi = 0.0;
	int argmax_xi = -1;
	int32_t pos_signs = (1 << n_trans) - 1;
	int32_t signs = pos_signs, alt_signs = pos_signs;
	uint32_t new_tf = 0;

	double orig_tp[4];

	msp_lcstr_func _plfunc = mstep_pen1_lagrange_cstr_func;
	msp_lcstr_derv _plderv = mstep_pen1_lagrange_cstr_derv;

	/*
	if (pen2) {
		_plfunc = mstep_pen2_lagrange_cstr_func;
		_plderv = mstep_pen2_lagrange_cstr_derv;
	}
	*/

	memcpy(orig_tp, tp, n_trans * sizeof(double));

	double sum_xi = 0.0;
	int n_tie = 0; 
	int32_t tie_flag = 0;
	for (int j = 0; j < n_trans; j++) {
		double rho = md->eta[j] / log1p(1/md->gamma);
		sum_xi += exp_trans_p[j];
		/* put penalty on the average of two reciprocal transitions */
		_k_penalty += rho * log1p(EXP(tp[j]) / md->gamma);

		if (exp_trans_p[j] > max_xi) {
			max_xi = exp_trans_p[j];
			argmax_xi = j;
		}
	}


	if (initial_p)
		*initial_p = sum_xi;

	if (n_trans == 0) {
		tf->n_nonzero_trans = 0;
		tf->trans_flag = 0;

		return d->opt->penalty_eta; 
	}
	else if (n_trans == 1) {
		tf->n_nonzero_trans = 1;
		tf->trans_flag = trans_flag; 

		tp[0] = 0.0;

		++d->n_uniq_trans_kmer;
		return _k_penalty;
	}
	else {
		for (int j = 0; j < n_trans; ++j) {
			if (exp_trans_p[j] == max_xi) {
				tie_flag |= (1 << j);
				++n_tie;
			}
		}
	}

	d->n_composite_trans += n_trans;

	/* solve for lambda in the Lagrange multiplier that 
	 * imposes constraint on transition probabilities,
	 * i.e. sum_j (p.ij = 1) for any i. */
	md->xi = exp_trans_p;
	md->n_xi = n_trans;

	if (argmax_xi >= 0 && n_tie <= 1) 
		alt_signs = pos_signs & ~(1U << argmax_xi);
	else
		alt_signs = ~(tie_flag) & ((1 << n_trans) - 1);

	/* compute support for lambda */
	double lamb_lb = mstep_pen1_lambda_support(md);

	double l0 = 0.0;

	md->signs = pos_signs;
	if (max_xi - thresh < 0) {
		/* negative lambda */

		/* EXPERIMENTAL! FIXME: do not apply the penalty */
		/*
		tf->n_nonzero_trans = n_trans;
		tf->trans_flag = trans_flag; 

		for (int i = 0; i < n_trans; ++i) {
			tp[i] = LOG(exp_trans_p[i] / sum_xi);
		}

		return _k_penalty;
		*/

		l0 = lamb_lb + 1e-6;
		/*
		if (pen2) {
			signs = mstep_pen2_determine_signs(lamb_lb, md);
		}
		else {
			signs = mstep_pen1_determine_signs(lamb_lb, md);
		}
		*/
		signs = mstep_pen1_determine_signs(lamb_lb, md);
		md->signs = signs;
	}
	else {
		l0 = max_xi - thresh;
		signs = pos_signs;
		md->signs = pos_signs;
	}

	double lambda = mstep_newton(_plfunc, _plderv, l0, 1e-10, 100, md);
	double cstr = _plfunc(lambda, md);

	if (isnan(lambda) || cstr < -1e-1) {
		/* Try bisection method */
		mstep_data_t d2 = {.gamma = md->gamma, .eta = md->eta,
			.xi = exp_trans_p, .signs = alt_signs, .n_xi = n_trans};

		lambda = cmstep_try_bisection(&d2, lamb_lb, pos_signs, &signs,
				_plfunc, _plderv);

		/* still can't find solution to lambda ... ? */
		if (isnan(lambda)) {
			d2.signs = pos_signs;

			/* guessing the lower bound for bisection algorithm */
			double _ub = LOG(sum_xi), _lb = LOG(1e-12);
			double _step = (_ub - _lb) / 100;
			for(; _lb < _ub; _lb += _step) {
				double _cstr = _plfunc(EXP(_lb), &d2);
				if (_cstr > 0) break;
			}
			/* trying bisection for positive lambda */
			lambda = mstep_bisection(_plfunc, 
					EXP(_lb), sum_xi,
					1e-10, &d2);
		}

		/* okay... then it must be we need to set another sign to negative,
		 * likely the one with the largest recip_tp. */
		if (isnan(lambda)) {
			mstep_data_t d2 = {.gamma = md->gamma, .eta = md->eta,
				.xi = exp_trans_p, .signs = signs, .n_xi = n_trans};

			for (int j = 0; j < n_trans && isnan(lambda); j++) {
				if (j == argmax_xi) continue;

				d2.signs = pos_signs & ~(1U << j);
				lambda = cmstep_try_bisection(&d2, lamb_lb, pos_signs, &signs,
						_plfunc, _plderv);
			}
		}
	}


	int n_nonzero_trans = 0;
	if (!isnan(lambda)) {

#ifdef PMR_CHECK_ASCENT
		double pll_m = 0.0, pll_mp1 = 0.0;
#endif

		double sum_tp = 0.0;
		for (int base = 0; base < n_trans; base++) {
			double rho = md->eta[base] / log1p(1/md->gamma);
			/* compute the MPLE for transition prob. */
			double _xi = exp_trans_p[base];
			double b = md->gamma * lambda - _xi + rho; 
			double sign = (double) ((((signs >> base) & 1) << 1) - 1);

#ifdef PMR_CHECK_ASCENT
			//double _mle_tp = _xi / sum_xi;
			//pll_m += _mle_tp * _xi - thresh * log(1+_mle_tp/md->gamma);
			pll_m += tp[base] * _xi - rho * log(1+EXP(tp[base])/md->gamma);
#endif

			tp[base] = (-b + sign * 
				sqrt(b*b + 4 * lambda * md->gamma * _xi)) /
				(2 * lambda);

			sum_tp += tp[base];
		}

		/* normalization */
		int new_n_trans = 0;
		k_trans_packed >>= 3;
		for (int i = 0; i < n_trans; k_trans_packed >>=2, i++) {
			kbits_t base = k_trans_packed & 3;

			/*
			printf("%lg[%lg]->%lg ", exp_trans_p[i], 
					recip_tp[i], tp[i] / sum_tp);
			*/
			
			tp[i] = LOG(tp[i] / sum_tp);

#ifdef PMR_CHECK_ASCENT
			pll_mp1 += tp[i] * exp_trans_p[i] - rho * log(1+EXP(tp[i])/md->gamma);
#endif

			if (!IS_ZERO(tp[i])) {
				++new_n_trans;
				new_tf |= (1 << base);
				++n_nonzero_trans;
			}
		}

#ifdef PMR_CHECK_ASCENT
		if ((lambda < 0) && (n_trans > 1) //&& (pll_m > pll_mp1)
				&& !isnan(pll_m) && (pll_m - pll_mp1 > 1e-4)
				&& (n_trans == new_n_trans)) {
			printf("penalized likelihood: (m) = %lg, (m+1) = %lg\n",
					pll_m, pll_mp1);
		}
#endif

	}
	else {
		++cannot_solve;
		fprintf(stderr, " ## cannot solve: ");
		for (int i = 0; i < n_trans; ++i) {
			fprintf(stderr, "%lg ", exp_trans_p[i]);
		}
		fprintf(stderr, "\n");
		//fprintf(stderr, "cannot solve: %d\n", cannot_solve);
		//getchar();

		/* hardcoded truth!! */
		/*
		tp[0] = 0.0;
		tp[1] = NAN;
		new_tf = 1;
		n_nonzero_trans = 1;

		for (int i = 0; i < n_trans; k_trans_packed >>=2, i++) {
			kbits_t base = k_trans_packed & 3;

			printf("%lg[%lg]-> x ", exp_trans_p[i], 
					recip_tp[i]);
		}
		*/
	}

	if (max_xi - thresh < 0) {
		d->n_trans_removed_below_rho += (tf->n_nonzero_trans - n_nonzero_trans);
	}

	tf->trans_flag = new_tf;
	tf->n_nonzero_trans = n_nonzero_trans;

	return _k_penalty;
}

void mstep_alt_emission(data *d)
{
	const size_t qmax = (d->opt->max_qual_score + 1);
	const size_t qmax2 = qmax << 1;

	size_t nb_ctx = 0;
	int kmer_len = d->opt->kmer_length;
	int _tw = kmer_len + d->n_emit_windows;

	for (int t = 0; t < _tw; ++t) {
		size_t base_toff = t * 20;
		size_t qual_toff = t * qmax2;
		if (t < kmer_len) {
			nb_ctx = BASE_N; 
		}
		else {
			int w = t - kmer_len;

			base_toff = kmer_len * 20 + w * d->base_ctx_radix;
			nb_ctx = d->n_base_context;

			qual_toff = kmer_len * qmax2 + w * d->qual_ctx_radix;
		}

		// base emission
		for (int i = 0; i < nb_ctx; i++) {
			double *ctx_base_emit_exp = &d->base_emit_exp[base_toff + i * 5];
			double *ctx_base_emit_p = &d->base_emission_p[base_toff + i * 5];

			int base = i & 3;

			double sum_emission_p = 0.0;
			for (int j = 0; j <= BASE_N; j++) {
				if (base != j) sum_emission_p += ctx_base_emit_exp[j];
			}

			/* avoid NAN */
			if (sum_emission_p == 0.0) sum_emission_p = 1.0;

			for (int j = 0; j <= BASE_N; j++) {
				if (j == base) 
					ctx_base_emit_p[j] = 0.0;
				else
					ctx_base_emit_p[j] = LOG(ctx_base_emit_exp[j] / sum_emission_p);
			}
		}

		// quality score emission
		for (int j = 0; j < qmax; ++j) {
			double qj_error = d->qual_emit_exp[qual_toff + j];
			double sum_qual_exp = qj_error + 
				d->qual_emit_exp[qual_toff + qmax + j];

			// error rate
			d->qual_emission_p[qual_toff + j] = LOG(qj_error / sum_qual_exp);
			// 1 - error rate
			d->qual_emission_p[qual_toff + qmax + j] = LOG(1 - (qj_error / sum_qual_exp));
		}

		if (t >= kmer_len) {
			fprintf(stderr, "[Estimated Qual.] <w=%2d>\n", t - kmer_len);
			for (int j = 1; j < qmax; ++j) {
				fprintf(stderr, "%2d %8.3e | ", j, 
						EXP(d->qual_emission_p[qual_toff + j]));
				if (j % 5 == 0) fprintf(stderr, "\n");
			}
			fprintf(stderr, "\n");
		}
	}
}

void mstep_emission(data *d)
{
#if PMR_USE_PRIMITIVE_EMISSION
	INFO("Previous per base error rate: %lg\n", d->base_error_rate);
	d->base_error_rate = (d->error_rate_emit_exp[1] / 
		(d->error_rate_emit_exp[0] + d->error_rate_emit_exp[1]));
	INFO("Estimated per base error rate: %lg\n", d->base_error_rate);
#else
	const size_t qmax = (d->opt->max_qual_score + 1);
	const size_t qmax2 = qmax << 1;

	size_t nb_ctx = 0, nq_ctx = 0;
	int kmer_len = d->opt->kmer_length;
	int _tw = kmer_len + d->n_emit_windows;

	const size_t n_base_emit = 20 * d->opt->kmer_length + 
		d->n_emit_windows * d->base_ctx_radix;
	const size_t n_qual_emit = qmax2 * d->opt->kmer_length +
		d->n_emit_windows * d->qual_ctx_radix;

	for (int i = 1; i < d->opt->n_threads; ++i) {
		// aggregate thread-specific counts
		for (int j = 0; j < n_base_emit; ++j) {
			d->base_emit_exp[j] += d->base_emit_exp[i * n_base_emit + j];
		}

		for (int j = 0; j < n_qual_emit; ++j) {
			d->qual_emit_exp[j] += d->qual_emit_exp[i * n_qual_emit + j];
		}
	}

#if !PMR_NO_DEP_EMISSION
	// prepare a table to report per window, average error rate.
	double *tbl_avg_erate = calloc(d->n_emit_windows * 6, sizeof(double));
	double *_ptbl = tbl_avg_erate;
#endif

	for (int t = 0; t < _tw; ++t) {
		size_t base_toff = t * 20;
		size_t qual_toff = t * qmax2;
		if (t < kmer_len) {
			nb_ctx = BASE_N; 
			nq_ctx = 2;
		}
		else {
			int w = t - kmer_len;

			base_toff = kmer_len * 20 + w * d->base_ctx_radix;
			nb_ctx = d->n_base_context;

			qual_toff = kmer_len * qmax2 + w * d->qual_ctx_radix;
			nq_ctx = d->n_qual_context;
		}

#if !PMR_NO_DEP_EMISSION
		double _sum_total_k = 0.0, _sum_err_free_k = 0.0;
#endif
		for (int i = 0; i < nb_ctx; i++) {
			double *ctx_base_emit_exp = &d->base_emit_exp[base_toff + i * 5];
			double *ctx_base_emit_p = &d->base_emission_p[base_toff + i * 5];

#if !PMR_NO_DEP_EMISSION
			int hd_pwr = i / 4;
			int base = i & 3;
#endif

			double sum_emission_p = 0.0;
			for (int j = 0; j <= BASE_N; j++) {
				sum_emission_p += ctx_base_emit_exp[j];

#if !PMR_NO_DEP_EMISSION
				_sum_total_k += ctx_base_emit_exp[j];
				if (t >= kmer_len && j == base) _sum_err_free_k += ctx_base_emit_exp[j];
#endif
			}

#if !PMR_NO_DEP_EMISSION
			if (t >= kmer_len && base == 3) {
				// last base
				*_ptbl = 1.0 - _sum_err_free_k / _sum_total_k;
				++_ptbl;
				_sum_err_free_k = 0.0;
				_sum_total_k = 0.0;
			}
#endif
			
			/* avoid NAN */
			if (sum_emission_p == 0.0) sum_emission_p = 1.0;

			for (int j = 0; j <= BASE_N; j++) {
				ctx_base_emit_p[j] = LOG(ctx_base_emit_exp[j] / sum_emission_p);
			}
		}

		for (int i = 0; i < nq_ctx; i++) {
			double *ctx_qual_emit_exp = &d->qual_emit_exp[qual_toff + i * qmax];
			double *ctx_qual_emit_p = &d->qual_emission_p[qual_toff + i * qmax];

			double sum_qual_exp = 0.0;
			for (int j = 0; j < qmax; j++) {
				sum_qual_exp += ctx_qual_emit_exp[j];
			}

			if (sum_qual_exp == 0.0) sum_qual_exp = 1.0;

			for (int j = 0; j < qmax; j++) {
				ctx_qual_emit_p[j] = LOG(ctx_qual_emit_exp[j] / 
						sum_qual_exp);
			}
		}
	}

#if !PMR_NO_DEP_EMISSION
	fprintf(stderr, 
			"\n==================[ Base emission / Avg. error rates ]==================\n"
			"       t / pfx.err        0        1      2-3      4-7     8-15      16-\n"
			"------------------------------------------------------------------------\n");

	for (int w = 0; w < d->n_emit_windows; ++w) {
		int w_end = kmer_len + (w+1) * d->opt->emit_window_size;
		if (w_end > d->max_read_length) w_end = d->max_read_length;

		fprintf(stderr, " %3d-%3d           ",
				kmer_len + w * d->opt->emit_window_size + 1,
				w_end);

		for (int hd = 0; hd < 6; ++hd) {
			fprintf(stderr, "%8.5f ", tbl_avg_erate[w * 6 + hd]);
		}
		fprintf(stderr, "\n");
	}

	fprintf(stderr, "========================================================================\n\n");
	free(tbl_avg_erate);
#endif
#endif
}

void mstep_initial(data *d)
{
	/* placeholder for now */
}

#if 0
void mstep_initial(data *d)
{
	if (d->opt->disable_penalized_em) {
		dictentry *de = d->kmer_dict->table;
		for (int i = d->kmer_dict->used; i > 0; de++) {
			if (de->value == NULL) continue;

			i--;
			kmer_t *pk = de->value;
			pk->init_p = LOG(pk->exp_count / d->total_kmer_exp_counts);
		}
	}
	else {
		double *kmer_exp_counts = malloc(sizeof(*kmer_exp_counts) *
				d->kmer_dict->used);
		double thresh = d->opt->penalty_eta / log1p(1/d->opt->penalty_gamma);
		double l0 = d->total_kmer_exp_counts - thresh * d->kmer_dict->used;

		mstep_data_t md = {
			.gamma = d->opt->penalty_gamma,
			/* use eta to not confuse with the lambda in 
			 * Lagrange multiplier */
			.eta = d->opt->penalty_eta, 
			.xi = kmer_exp_counts,
			.signs = 0x0,
			.n_xi = d->kmer_dict->used
		};

		double loglik_penalty = 0.0;

		dictentry *de = d->kmer_dict->table;
		register int imax = d->kmer_dict->used;
		for (int i = 0; i < imax; de++) {
			if (de->value == NULL) continue;

			kmer_t *pk = de->value;
			if (pk->exp_count == 0.0 || isnan(pk->init_p)) {
				pk->transition_p = d->transition_p;
				pk->trans_flag = 0;
			}
			loglik_penalty += thresh * log1p(EXP(pk->init_p) / md.gamma);
			kmer_exp_counts[i++] = pk->exp_count;
		}

		d->loglik -= loglik_penalty;

		double lambda = mstep_newton(mstep_lagrange_cstr_func, 
				mstep_lagrange_cstr_derv, l0, 
				1e-10, 100, &md);
		
		if (!isnan(lambda)) {
			dictentry *de = d->kmer_dict->table;
			int i = d->kmer_dict->used, j = 0;
			for (; i > 0; de++) {
				if (de->value == NULL) continue;

				i--;
				kmer_t *pk = (kmer_t *) de->value;
				double _exp_cnt = kmer_exp_counts[j++];
				
				double b = md.gamma * lambda - _exp_cnt + thresh;
				pk->init_p = LOG((-b + sqrt(b*b + 4 * lambda * 
								md.gamma * _exp_cnt)) / (2 * lambda));
			}
		}

		free(kmer_exp_counts);
	}
}
#endif

double mstep_pen1_lambda_support(void *fdata)
{
	mstep_data_t *md = fdata;
	register double _gamma = md->gamma;
	register int n_xi = md->n_xi;

	register double max_lamb_lb = -INFINITY;
	for (int j = 0; j < n_xi; j++){
		register double _rho = md->eta[j] / log1p(1.0/_gamma);
		register double xi_j = md->xi[j];
		register double _lb = (-(xi_j + _rho) + 2 * sqrt(xi_j * _rho)) / _gamma;
		if (_lb > max_lamb_lb) 
			max_lamb_lb = _lb;
	}

	return max_lamb_lb;
}

#if 0
double mstep_pen2_lambda_support(void *fdata)
{
	mstep_data_t *md = fdata;
	double _gamma = md->gamma;
	double *_xi = md->xi;
	register double _nu = 2 * md->eta / log(1 + 1/_gamma);
	register int n_xi = md->n_xi;

	register double max_lamb_lb = -INFINITY;
	for (int j = 0; j < n_xi; j++){
		double xi_j = _xi[j];
		double _mu = 2 * _gamma + md->recip_tp[j];
		register double _lb = (-(xi_j + _nu) + 2 * sqrt(xi_j * _nu)) / 
			_mu;
		if (_lb > max_lamb_lb) 
			max_lamb_lb = _lb;
	}

	return max_lamb_lb;
}
#endif

int32_t mstep_pen1_determine_signs(double lamb_lb, void *fdata)
{
	register int32_t signs_v = 0xF;
	register mstep_data_t *md = fdata;
	register double *pxi = md->xi;
	register int n_xi = md->n_xi;

	mstep_data_t param = {.gamma = md->gamma, .eta = md->eta,
		.xi = md->xi, .signs = signs_v};

	int argmax_idx = -1;
	double max_xi = 0.0;
	for (int j = 0; j < n_xi; j++) {
		if (pxi[j] > max_xi) {
			max_xi = pxi[j];
			argmax_idx = j;
		}
	}

	if (argmax_idx < 0) return signs_v;

	register int32_t alt_signs_v = signs_v & ~(1 << argmax_idx);

	/* FIXME: not thread-safe */
	double l0 = lamb_lb + EPSILON;
	double l1 = l0 - mstep_pen1_lagrange_cstr_func(l0, &param) / 
		mstep_pen1_lagrange_cstr_derv(l0, &param);

	param.signs = alt_signs_v;

	double l1_alt = l0 - mstep_pen1_lagrange_cstr_func(l0, &param) / 
		mstep_pen1_lagrange_cstr_derv(l0, &param);

	if (lamb_lb < 0 && l1 > 0) 
		return signs_v;

	if (l1_alt > lamb_lb) {
		return alt_signs_v;
	}
	else {
		return signs_v;
	}
}

#if 0
int32_t mstep_pen2_determine_signs(double lamb_lb, void *fdata)
{
	register mstep_data_t *md = fdata;
	register double *pxi = md->xi;
	register int n_xi = md->n_xi;
	register int32_t signs_v = (1 << n_xi) - 1; 

	mstep_data_t param = {.gamma = md->gamma, .eta = md->eta,
		.xi = md->xi, .signs = signs_v, .n_xi = n_xi,
		.recip_tp = md->recip_tp};

	int argmax_idx = -1;
	double max_xi = 0.0;
	for (int j = 0; j < n_xi; j++) {
		if (pxi[j] + md->recip_tp[j] > max_xi) {
			max_xi = pxi[j] + md->recip_tp[j];
			argmax_idx = j;
		}
	}

	if (argmax_idx < 0) return signs_v;

	register int32_t alt_signs_v = signs_v & ~(1 << argmax_idx);

	/* FIXME: not thread-safe */
	double l0 = lamb_lb + EPSILON;
	double l1 = l0 - mstep_pen2_lagrange_cstr_func(l0, &param) / 
		mstep_pen2_lagrange_cstr_derv(l0, &param);

	param.signs = alt_signs_v;

	double l1_alt = l0 - mstep_pen2_lagrange_cstr_func(l0, &param) / 
		mstep_pen2_lagrange_cstr_derv(l0, &param);

	if (lamb_lb < 0 && l1 > 0) 
		return signs_v;

	if (l1_alt > lamb_lb) {
		return alt_signs_v;
	}
	else {
		return signs_v;
	}
}


double mstep_pen2_lagrange_cstr_func(double lambda, void *fdata)
{
	register mstep_data_t *md = fdata;
	register double gam = md->gamma;
	register double thresh = md->eta / log1p(1.0/gam);
	register double *pxi = md->xi;
	register int32_t sgns = md->signs;
	register int n_xi = md->n_xi;
	int sgn_off = n_xi > 4;

	register double sum_root = 0.0;
	for (int j = 0; j < n_xi; j++) {
		register double xij = pxi[j];

		double b = 2 * gam * lambda + lambda * md->recip_tp[j] +
			2 * thresh - xij;
		/* minus c */
		double mc = xij * (2 * gam + md->recip_tp[j]);

		double delta = b * b + 4 * lambda * mc;
		/* sanity check */
		if (delta < 0) return NAN;

		delta = sqrt(delta);

		double sign = (double) (((((sgns >> j) & 1) + sgn_off) << 1) - 1);
		sum_root += (-b + sign * delta)/(2 * lambda);
	}

	return sum_root - 1;
}

double mstep_pen2_lagrange_cstr_derv(double lambda, void *fdata)
{
	register mstep_data_t *md = fdata;
	register double gam = md->gamma;
	register double thresh = md->eta / log1p(1.0/gam);
	register double *pxi = md->xi;
	register int32_t sgns = md->signs;
	register int n_xi = md->n_xi;
	int sgn_off = n_xi > 4;

	register double derv = 0.0;
	for (int j = 0; j < n_xi; j++) {
		register double xij = pxi[j];

		double b = 2 * gam * lambda + lambda * md->recip_tp[j] +
			2 * thresh - xij;
		/* minus c */
		double mc = xij * (2 * gam + md->recip_tp[j]);
		/* with respect to lambda */
		double bprime = 2 * gam + md->recip_tp[j];

		double delta = b * b + 4 * lambda * mc;
		/* sanity check */
		if (delta < 0) return NAN;

		delta = sqrt(delta);
		double sign = (double) (((((sgns >> j) & 1) + sgn_off) << 1) - 1);
		derv += (lambda * (-bprime + sign / delta * 
					(b * bprime + 2 * mc))) - (-b + sign * delta);
	}

	return derv / (2 * lambda * lambda);
}
#endif

double mstep_pen1_lagrange_cstr_func(double lambda, void *fdata)
{
	register mstep_data_t *md = fdata;
	register double gam = md->gamma;
	register int32_t sgns = md->signs;
	register int n_xi = md->n_xi;
	int sgn_off = n_xi > 4;

	register double sum_root = 0.0;
	for (int j = 0; j < n_xi; j++) {
		register double thresh = md->eta[j] / log1p(1.0/gam);
		register double xij = md->xi[j];

		double b = gam * lambda - xij + thresh;
		double mc = gam * xij;

		double delta = b * b + 4 * lambda * mc;
		/* sanity check */
		if (delta < 0) return NAN;

		delta = sqrt(delta);

		double sign = (double) (((((sgns >> j) & 1) + sgn_off) << 1) - 1);
		sum_root += (-b + sign * delta)/(2 * lambda);
	}

	return sum_root - 1;
}


double mstep_pen1_lagrange_cstr_derv(double lambda, void *fdata)
{
	register mstep_data_t *md = fdata;
	register double gam = md->gamma;
	register int32_t sgns = md->signs;
	register int n_xi = md->n_xi;
	int sgn_off = n_xi > 4;

	register double derv = 0.0;
	for (int j = 0; j < n_xi; j++) {
		register double thresh = md->eta[j] / log1p(1.0/gam);
		register double xij = md->xi[j];

		double b = gam * lambda - xij + thresh;
		double mc = gam * xij;

		double delta = b * b + 4 * lambda * mc;
		/* sanity check */
		if (delta < 0) return NAN;

		delta = sqrt(delta);
		double sign = (double) (((((sgns >> j) & 1) + sgn_off) << 1) - 1);
		derv += (lambda * (-gam + sign / delta * gam * (b + 2 * xij))) - 
			(-b + sign * delta);
	}

	return derv / (2 * lambda * lambda);
}

/* ----- Root finding algorithms ----- */
double mstep_newton(double (*fx)(double x, void *data), 
		double (*fderv)(double x, void *data), double x0, double eps,
		int maxiter, void *fdata)
{
	int iter = 1;
	double x1 = x0 - fx(x0, fdata) / fderv(x0, fdata);

	if (isnan(x1)) return NAN;

	while (fabs(x1 - x0) > eps && iter < maxiter) {
		x0 = x1;
		x1 = x0 - fx(x0, fdata) / fderv(x0, fdata);

		if (isnan(x1)) return NAN;

		iter++;
	}

	return x1;
}

/**
 * Finding root of f(x) using bisection method.
 *
 * @param fx function f(x)
 * @param a lower bound
 * @param b upper bound
 * @param eps accepted error
 * @param fdata extra data pass to the function pointer
 */
double mstep_bisection(double (*fx)(double x, void *data), double a, double b, 
		double eps, void *fdata)
{
	double fa = fx(a, fdata),
		   fb = fx(b, fdata);

	double mid, fmid;

	if (fa == 0) return a;
	if (fb == 0) return b;

	if (a > b || (fa * fb) > 0)
		return NAN;

	mid = (a+b) / 2;

	if (b - a < eps)
		return mid;

	fmid = fx(mid, fdata);

	if ((fa < 0 && fmid >= 0) || (fa > 0 && fmid <= 0))
		return mstep_bisection(fx, a, mid, eps, fdata);
	else
		return mstep_bisection(fx, mid, b, eps, fdata);
}	/* bisection */
