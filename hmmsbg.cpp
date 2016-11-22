#include <vector>
#include <algorithm>
#include <iterator>
#include <memory>
#include <cstdio>
#include <unordered_set>
#include <unordered_map>
#include <string>

#include <signal.h>
#include <omp.h>

#include "hmmsbg.h"

extern "C" {
#include "message.h"
#include "iodata.h"
#include "read.h"
#include "numeric.h"
#include "trans.h"
#include "bitarray.h"
#include "atomic.h"
}

static void hmmsbg_coalesce_lbl2_trans(data *d);
static uint64_t hmmsbg_determine_ctx(data *d, read_t *read, kmer_t *pk, int t, 
		bool rev_strand = false);
static uint64_t hmmsbg_load_ctx_trans_p(data *d, read_t *read, kmer_t *pk, int t, 
		double *tp);
static void hmmsbg_setup_nbhd_ptrs(state_nbhd_t *nbhd, int n, 
		char *pmem);
static int hmmsbg_load_preconstructed_nbhd(data *d, read_t *read,
		kbits_t obs_1st_kid, kbits_t obs_n_flag, const char *conv_q, 
		state_nbhd_t *first_nbhd, int kmer_len, int compute_emit_dens, 
		dict *kdict, mempool_t *mp, const double *uniq_tp);
static double hmmsbg_update_lbl2_trans_p(data *d, kmer_t *pk, double *_exp_cnt_lbl1);
static void hmmsbg_reciprocal_exp_trans(data *d, kmer_t *pk, double *recip_exp_trans);
static void hmmsbg_debug_console(data *d);
static void hmmsbg_debug_kmer(data *d, kbits_t kid);

inline void hmmsbg_setup_nbhd_ptrs(state_nbhd_t *nbhd, int n, 
		char *pmem)
{
	memset(pmem, 0, hmmsbg_nbhd_alloc_size(n));

	/* FIXME: consider lookup table and struct initialization */
	nbhd->kmer_ptrs = (kmer_t **) pmem;
	nbhd->kmer_trans_prob = (double **) (pmem + n * sizeof(*nbhd->kmer_ptrs));

	nbhd->ba_kmer_trans_flag = (kbits_t *) (
			(char *) nbhd->kmer_trans_prob + n * 
			sizeof(*nbhd->kmer_trans_prob));
	nbhd->ba_distinct_sfx = nbhd->ba_kmer_trans_flag + BITS_TO_U64(n << 2);
	nbhd->ba_hamming_dist = nbhd->ba_distinct_sfx + BITS_TO_U64(n << 1);
	nbhd->ba_pfx_sfx_flags = nbhd->ba_hamming_dist + BITS_TO_U64(n << 3);

	nbhd->suffixes_order = (int64_t *) (nbhd->ba_pfx_sfx_flags + 
			BITS_TO_U64(n << 3));
	nbhd->states_sorted_pfx = (kbits_t *) (nbhd->suffixes_order + n);
	nbhd->states_sorted_sfx = nbhd->states_sorted_pfx + n;

	nbhd->alphas = (double *) (nbhd->states_sorted_sfx + n);
	nbhd->betas = nbhd->alphas + n;
	nbhd->emit_dens = nbhd->betas + n;
	nbhd->base_emit_ctx = (int64_t *) (nbhd->emit_dens + n);
	nbhd->qual_emit_ctx = nbhd->base_emit_ctx + n;
	nbhd->kmer_ups_ctx = (uint64_t *) nbhd->qual_emit_ctx + n;
	nbhd->jump_ahead_countdown = (int64_t *) nbhd->kmer_ups_ctx + n;
}

double hmmsbg_compute_alpha(data *d, kbits_t ups_trans_packed, 
		int ups_kmer_idx, state_nbhd_t *nbhd, int64_t *countdown,
		// kmer_t **t_kmer_ptrs, double *t_state_alphas, double **t_kmer_trans_p, 
		double emit_dens, int dns_base)
{
	double summands[4] = {NAN};
	int n_ups_states = ups_trans_packed & 7;
	kbits_t utp = ups_trans_packed >> 3;
	dbl_t max_summand = {.dbl = NAN};
	int argmax_idx = 0;

#if PMR_JUMP
	int jump_ahead_countdown = 0;
#endif

	int n = 0;
	for (int i = 0; i < n_ups_states; i++, utp >>= 2) {
		kbits_t ups_base = utp & 3;
		register int _idx = ups_kmer_idx + i;
		
		kmer_t *pk = nbhd->kmer_ptrs[_idx];
		register kbits_t _tf = pk->trans_flag;
		if ((_tf & (1 << dns_base)) == 0)
			continue;

		register int dns_bidx = _mm_popcnt_u64(_tf & ((1 << dns_base) - 1));
		// if (IS_ZERO(nbhd->kmer_trans_prob[_idx][dns_bidx])) continue;

#if PMR_JUMP
		// despite the upstream source, the count down should be exactly the
		// same, hence we can return any one of those.
		if (pk->label == 1 && pk->has_multi_copies == 1 &&
				nbhd->jump_ahead_countdown[_idx] == 0) {
			uint32_t n_ctx_i = (d->nonuniq_kmers[pk->index].n_ctx >> 
					(dns_bidx * PMR_MAX_CTX_BITS)) & PMR_MAX_NCTX;
			uint32_t ctx_i = (nbhd->kmer_ups_ctx[_idx] >>
					(dns_bidx * PMR_MAX_CTX_BITS)) & PMR_MAX_NCTX;

			if (n_ctx_i == ctx_i) jump_ahead_countdown = 0;
			else 
				jump_ahead_countdown = (d->nonuniq_kmers[pk->index].ctx_distance >>
						(dns_bidx * PMR_MAX_CTX_BITS)) & PMR_MAX_NCTX;
		}
		else if (nbhd->jump_ahead_countdown[_idx] > 0) {
			jump_ahead_countdown = nbhd->jump_ahead_countdown[_idx] - 1;
		}
#endif

		register dbl_t _summand = {.dbl = PRODUCT(nbhd->alphas[_idx], 
					nbhd->kmer_trans_prob[_idx][dns_bidx])};

		summands[n] = _summand.dbl;

		/* (greater or equal) ge_mask is set to 0, 
		 * if max_summand < _summand */
		uint64_t ge_mask = (uint64_t) (isnan(max_summand.dbl) || 
				(max_summand.dbl < _summand.dbl)) - 1UL;
		uint64_t inv_ge_mask = ~ge_mask;

		max_summand.u64 = (ge_mask & max_summand.u64) | 
			(inv_ge_mask & _summand.u64);
		argmax_idx = (ge_mask & argmax_idx) | (inv_ge_mask & n);

		++n;
	}

	/* compute the sum of all summands gives us alpha */
	double _expsum = 0.0;
	for (int i = 0; i < n; i++) {
		register uint64_t argmax_mask = ((i == argmax_idx) - 1UL);
		register dbl_t _summand = {
			.dbl = EXP(summands[i] - max_summand.dbl)};
		_summand.u64 &= argmax_mask;

		_expsum += _summand.dbl; 
	}

#if PMR_JUMP
	*countdown = jump_ahead_countdown;
#endif

	return PRODUCT(emit_dens,
			PRODUCT(max_summand.dbl, log1p(_expsum)));
}

double hmmsbg_compute_delta(data *d, kbits_t ups_trans_packed, 
		int ups_kmer_idx, state_nbhd_t *nbhd, 
		int64_t *countdown, double emit_dens, int dns_base, double *delta)
{
	dbl_t max_tail = {.dbl = NAN};
	int argmax_tail = 0;

	int n_ups_states = ups_trans_packed & 7;
	kbits_t utp = ups_trans_packed >> 3;

#if PMR_JUMP
	int jump_ahead_countdown = 0;
#endif

	int n = 0;
	for (int i = 0; i < n_ups_states; i++, utp >>= 2) {
		kbits_t ups_base = utp & 3;
		register int _idx = ups_kmer_idx + i;
		
		kmer_t *pk = nbhd->kmer_ptrs[_idx];
		register kbits_t _tf = pk->trans_flag;
		if ((_tf & (1 << dns_base)) == 0)
			continue;

		register int dns_bidx = _mm_popcnt_u64(_tf & ((1 << dns_base) - 1));
		// if (IS_ZERO(nbhd->kmer_trans_prob[_idx][dns_bidx])) continue;

#if PMR_JUMP
		// despite the upstream source, the count down should be exactly the
		// same, hence we can return any one of those.
		if (pk->label == 1 && pk->has_multi_copies == 1 &&
				nbhd->jump_ahead_countdown[_idx] == 0) {
			uint32_t n_ctx_i = (d->nonuniq_kmers[pk->index].n_ctx >> 
					(dns_bidx * PMR_MAX_CTX_BITS)) & PMR_MAX_NCTX;
			uint32_t ctx_i = (nbhd->kmer_ups_ctx[_idx] >>
					(dns_bidx * PMR_MAX_CTX_BITS)) & PMR_MAX_NCTX;

			if (n_ctx_i == ctx_i) jump_ahead_countdown = 0;
			else 
				jump_ahead_countdown = (d->nonuniq_kmers[pk->index].ctx_distance >>
						(dns_bidx * PMR_MAX_CTX_BITS)) & PMR_MAX_NCTX;
		}
		else if (nbhd->jump_ahead_countdown[_idx] > 0) {
			jump_ahead_countdown = nbhd->jump_ahead_countdown[_idx] - 1;
		}
#endif

		register dbl_t _tail = {.dbl = PRODUCT(nbhd->alphas[_idx], 
					nbhd->kmer_trans_prob[_idx][dns_bidx])};

		/* (greater or equal) ge_mask is set to 0, 
		 * if max_summand < _summand */
		uint64_t ge_mask = (uint64_t) (isnan(max_tail.dbl) || 
				(max_tail.dbl < _tail.dbl)) - 1UL;
		uint64_t inv_ge_mask = ~ge_mask;

		max_tail.u64 = (ge_mask & max_tail.u64) | 
			(inv_ge_mask & _tail.u64);
		argmax_tail = (ge_mask & argmax_tail) | (inv_ge_mask & _idx);

		++n;
	}

#if PMR_JUMP
	*countdown = jump_ahead_countdown;
#endif

	*delta = PRODUCT(max_tail.dbl, emit_dens);
	return argmax_tail;	
}

int hmmsbg_load_preconstructed_nbhd(data *d, read_t *read, 
		kbits_t obs_1st_kid, kbits_t obs_n_flag, const char *conv_q, 
		state_nbhd_t *first_nbhd, int kmer_len, int compute_emit_dens, 
		dict *kdict, mempool_t *mp, double *uniq_tp)
{
	size_t t_nbhd_size = 1;
	if (d->preconstructed_nbhds == NULL ||
			d->preconstructed_nbhds[read->id].nbhd == NULL) {
		// will use the observed kmer
		char *pmem = (char *) mempool_nalloc(mp, hmmsbg_nbhd_alloc_size(t_nbhd_size), 16);
		hmm_setup_nbhd_ptrs(first_nbhd, t_nbhd_size, pmem);

		kmer_t *nb_k = (kmer_t *) &(dict_find(d->kmer_dict, 
					kbits_cast_to_ptr(obs_1st_kid))->value);
		first_nbhd->states_sorted_sfx[0] = obs_1st_kid;
		first_nbhd->kmer_ptrs[0] = nb_k; 

		bitarray_set_pwr2(first_nbhd->ba_hamming_dist, 0, 0, 3);
		bitarray_set_pwr2(first_nbhd->ba_kmer_trans_flag, 0, 
				nb_k->trans_flag, 2);

		return 1;
	}

	int tmin_shifted = d->read_start_positions[read->id] > 0;
	kmer_nbhd_t knbhd = d->preconstructed_nbhds[read->id];

	if (tmin_shifted || knbhd.n == 0) {
		// use observed kmer
		knbhd.n = 1;
		knbhd.nbhd[0] = (kmer_t *) &(dict_find(d->kmer_dict, 
					kbits_cast_to_ptr(obs_1st_kid))->value);
	}

	t_nbhd_size = knbhd.n;
	char *pmem = (char *) mempool_nalloc(mp, hmmsbg_nbhd_alloc_size(t_nbhd_size), 16);
	hmmsbg_setup_nbhd_ptrs(first_nbhd, t_nbhd_size, pmem);

	first_nbhd->size = t_nbhd_size;
	for (int i = 0; i < t_nbhd_size; i++) {
		register kmer_t *nb_k = knbhd.nbhd[i];
		if (compute_emit_dens) {
			double emit_prob = hmm_compute_kmer_emission(d, obs_1st_kid, 
					obs_n_flag, nb_k, kmer_len, conv_q);
			double *init_p = hmmsbg_initial_p(d, nb_k);
			if (nb_k->has_multi_copies && nb_k->label == 0) {
				first_nbhd->alphas[i] = PRODUCT(emit_prob, 
						init_p[d->nonuniq_kmers[nb_k->index].n_ctx + 1]); 
			}
			else
				first_nbhd->alphas[i] = PRODUCT(emit_prob, *init_p); 
		}
		else {
			first_nbhd->alphas[i] = NAN;
		}

		kbits_t kid = kbits_id_of(nb_k);

		first_nbhd->states_sorted_sfx[i] = kid;
		first_nbhd->kmer_ptrs[i] = nb_k; 
		first_nbhd->kmer_trans_prob[i] = hmmsbg_locate_trans_p(
				d, nb_k, uniq_tp);

		if (nb_k->has_multi_copies) {
			kmer_ext_t kext = d->nonuniq_kmers[nb_k->index];

			if (nb_k->label == 0) {
				first_nbhd->kmer_trans_prob[i] += nb_k->n_trans * kext.n_ctx;
				first_nbhd->kmer_ups_ctx[i] = kext.n_ctx;
			}
			else {
				if (nb_k->n_trans == 1) {
					uint64_t kups_ctx = hmmsbg_determine_ctx(d, read, nb_k, 
							d->read_start_positions[read->id]);
					// if kups_ctx == n_ctx, it suggests ambiguous context
					first_nbhd->kmer_ups_ctx[i] = kups_ctx;
					//first_nbhd->kmer_trans_prob[i] += nb_k->n_trans * kups_ctx;

					// quick renormalization
					first_nbhd->kmer_trans_prob[i] = uniq_tp;
				}
				else {
					uint64_t kups_ctx = 0;
					// label-2 (implied by else) and non-unique transitions
					double *ctx_tp = (double *) mempool_alloc(mp, 
							nb_k->n_trans * sizeof(double));

					kups_ctx = hmmsbg_load_ctx_trans_p(d, read, nb_k, 
							d->read_start_positions[read->id], ctx_tp);
					renormalize_trans_p(nb_k->n_trans, ctx_tp);

					first_nbhd->kmer_ups_ctx[i] = kups_ctx;
					first_nbhd->kmer_trans_prob[i] = ctx_tp;
				}
			}
		}

		bitarray_set_pwr2(first_nbhd->ba_hamming_dist, i, 
				kmer_n_hamming_dist(obs_1st_kid, kid, obs_n_flag), 3);
		bitarray_set_pwr2(first_nbhd->ba_kmer_trans_flag, i, 
				nb_k->trans_flag, 2);
	}

	return t_nbhd_size;
}

void hmmsbg_reciprocal_exp_trans(data *d, kmer_t *pk, double *recip_exp_trans)
{

	int kmer_len = d->opt->kmer_length;
	int shift = (kmer_len - 1) << 1;
	kbits_t k_trans_packed = trans_flag2packed(pk->trans_flag);
	kbits_t n_trans = k_trans_packed & 7;
	int i = 0;

	kbits_t kid = kbits_id_of(pk);

	for (k_trans_packed >>= 3; i < n_trans; i++, k_trans_packed >>= 2) {
		kbits_t dns_base = k_trans_packed & 3;
		kbits_t rc_dns_base = complementary_base(kid & 3);

		kbits_t dns_kid = kbits_suffix(kid) | (dns_base << shift);
		kbits_t rc_dns_kid = kmer_reverse_complement(dns_kid, kmer_len);
		kmer_t *rc_pk = (kmer_t *) &(dict_find(d->kmer_dict, 
				kbits_cast_to_ptr(rc_dns_kid))->value);

		if ((rc_pk->trans_flag & (1 << rc_dns_base)) == 0) {
			/* no such transition on reverse complement kmer */
			recip_exp_trans[i] = 0.0;
		}
		else {
			int rc_nbidx = base_to_transidx(rc_pk->trans_flag, rc_dns_base);
			double *rc_exp_trans_p = hmmsbg_locate_exp_trans(d, rc_pk);
			recip_exp_trans[i] = rc_exp_trans_p[rc_nbidx];
		}
	}
}


void hmmsbg_multictx_lbl2_exp_trans(data *d, kmer_t *pk, int c, double *recip_exp_trans)
{
	int kmer_len = d->opt->kmer_length;
	int shift = (kmer_len - 1) << 1;
	kbits_t kid = kbits_id_of(pk);
	kbits_t tpacked = trans_flag2packed(pk->trans_flag) >> 3;
	kbits_t b1_rc = complementary_base(kid & 3);
	for (int i = 0; i < pk->n_trans; ++i, tpacked >>= 2) {
		kbits_t _k2 = kmer_reverse_complement(
				(kid >> 2) | ((tpacked & 3) << shift), kmer_len);
		kmer_t *pk2 = (kmer_t *) &(dict_find(d->kmer_dict,
					kbits_cast_to_ptr(_k2))->value);
		kmer_ext_t *pkext = &d->nonuniq_kmers[pk2->index];
		double *pk2_exp_trans = hmmsbg_locate_exp_trans(d, pk2);

		if ((pk2->trans_flag & (1 << b1_rc)) == 0) {
			recip_exp_trans[i] = 0.0;
			continue;
		}

		uint32_t tidx = base_to_transidx(pk2->trans_flag, (uint32_t) b1_rc);
		if (pk2->n_trans == 1) {
			recip_exp_trans[i] = pk2_exp_trans[c];
		}
		else {
			if (pk2->has_multi_copies) {
				int _toff = 0;
				for (int b = 0; b < tidx; ++b) 
					_toff += (((pkext->n_ctx >> (b * PMR_MAX_CTX_BITS)) & PMR_MAX_NCTX) + 1);
				recip_exp_trans[i] = pk2_exp_trans[_toff + c];
			}
			else {
				recip_exp_trans[i] = pk2_exp_trans[tidx];
			}
		}
	}
}

double hmmsbg_update_lbl2_trans_p(data *d, kmer_t *pk, double *_exp_cnt_lbl1)
{
	//kmer_ext_t *pkext = &d->nonuniq_kmers[pk->index];
	double _uniq_tp = 0.0;

	kbits_t kid = kbits_id_of(pk);
	uint32_t b1_rc = (uint32_t) complementary_base(kid & 3);
	int kmer_len = d->opt->kmer_length;
	int shift = (kmer_len - 1) << 1;

	double k_cexp_cnt = 0.0, rev_cexp_cnt = 0.0;

	uint32_t n_ctx = (pk->has_multi_copies) ? 
		d->nonuniq_kmers[pk->index].n_ctx : 0;

	double *cur_tp = hmmsbg_locate_trans_p(d, pk, &_uniq_tp);
	double mu = *hmmsbg_initial_p(d, pk);

	int par_off = 0;
	kbits_t tpacked = trans_flag2packed(pk->trans_flag) >> 3;
	double sum_tp = 0.0;
	for (int i = 0; i < pk->n_trans; ++i, tpacked >>= 2) {
		kbits_t base = tpacked & 3;
		/*
		uint32_t n_ctx_i = pk->has_multi_copies ?
			(n_ctx >> (i * 4)) & 0xf : 1;
		*/
		uint32_t n_ctx_i = (n_ctx >> (i * PMR_MAX_CTX_BITS)) & PMR_MAX_NCTX;

		// single-context
		// if (n_ctx_i == 0) n_ctx_i = 1;

		kbits_t kid_lbl1 = kmer_reverse_complement(
				(kid >> 2) | (base << shift), kmer_len);
		kmer_t *pk_lbl1 = (kmer_t *) &(dict_find(d->kmer_dict, 
					kbits_cast_to_ptr(kid_lbl1))->value);
		int tidx = base_to_transidx(pk_lbl1->trans_flag, b1_rc);

		double *pmu_lbl1 = hmmsbg_initial_p(d, pk_lbl1);
		if (pk_lbl1->trans_flag & (1 << b1_rc)) {
			double *_tp_lbl1 = hmmsbg_locate_trans_p(d, pk_lbl1, &_uniq_tp);

			// for (int c = 0; c < n_ctx_i; ++c) {
			for (int c = 0; c <= n_ctx_i; ++c) {
				double _tp = _tp_lbl1[c * pk_lbl1->n_trans + tidx] - mu + pmu_lbl1[c];

				k_cexp_cnt += *(hmmsbg_locate_exp_trans(d, pk_lbl1) + 
					(c * pk_lbl1->n_trans) + tidx);
				
				// rev_cexp_cnt += EXP(_tp) *
				/*
				if (c == n_ctx_i) {
		 			k_cexp_cnt += *(hmmsbg_locate_exp_trans(d, pk_lbl1) + 
						(n_ctx_i * pk_lbl1->n_trans) + tidx);
					rev_cexp_cnt += EXP(_tp) * 
						_exp_cnt_lbl1[hmmsbg_initial_index(d, pk_lbl1)];
				}
				*/

				sum_tp += EXP(_tp);
				if (pk->n_trans > 1 || pk->has_multi_copies) 
					cur_tp[par_off + c] = EXP(_tp); 
			}
		}

		par_off += (n_ctx_i + 1);
		// par_off += (n_ctx_i + pk_lbl1->has_multi_copies);
	}

	if (pk->n_trans == 1 && pk->has_multi_copies == 0) {
		if (sum_tp == 0.0) {
			pk->n_trans = 0;
			pk->trans_flag = 0;
			pk->index = 0;
		}
		return k_cexp_cnt;
	}

	// renormalization
	tpacked = trans_flag2packed(pk->trans_flag) >> 3;
	par_off = 0;
	for (int i = 0; i < pk->n_trans; ++i, tpacked >>= 2) {
		/*
		uint32_t n_ctx_i = 0, n_params = 0;
		if (pk->has_multi_copies) {
			n_ctx_i = (n_ctx >> (i * 4)) & 0xf;
			n_params = n_ctx_i + 1;
			if (n_ctx_i == 0) n_ctx_i = 1;
		}
		else {
			n_ctx_i = 1;
			n_params = 1;
		}
		*/
		uint32_t n_ctx_i = (n_ctx >> (i * PMR_MAX_CTX_BITS)) & PMR_MAX_NCTX;

		//for (int c = 0; c < n_ctx_i; ++c) {
		for (int c = 0; c <= n_ctx_i; ++c) {
			cur_tp[par_off + c] = LOG(cur_tp[par_off + c] / sum_tp);
		}

		// par_off += n_params; 
		par_off += (n_ctx_i + 1);
	}

	return k_cexp_cnt;
}

/* Record the observed contexts for label-1, non-unique kmers */
void hmmsbg_preindex_contexts(data *d, read_t *read, void *fdata)
{
	const int kmer_len = d->opt->kmer_length;
	const int shift = (kmer_len - 1) << 1;
	int tmax = read->length - kmer_len;
	char *rseq = read->sequence;

	std::unordered_map<kbits_t, std::unordered_set<std::string>> &ups_ctx =
		*(std::unordered_map<kbits_t, std::unordered_set<std::string>> *) fdata;

	int tmin = d->opt->max_backtrack_dist >> 1;

	kbits_t obs_kid = kmer_seq_to_id(&rseq[tmin], kmer_len);
	kmer_t *pk = kptr_locate(d, obs_kid);
	for (int t = tmin; t <= tmax; ++t) {
		int tnext_kp = t + kmer_len;
		int ctxlen = MIN(d->opt->max_backtrack_dist, t);

		if (pk->label == 0 && pk->n_trans > 1) {
			std::string _ctx(&read->sequence[t - ctxlen], ctxlen);
			ups_ctx[obs_kid].insert(std::move(_ctx));
		}

		if (t < tmax) {
			obs_kid = (obs_kid >> 2) | 
				((kbits_t) base_to_bits(rseq[tnext_kp]) << shift);
			pk = kptr_locate(d, obs_kid);
		}
	}
}

/* 
 * Before calling this function,
 * kmer_t.uniq_trans is the indicator used to determine whether the kmer
 * associated variables are located by kmer_t.index or
 * d->nonuniq_kmers[kmer_t.index].poff;
 *
 * After calling this function, one also needs to examine
 * kmer_t.has_multi_copies.
 */

int hmmsbg_identify_contexts(data *d)
{
	const int kmer_len = d->opt->kmer_length;
	const int shift = (kmer_len - 1) << 1;
	const double thresh = d->opt->penalty_eta / log1p(1/d->opt->penalty_gamma);
	int64_t n = d->kmer_dict->sizemask + 1;

	std::vector<kmer_ext_t> _multictx_kmers;
	std::vector<kmer_t *> _multictx_ptrs;
	std::vector<double> _multictx_trans_p;
	std::vector<kbits_t> _multictx_contexts;
	std::vector<double> _multictx_init_p;

	_multictx_kmers.reserve(d->n_nonuniq_kmers);
	_multictx_ptrs.reserve(d->n_nonuniq_kmers);
	_multictx_trans_p.reserve(d->n_nonuniq_trans * 
			d->opt->n_max_kmer_contexts);
	_multictx_contexts.reserve(d->n_nonuniq_kmers *
			d->opt->n_max_kmer_contexts);
	_multictx_init_p.reserve(d->n_nonuniq_kmers *
			d->opt->n_max_kmer_contexts);

	// idx_multi_kmer corresponds to a kmer with multiple outgoing transitions,
	// or with multiple upstream contexts.
	uint32_t idx_multi_kmer = 0;
	uint32_t idx_params = 0, idx_contexts = 0;
	
	// number of total transitions for all single context kmers,
	// this does not account for the number of unique transitions.
	int n_sctx_nonuniq_trans = 0, n_sctx_nonuniq_kmers = 0;
	int n_sctx_nonzero_kmers = 0;

	std::unordered_map<kbits_t, 
		std::unordered_set<std::string>> k_ctx_index;
	/* pre-indexing all upstream contexts for label-1, non-unique kmers */
	io_iterate_reads_fastq(d, 1, d->opt->is_tty, hmmsbg_preindex_contexts, NULL,
			(void *) &k_ctx_index);

	/* first pass, identify all kmers with multiple contexts, 
	 * and mark the corresponding kmers on the opposite strand */
	for (int64_t i = 0; i < n; ++i) {
		dictentry *de = &d->kmer_dict->table[i];
		if (de->value == NULL) continue;

		kmer_t *pk = (kmer_t *) &(de->value);

		n_sctx_nonzero_kmers += (pk->n_trans > 0);
		if (pk->n_trans > 1) {
			n_sctx_nonuniq_trans += pk->n_trans;
			++n_sctx_nonuniq_kmers;
		}

		if (pk->label == 0 && pk->n_trans > 1) {
			kbits_t kid = kbits_id_of(pk);

			uint32_t b1_rc = (uint32_t) complementary_base(kid & 3);

			double kcvg = d->total_kmer_exp_counts *
				EXP(d->avg_initial_p[kmer_initial_index(pk, d->uniq_kmer_init_off)]);

			std::vector<kbits_t> ups_ctx;
			auto &ctx_idx = k_ctx_index[kid];	

			int n_max_ctxs = std::min((int) ctx_idx.size(),
					std::min(d->opt->n_max_kmer_contexts, 
						(int) (kcvg / 2.0 / thresh)));
			// if (n_max_ctxs >= 127) n_max_ctxs = 127;
			
			if (n_max_ctxs <= 1) continue;

			uint32_t ctx_dist = hmmsbg_backtrack(d, kid, 
					n_max_ctxs,
					d->opt->max_backtrack_dist, ups_ctx,
					ctx_idx);

			// 16 bits is enough here, since we are not multiplexing the 16
			// bits (into 4 groups).
			uint32_t n_ctx = ctx_dist & UINT16_MAX;
			ctx_dist >>= 16;

			// unique context within the backtrack distance.
			if (n_ctx <= 1) continue;

			n_sctx_nonuniq_trans -= pk->n_trans;
			--n_sctx_nonzero_kmers;
			--n_sctx_nonuniq_kmers;

			pk->has_multi_copies = 1;
			kmer_ext_t kext = {d->nonuniq_kmers[pk->index].poff, 
				idx_contexts, n_ctx, ctx_dist};

			_multictx_kmers.push_back(kext);
			_multictx_ptrs.push_back(pk);

			std::copy(ups_ctx.begin(), ups_ctx.end(), 
					std::back_inserter(_multictx_contexts));

			// pk->index = idx_multi_kmer++;

			// idx_params += pk->n_trans * (n_ctx + 1);
			/*
			double *cur_tp = d->transition_p + d->nonuniq_kmers[pk->index].poff;
			for (int j = 0; j <= n_ctx; ++j)
				std::copy(cur_tp, cur_tp + pk->n_trans, 
						std::back_inserter(_multictx_trans_p));
			*/

			// process the label 2 kmers linked with current kmer
			kbits_t tpacked = trans_flag2packed(pk->trans_flag);
			int n = tpacked & 7;
			for (tpacked >>= 3; n > 0; --n, tpacked >>= 2) {
				kbits_t klbl2 = kmer_reverse_complement(
						(kid >> 2) | ((tpacked & 3) << shift), kmer_len);

				kmer_t *pk_lbl2 = (kmer_t *) &(dict_find(d->kmer_dict, 
							kbits_cast_to_ptr(klbl2))->value);

				if (pk_lbl2->n_trans == 1) {
					// simpler case, the label-2 kmer has a unique transition.
					kmer_ext_t kext = {0, idx_contexts, n_ctx, ctx_dist};
					_multictx_kmers.push_back(kext);
					_multictx_ptrs.push_back(pk_lbl2);

					--n_sctx_nonzero_kmers;
					pk_lbl2->has_multi_copies = 1;
					// pk_lbl2->index = idx_multi_kmer++;

					// idx_params += (n_ctx + 1);
				}
				else if (pk_lbl2->n_trans > 1) {
					// label-2 kmer also has multiple transitions
					if (pk_lbl2->trans_flag & (1 << b1_rc)) {
						int idx = base_to_transidx(pk_lbl2->trans_flag, b1_rc);
						double *_tp = NULL;
						// first encounter
						if (pk_lbl2->has_multi_copies == 0) {
							// kext.ctxoff is not used for label-2 kmers, will
							// use it to store the old pk_lbl2->index.
							kmer_ext_t kext = {
								d->nonuniq_kmers[pk_lbl2->index].poff, 0,
								n_ctx << (idx * PMR_MAX_CTX_BITS), 
								ctx_dist << (idx * PMR_MAX_CTX_BITS)
							};

							// rewrite poff to store the index of kext in
							// _multictx_kmers
							d->nonuniq_kmers[pk_lbl2->index].poff = _multictx_kmers.size();

							_multictx_kmers.push_back(kext);
							_multictx_ptrs.push_back(pk_lbl2);

							/*
							_tp = d->transition_p + 
								d->nonuniq_kmers[pk_lbl2->index].poff;
							*/

							n_sctx_nonuniq_trans -= pk_lbl2->n_trans;
							--n_sctx_nonuniq_kmers;
							--n_sctx_nonzero_kmers;

							pk_lbl2->has_multi_copies = 1;

							// pk_lbl2->index = idx_multi_kmer++;
						}
						else {
							// update the existing contexts
							kmer_ext_t &kext = _multictx_kmers[d->nonuniq_kmers[pk_lbl2->index].poff];
							kext.n_ctx |= (n_ctx << (idx * PMR_MAX_CTX_BITS));
							kext.ctx_distance |= (ctx_dist << (idx * PMR_MAX_CTX_BITS));

							/*
							_tp = d->transition_p + 
								d->nonuniq_kmers[kext.ctxoff].poff;
							*/
						}

						/*
						for (int b = 0; b <= n_ctx; ++b)
							_multictx_trans_p.push_back(_tp[idx]);
						idx_params += (n_ctx + 1);
						*/
					}
				}
			}

			idx_contexts += n_ctx;
		}
	}

	int n_multictx_kmers = _multictx_kmers.size();
	int _uoff = n_sctx_nonuniq_kmers;

	double *_avg_init_p = new double[d->n_nonzero_kmers + 1];
#if PMR_EXTEND_READS
	double *_prob_kmer_ext = new double[d->n_nonzero_kmers + 1];
#endif
	uint32_t idx_uniq_kmer = 0; 
	/* second pass, recompute the parameter indices etc */
	for (int64_t i = 0; i < n; ++i) {
		dictentry *de = &d->kmer_dict->table[i];
		if (de->value == NULL) continue;

		kmer_t *pk = (kmer_t *) &(de->value);

		if (pk->has_multi_copies || pk->n_trans == 0) continue;

		int _poff = kmer_initial_index(pk, d->uniq_kmer_init_off);
		double _init_p = d->avg_initial_p[_poff];
#if PMR_EXTEND_READS
		double _p_ext = d->prob_kmer_extension[_poff];
#endif
		if (pk->n_trans == 1) {
			pk->index = idx_uniq_kmer++;
			_avg_init_p[_uoff + pk->index + 1] = _init_p;
#if PMR_EXTEND_READS
			_prob_kmer_ext[_uoff + pk->index + 1] = _p_ext;
#endif
		}
		else if (pk->n_trans > 1) {
			double *cur_tp = d->transition_p + d->nonuniq_kmers[pk->index].poff;

			pk->index = idx_multi_kmer++;
			_avg_init_p[pk->index + 1] = _init_p;
#if PMR_EXTEND_READS
			_prob_kmer_ext[pk->index + 1] = _p_ext;
#endif
			
			kmer_ext_t kext = {idx_params};
			_multictx_kmers.push_back(kext);
			std::copy(cur_tp, cur_tp + pk->n_trans, 
					std::back_inserter(_multictx_trans_p));
			idx_params += pk->n_trans;
		}
	}

	int _moff = n_sctx_nonzero_kmers - idx_multi_kmer;
	int idx_multi_ctx = 0, idx_mctx_init_p = 0;
	for (auto it = _multictx_ptrs.begin(); it != _multictx_ptrs.end(); 
				++it, ++idx_multi_ctx) {
		kmer_t *pk = *it;
		kmer_ext_t &kext = _multictx_kmers[idx_multi_ctx];

		int _poff = kmer_initial_index(pk, d->uniq_kmer_init_off);
		double _init_p = d->avg_initial_p[_poff];
		pk->index = idx_multi_kmer++; 

		_avg_init_p[_moff + pk->index + 1] = (pk->label == 0) ?
			(double) idx_mctx_init_p : _init_p; 
#if PMR_EXTEND_READS
		_prob_kmer_ext[_moff + pk->index + 1] = d->prob_kmer_extension[_poff];
#endif

		if (pk->label == 0) {
			double *_tp = d->transition_p + kext.poff;

			for (int j = 0; j <= kext.n_ctx; ++j) {
				std::copy(_tp, _tp + pk->n_trans, 
						std::back_inserter(_multictx_trans_p));

				// split pk into n_ctx copies
				_multictx_init_p.push_back(_init_p - LOG(kext.n_ctx + 1));
			}

			_multictx_init_p.push_back(_init_p);

			kext.poff = idx_params;
			idx_params += pk->n_trans * (kext.n_ctx + 1);
			idx_mctx_init_p += (kext.n_ctx + 2);
		}
		else if (pk->unique_trans) {
			for (int j = 0; j < kext.n_ctx; ++j)
				//_multictx_trans_p.push_back(0.0);
				_multictx_trans_p.push_back(-LOG(kext.n_ctx));

			_multictx_trans_p.push_back(0.0);

			kext.poff = idx_params;
			idx_params += (kext.n_ctx + 1);
		}
		else {
			double *_tp = d->transition_p + kext.poff;
			kext.poff = idx_params;

			for (int b = 0; b < pk->n_trans; ++b) {
				uint32_t n_ctx_b = (kext.n_ctx >> (b * PMR_MAX_CTX_BITS)) & PMR_MAX_NCTX;
				double _tpb = (n_ctx_b > 0) ? 
					_tp[b] - LOG(n_ctx_b) : _tp[b];

				// equivalent to splitting the trans. prob. into multiple
				// copies, and immediately coalesce them.
				/*
				for (int j = 0; j <= n_ctx_b; ++j)
					_multictx_trans_p.push_back(_tp[b]);
				*/
				for (int j = 0; j < n_ctx_b; ++j)
					_multictx_trans_p.push_back(_tpb);

				// for ambiguous context, use the sum for individual copies.
				_multictx_trans_p.push_back(_tp[b]);

				idx_params += (n_ctx_b + 1);
			}
		}
	}

	// overwrite old data
	d->nonuniq_kmers = (kmer_ext_t *) realloc(d->nonuniq_kmers, 
			_multictx_kmers.size() * sizeof(*d->nonuniq_kmers));
	memcpy(d->nonuniq_kmers, &_multictx_kmers[n_multictx_kmers], 
			(_multictx_kmers.size() - n_multictx_kmers) * sizeof(*d->nonuniq_kmers));
	memcpy(&d->nonuniq_kmers[_multictx_kmers.size() - n_multictx_kmers], 
			&_multictx_kmers[0], n_multictx_kmers * sizeof(*d->nonuniq_kmers));

	d->transition_p = (double *) realloc(d->transition_p,
			_multictx_trans_p.size() * sizeof(*d->transition_p));
	memcpy(d->transition_p, &_multictx_trans_p[0],
			_multictx_trans_p.size() * sizeof(*d->transition_p));

	size_t n_tp = _multictx_trans_p.size();
	d->exp_trans_p = (double *) realloc(d->exp_trans_p,
			(n_tp + idx_uniq_kmer) * sizeof(*d->exp_trans_p));
	memset(d->exp_trans_p, 0, 
			(n_tp + idx_uniq_kmer) * sizeof(*d->exp_trans_p));

	d->kmer_ups_contexts = (kbits_t *) malloc(
			_multictx_contexts.size() * sizeof(*d->kmer_ups_contexts));
	memcpy(d->kmer_ups_contexts, &_multictx_contexts[0], 
			_multictx_contexts.size() * sizeof(*d->kmer_ups_contexts));

	memcpy(d->avg_initial_p, _avg_init_p,
			(d->n_nonzero_kmers + 1) * sizeof(*d->avg_initial_p));
#if PMR_EXTEND_READS
	memcpy(d->prob_kmer_extension, _prob_kmer_ext,
			(d->n_nonzero_kmers + 1) * sizeof(*d->prob_kmer_extension));
#endif

	d->multictx_initial_p = new double[_multictx_init_p.size()];
	memcpy(d->multictx_initial_p, &_multictx_init_p[0],
			_multictx_init_p.size() * sizeof(*d->multictx_initial_p));

	d->uniq_kmer_init_off = _uoff;
	d->multictx_kmer_init_off = _moff;
	d->uniq_kmer_param_off = n_tp; 
	d->n_total_trans = n_tp + idx_uniq_kmer;

	delete [] _avg_init_p;
#if PMR_EXTEND_READS
	delete [] _prob_kmer_ext;
#endif
}

/** Perform kmer backtracking on the de Brujin graph,
 * returns an uint32_t whose lowest 16-bits are the # of contexts,
 * and the highest 16-bits are the distance to the context */
uint32_t hmmsbg_backtrack(data *d, kbits_t kid, int max_paths, int max_len,
		std::vector<kbits_t> &vctx, 
		std::unordered_set<std::string> &k_ups_ctxs)
{
	int actg_trans_count[4];
	int cum_actg_trans_count[4];
	int total_hmm_trans = 0;

	const int kmer_len = d->opt->kmer_length;
	const kbits_t klen_flag = (1UL << (kmer_len << 1)) - 1;
	const int shift = (kmer_len - 1) << 1;

	char tmp_kmer[kmer_len + 1];

	mempool_t *mp = mempool_create(MEMPOOL_DEFAULT_SIZE);

	state_nbhd_t *T_nbhds = (state_nbhd_t *) mempool_nalloc(mp,
			sizeof(*T_nbhds) * (max_len + 1), 16);

	kid = kmer_reverse_complement(kid, kmer_len);

	state_nbhd_t *first_nbhd = &T_nbhds[0];
	first_nbhd->size = 1;
	char *pmem = (char *) mempool_nalloc(mp, nbhd_alloc_size(1), 16);
	hmm_setup_nbhd_ptrs(first_nbhd, 1, pmem);

	kmer_t *nb_k = (kmer_t *) &(dict_find(d->kmer_dict, 
				kbits_cast_to_ptr(kid))->value);
	first_nbhd->states_sorted_sfx[0] = kid;
	first_nbhd->kmer_ptrs[0] = nb_k; 

	bitarray_set_pwr2(first_nbhd->ba_hamming_dist, 0, 0x7f, 3);
	bitarray_set_pwr2(first_nbhd->ba_kmer_trans_flag, 0, nb_k->trans_flag, 2);

	ksubstr_t dmax = 0xff;

	//int uniq_ups_context = 1;
	int argmax_ctx_size = -1, max_ctx_size = 0;
	/* initial state distribution */
	int t = 0;
	for (; t < max_len; ++t) {
		state_nbhd_t *curr_nbhd = &T_nbhds[t],
					 *next_nbhd = &T_nbhds[t+1];

		int t_nbhd_size = curr_nbhd->size;

		kbits_t *t_states_sorted_sfx = T_nbhds[t].states_sorted_sfx;

		memset(actg_trans_count, 0, sizeof(int) << 2);

		int n_distinct_suffixes = 0;
		int tnext_nbhd_size = 0;
		int tnext_valid_ctxs = 0;

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

		//uniq_ups_context &= (tnext_nbhd_size == 1);

		if (tnext_nbhd_size == 0) { 
			// use curr_nbhd as the contexts
			break;
		}

		/* allocate memory given tnext_nbhd_size */
		char *pmem = (char *) mempool_nalloc(mp,
				nbhd_alloc_size(tnext_nbhd_size), 16);

		next_nbhd->size = tnext_nbhd_size;
		hmm_setup_nbhd_ptrs(next_nbhd, tnext_nbhd_size, pmem);

		int kmer_idx = 0, index_pfx = 0;
		for (int i = 0; i < n_distinct_suffixes; i++) {

			kbits_t _repr_kid = t_states_sorted_sfx[kmer_idx];
			kbits_t common_sfx = kbits_suffix(_repr_kid);

			/* v 4 bits                   v 4 bits (higher)
			 *
			 * A  -                       - A
			 * C  |  -- k-1 substring --  | C
			 * T  |                       | T
			 * G  -                       - G
			 */
			kbits_t _tf = bitarray_get_pwr2(curr_nbhd->ba_pfx_sfx_flags, i, 3);
			kbits_t ups_trans_packed = trans_flag2packed(_tf & 15);
			/* which kmers/nucleotides current suffix transition into */
			kbits_t dns_trans_packed = trans_flag2packed(_tf >> 4);

			int n_common_sfx_kmers = ups_trans_packed & 7; 

			int j = dns_trans_packed & 7;
			for (dns_trans_packed >>= 3; j > 0; --j, dns_trans_packed >>= 2) {
				kbits_t dns_base = dns_trans_packed & 3;
				kbits_t tnext_kid = common_sfx | (dns_base << shift);

				/* vetting if the upstream context is valid or not */
				int ctx_len = ((t+1) < kmer_len) ? t + 1 : kmer_len;
				// kbits_t ctx_mask = (t < kmer_len) ? (1UL << t) - 1 : UINT64_MAX;
				// kbits_t tnext_ctx = kmer_reverse_complement(tnext_kid, kmer_len);

				std::string put_ctx(
						kmer_id_to_seq(kmer_reverse_complement(tnext_kid, kmer_len), 
							ctx_len, tmp_kmer));

				auto hit = std::find_if(
						k_ups_ctxs.begin(), k_ups_ctxs.end(), 
						[&](const std::string &uctx) {
							if (uctx.size() <= t) return false;
							return uctx.compare(uctx.size() - (t + 1), ctx_len, put_ctx) == 0;
						});
				bool ctx_valid = (hit != k_ups_ctxs.end());

				/* compute the position sorted by suffix */
				int index_sfx = cum_actg_trans_count[dns_base] - 
					actg_trans_count[dns_base];
				actg_trans_count[dns_base]--;

				kmer_t *tnext_kmer = (kmer_t *) &(dict_find(d->kmer_dict, 
						kbits_cast_to_ptr(tnext_kid))->value);
				next_nbhd->kmer_ptrs[index_sfx] = tnext_kmer;

				/* set up transition flag */
				if (ctx_valid) {
					++tnext_valid_ctxs;
					bitarray_set_pwr2(next_nbhd->ba_kmer_trans_flag,
							index_sfx, tnext_kmer->trans_flag, 2);
					// if the context is not valid, tnext_kid is a "chimeric"
					// kmer, hence there's no point to move forward.
				}
				bitarray_set_pwr2(next_nbhd->ba_hamming_dist, index_sfx, 
						0x7f, 3);

				//next_nbhd->suffixes_order[index_pfx] = index_sfx;
				next_nbhd->states_sorted_sfx[index_sfx] = tnext_kid;
				next_nbhd->alphas[index_sfx] = (double) ctx_valid;
				//index_pfx++;
			}

			kmer_idx += n_common_sfx_kmers;
		}

		if (tnext_valid_ctxs > max_ctx_size) {
			max_ctx_size = tnext_valid_ctxs;
			argmax_ctx_size = t;
		}

		if (tnext_valid_ctxs >= max_paths) break;
	}

	if (t == 0) {
		mempool_destroy(mp);
		return 0;
	}

	if (t == max_len && argmax_ctx_size >= 0)
		t = argmax_ctx_size + 1;

	state_nbhd_t *t_nbhd = &T_nbhds[t];

	/*
	fprintf(stderr, "Backtracking from %lu [%d steps] (%d states):\n", 
			kid, t, t_nbhd->size);
	*/

	if (t_nbhd->size <= 0) return 0;
	vctx.reserve(t_nbhd->size);

	// char ctx_kmer[kmer_len + 1];
	for (int i = 0; i < t_nbhd->size; ++i) {
		// skip "invalid" upstream contexts without support by reads.
		if (t_nbhd->alphas[i] == 0.0) continue;

		// if the number of steps backtracked is smaller than k, only record the
		// non-overlapped prefix.
		kbits_t ctx_kid = kmer_reverse_complement(t_nbhd->states_sorted_sfx[i], 
				kmer_len);

		// if the context is found within a kmer overlapping with the current
		// kmer, only store the non-overlapping prefix.
		if (t < kmer_len) {
			ctx_kid = kbits_prefix(ctx_kid, t);
		}

		/*
		double kcexp_cnt = EXP(d->avg_initial_p[kmer_initial_index(
					t_nbhd->kmer_ptrs[i], d->uniq_kmer_init_off)]) * d->total_kmer_exp_counts;

			kmer_id_to_seq(ctx_kid, t, ctx_kmer);
		}
		else {
			kmer_id_to_seq(ctx_kid, kmer_len, ctx_kmer);
		}

		fprintf(stderr, " -> [5'] %s (%g) [3']\n", ctx_kmer, kcexp_cnt);
		*/
		vctx.push_back(ctx_kid);
	}


	mempool_destroy(mp);

	return (vctx.size()) | (t << 16);
}

void HMMSBG(data *d)
{
	INFO_TEE(d->log_fp,
			"Splitting kmers into multiple copies using the de Bruijn graph...\n");

	// re-assign labels to kmers to optimize performance
	//hmmsbg_relabel_kmers(d);

	hmmsbg_identify_contexts(d);

	// hmm_init_model_params(d);

	// EM on the graph
	for (int iter = 0; iter < 4; ++iter) {
		d->loglik = 0.0;

		hmm_init_e_step(d);

		INFO_TEE(d->log_fp, "E-step / HMMSBG ... \n");
		io_iterate_reads_fastq(d, d->opt->n_threads, d->opt->is_tty, hmmsbg_e_step, NULL, NULL);

		if (d->opt->mstep_debug_console) {
			hmmsbg_debug_console(d);
		}
		hmmsbg_m_step(d);
	}

	hmmsbg_coalesce_lbl2_trans(d);

	// hmmsbg_debug_console(d);
	
	// Viterbi
	io_iterate_reads_fastq(d, d->opt->n_threads, d->opt->is_tty, hmmsbg_viterbi, NULL, NULL);
}

void hmmsbg_e_step(data *d, read_t *read, void *fdata)
{
	int actg_trans_count[4];
	int cum_actg_trans_count[4];

	const int kmer_len = d->opt->kmer_length;
	const int num_kmers = read->length - kmer_len + 1;
	const kbits_t klen_flag = (1UL << (kmer_len << 1)) - 1;
	const int shift = (kmer_len - 1) << 1;
	const int ctx_shift = (kmer_len - d->opt->base_emit_context) << 1;
	const int qmax = d->opt->max_qual_score + 1;
	const int qmax2 = qmax << 1;
	const int bwidth = d->opt->qscore_bin_width;
	double _uniq_tp = 0.0;
	const size_t max_approx_nbhd_size = d->opt->max_approx_nbhd_size;

	int tmin = d->read_start_positions[read->id];
	int tmax = read->length - kmer_len;
	if (tmin >= tmax) {
		return; 
	}

	mempool_t *mp = mempool_create(MEMPOOL_DEFAULT_SIZE);

	dict *kdict = d->kmer_dict;

	/* sequence */
	char *rseq = read->sequence;
	char *conv_q = convert_quality_scores(read->qscore, read->length, 
			d->opt->qual_score_offset, bwidth, mp);

	state_nbhd_t *T_nbhds = (state_nbhd_t *) mempool_nalloc(mp,
			sizeof(*T_nbhds) * num_kmers, 16);

	/* flag for N (not determined) base, first kmer only */
	kbits_t obs_n_flag = kmer_n_base_flag(rseq, kmer_len);
	kbits_t next_obs_n_flag = kmer_n_base_flag(rseq + 1, kmer_len);

	double tnext_max_alpha = NAN;
	int tnext_argmax_alpha = 0;

	/* convert string literal to numeric representation */
	kbits_t observed_kmers[BITS_TO_U64(read->length << 1)];
	read_seq_to_numeric(rseq, observed_kmers, read->length, kmer_len);


	/* set up first kmer neighborhood */
	kbits_t obs_kid = bitarray_get(observed_kmers, tmin << 1, kmer_len << 1);
	ksubstr_t obs_base = kmer_effective_base(rseq[tmin+kmer_len-1]);
	ksubstr_t obs_next_base = kmer_effective_base(rseq[tmin+kmer_len]);

	T_nbhds[tmin].size = hmmsbg_load_preconstructed_nbhd(
			d, read, obs_kid, 
			obs_n_flag, conv_q, &T_nbhds[tmin], kmer_len, 1, kdict, mp, 
			const_cast<double *>(&_uniq_tp));

	int wsize = d->opt->hd_window_size;
	uint8_t *win_dmax = &d->read_dmax[read->id * d->n_hd_windows];

	/* initial state distribution */
	const size_t sz_base_emit_p_off = kmer_len * 20;
	const size_t sz_qual_emit_p_off = kmer_len * qmax2;
	int t = tmin;
	for (; t < tmax; t++) {
		int tnext_kp = t + kmer_len;
		int w_emit = t / d->opt->emit_window_size;

#if !PMR_USE_PRIMITIVE_EMISSION
		double *tnext_base_emit_p = &d->base_emission_p[w_emit * 
			d->base_ctx_radix + sz_base_emit_p_off];
		double *tnext_qual_emit_p = &d->qual_emission_p[w_emit *
			d->qual_ctx_radix + sz_qual_emit_p_off];
#endif

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

		ksubstr_t dmax = win_dmax[t / wsize];

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
				hmmsbg_nbhd_alloc_size(tnext_nbhd_size), 16);

		next_nbhd->size = tnext_nbhd_size;
		hmmsbg_setup_nbhd_ptrs(next_nbhd, tnext_nbhd_size, pmem);

		int kmer_idx = 0, index_pfx = 0;

		tnext_max_alpha = NAN;
		tnext_argmax_alpha = 0;
		int tnext_nonzero_states = 0;

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

#if PMR_USE_PRIMITIVE_EMISSION
				double emit_dens = (dns_base == obs_next_base) ? 
					LOG(1.0-d->base_error_rate) : LOG(d->base_error_rate/3);
#else
				kbits_t emit_renorm_flag = (sfx_hd == dmax) ? 
					((1U << obs_next_base) & 0xF) : 15; //obs_next_kmer->emit_norm_flag;

				/* compute the emission probability */
				int bctx = calc_bctx(tnext_kid, sfx_hd, shift, tnext_kmer->single_stranded);
				// int qctx = calc_qctx(tnext_kid, obs_next_kid, conv_q + t + 1, kmer_len);
				int qctx = (tnext_kid >> shift) == (obs_next_kid >> shift);

				double emit_dens = hmm_compute_emit_dens(
						dns_base, obs_next_base, 
						(int) conv_q[tnext_kp], 
						emit_renorm_flag,
						tnext_base_emit_p + bctx * 5,
						tnext_qual_emit_p + qctx * qmax);
#endif

				int64_t countdown = 0;
				double alpha = hmmsbg_compute_alpha(d,
						ups_trans_packed, kmer_idx,
						curr_nbhd, &countdown, emit_dens, dns_base);

#if PMR_EXTEND_READS
				if (t == tmax - 1) {
					// for the last kmer, incorporate the probability of read
					// extension.
					double pext = d->prob_kmer_extension[
							hmmsbg_initial_index(d, tnext_kmer)];
					alpha = PRODUCT(alpha, pext);

					// initialize the backward algorithm recursion properly.
					next_nbhd->betas[index_sfx] = pext;
				}
#endif

				if (!isnan(alpha)) tnext_nonzero_states++;

				update_max_quantity(alpha, index_sfx,
						tnext_max_alpha, tnext_argmax_alpha);

				next_nbhd->kmer_ptrs[index_sfx] = tnext_kmer;
				next_nbhd->kmer_trans_prob[index_sfx] = 
					hmmsbg_locate_trans_p(d, tnext_kmer, &_uniq_tp);

#if PMR_JUMP
				next_nbhd->jump_ahead_countdown[index_sfx] = countdown;

				if (countdown > 0) {
					kbits_t jump_step = kmer_effective_base(rseq[t + kmer_len + 1]);
					// set the transition probability accordingly
					double *ctx_tp = (double *) mempool_alloc(mp, 
							tnext_kmer->n_trans * sizeof(double));

					kbits_t tpacked = trans_flag2packed(tnext_kmer->trans_flag) >> 3;
					for (int b = 0; b < tnext_kmer->n_trans; ++b, tpacked >>= 2) {
						kbits_t base = tpacked & 3;
						if (base == jump_step) ctx_tp[b] = 0.0;
						else ctx_tp[b] = NAN;
					}

					next_nbhd->kmer_trans_prob[index_sfx] = ctx_tp;
				}
				else if (tnext_kmer->has_multi_copies) {
#else
				if (tnext_kmer->has_multi_copies) {
#endif
					if (tnext_kmer->label == 0 || tnext_kmer->n_trans == 1) {
						uint64_t kups_ctx = hmmsbg_determine_ctx(d, read, tnext_kmer, t+1);
						// if kups_ctx == n_ctx, it suggests ambiguous context
						next_nbhd->kmer_ups_ctx[index_sfx] = kups_ctx;
						next_nbhd->kmer_trans_prob[index_sfx] += 
							tnext_kmer->n_trans * kups_ctx;

						// quick renormalization if the transition probability
						// is NOT zero.
						if (tnext_kmer->n_trans == 1) // && 
								//!IS_ZERO(next_nbhd->kmer_trans_prob[index_sfx][0]))
							next_nbhd->kmer_trans_prob[index_sfx] = &_uniq_tp;
					}
					else {
						uint64_t kups_ctx = 0;
						// label-2 (implied by else) and non-unique transitions
						double *ctx_tp = (double *) mempool_alloc(mp, 
								tnext_kmer->n_trans * sizeof(double));

						// kups_ctx is a bit-field consisting of 4-bit
						// unsigned integers.
						kups_ctx = hmmsbg_load_ctx_trans_p(d, read, tnext_kmer, t+1, ctx_tp);
						renormalize_trans_p(tnext_kmer->n_trans, ctx_tp);
						next_nbhd->kmer_ups_ctx[index_sfx] = kups_ctx;
						next_nbhd->kmer_trans_prob[index_sfx] = ctx_tp;
					}
				}

				/* set up transition flag */
				bitarray_set_pwr2(next_nbhd->ba_kmer_trans_flag,
						index_sfx, tnext_kmer->trans_flag, 2);
				/* set up Hamming distance */
				bitarray_set_pwr2(next_nbhd->ba_hamming_dist, 
						index_sfx, 
						sfx_hd + ((mismatch_flag >> dns_base) & 1), 3);

				next_nbhd->suffixes_order[index_pfx] = index_sfx;
				//next_nbhd->states_sorted_pfx[index_pfx] = tnext_kid;
				next_nbhd->states_sorted_sfx[index_sfx] = tnext_kid;
				next_nbhd->alphas[index_sfx] = alpha;
				next_nbhd->emit_dens[index_sfx] = emit_dens;
#if !PMR_USE_PRIMITIVE_EMISSION
				next_nbhd->base_emit_ctx[index_sfx] = bctx;
				next_nbhd->qual_emit_ctx[index_sfx] = qctx;
#endif
				index_pfx++;
			}

			kmer_idx += n_common_sfx_kmers;
		}

		obs_base = obs_next_base;
		obs_n_flag = next_obs_n_flag;
	}

	tmax = t;

	/* compute log-likelihood at the last position */
	double hmm_loglik = log_sum(T_nbhds[tmax].size, T_nbhds[tmax].alphas,
			tnext_argmax_alpha, tnext_max_alpha);

	if (!isnan(hmm_loglik) && tmax > tmin) {
		int wmax = (tmax-1) / d->opt->emit_window_size;

		int tmax_kp = tmax + kmer_len - 1;
		state_nbhd_t *tmax_nbhd = &T_nbhds[tmax];
		int tmax_nbhd_size = tmax_nbhd->size;

		atomic_add_dbl(&d->loglik, hmm_loglik);
#if !PMR_USE_PRIMITIVE_EMISSION
	//#pragma omp critical
	//{
		double *tmax_base_emit_exp = &d->base_emit_exp[wmax * 
			d->base_ctx_radix + sz_base_emit_p_off];
		double *tmax_qual_emit_exp = &d->qual_emit_exp[wmax *
			d->qual_ctx_radix + sz_qual_emit_p_off];

		/* --- update emission expected values --- */
		kbits_t xt = kmer_effective_base(rseq[tmax_kp]);
		int yt = (int) conv_q[tmax_kp];

		for (int i = 0; i < tmax_nbhd_size; i++) {
			double norm_exp_count = EXP(PRODUCT(tmax_nbhd->alphas[i],
					-hmm_loglik));

			int bctx = tmax_nbhd->base_emit_ctx[i]; 
			int qctx = tmax_nbhd->qual_emit_ctx[i]; 

			atomic_add_dbl(&tmax_base_emit_exp[bctx * 5 + xt],
					norm_exp_count);
			atomic_add_dbl(&tmax_qual_emit_exp[qctx * qmax + yt],
					norm_exp_count);
		}
	//}
#endif
	}
	else {
		mempool_destroy(mp);
		return;
	}

	state_nbhd_t *next_nbhd = &T_nbhds[tmax];
	/* set BETA's at the tmax to be 1.0 */
	/*
	memset(next_nbhd->betas, 0, 
			next_nbhd->size * sizeof(*next_nbhd->betas));
	*/

	/* ---- backward algorithm ---- */
	for (int t = tmax - 1; t >= tmin; t--) {
		state_nbhd_t *curr_nbhd = &T_nbhds[t],
					 *next_nbhd = &T_nbhds[t+1];

		const int t_kp = t + kmer_len - 1;
		int t_nbhd_size = curr_nbhd->size;

		/* expected counts and expected transitions, the gamma's and
		 * xi's in Rabiner's notation. */
		double *t_exp_counts = (double *) mempool_alloc(mp, sizeof(*t_exp_counts) *
				t_nbhd_size);
		double *t_exp_trans = (double *) mempool_alloc(mp, sizeof(*t_exp_trans) *
				curr_nbhd->n_trans);

		//memset(t_exp_counts, 0, sizeof(*t_exp_counts) * t_nbhd_size);
		//memset(t_exp_trans, 0, sizeof(*t_exp_trans) * curr_nbhd->n_trans);

		int n_distinct_suffixes = curr_nbhd->n_uniq_sfx;
		int t_kmer_idx = 0;
		int tnext_kmer_idx = 0;
		double *kmer_exp_trans = t_exp_trans;
#if PMR_USE_LOG_ADD
		double sum_exp_count = 0.0;	/* for normalization */
#else
		int argmax_gamma = -1;
		double max_gamma = NAN;
#endif
		for (int i = 0; i < n_distinct_suffixes; i++) {
			kbits_t _tf = bitarray_get_pwr2(curr_nbhd->ba_pfx_sfx_flags,
					i, 3);
			kbits_t ups_trans_packed = trans_flag2packed(_tf & 15);
			kbits_t dns_trans_packed = trans_flag2packed(_tf >> 4);

			register int j;
			for (j = ups_trans_packed & 7, ups_trans_packed >>= 3;
					j > 0; --j, ups_trans_packed >>= 2) {

				kbits_t k_tf = bitarray_get_pwr2(
						curr_nbhd->ba_kmer_trans_flag, t_kmer_idx, 2);
				int k_n_trans = trans_flag2packed(k_tf) & 7;
				//int k_n_trans = curr_nbhd->kmer_ptrs[t_kmer_idx]->n_trans;
				for (int k = 0; k < k_n_trans; k++) {
					kmer_exp_trans[k] = NAN;
				}

				kmer_t *pk = curr_nbhd->kmer_ptrs[t_kmer_idx];
				double _alpha = curr_nbhd->alphas[t_kmer_idx];
				double _beta = hmm_compute_beta(
						dns_trans_packed,
						tnext_kmer_idx,  /* prefix order */
						k_tf, //pk,
						next_nbhd->betas, 
						curr_nbhd->kmer_trans_prob[t_kmer_idx], 
						kmer_exp_trans, next_nbhd->suffixes_order,
						next_nbhd->emit_dens, _alpha);

				curr_nbhd->betas[t_kmer_idx] = _beta;

				double _gamma = PRODUCT(_alpha, _beta);
				t_exp_counts[t_kmer_idx] = _gamma;

#if PMR_USE_LOG_ADD
				sum_exp_count += EXP(_gamma);
#else
				update_max_quantity(_gamma, t_kmer_idx,
						max_gamma, argmax_gamma);
#endif
				kmer_exp_trans += k_n_trans;
				t_kmer_idx++;
			}
			
			tnext_kmer_idx += (dns_trans_packed & 7);
		}


#if PMR_USE_LOG_ADD
		double inv_sum_exp_count = -LOG(sum_exp_count);
#else
		double inv_sum_exp_count = -log_sum(curr_nbhd->size, t_exp_counts, 
				argmax_gamma, max_gamma);
#endif

#if PMR_USE_PRIMITIVE_EMISSION
		double *t_base_emit_exp = NULL;
		double *t_qual_emit_exp = NULL;
#else
		int widx = (t-1) / d->opt->emit_window_size;
		double *t_base_emit_exp = (t == 0) ? 
			d->base_emit_exp : 
			&d->base_emit_exp[widx * d->base_ctx_radix + sz_base_emit_p_off];

		double *t_qual_emit_exp = (t == 0) ? 
			d->qual_emit_exp : 
			&d->qual_emit_exp[widx * d->qual_ctx_radix + sz_qual_emit_p_off];

		// thread-specific offset
		t_base_emit_exp += omp_get_thread_num() * 
			(sz_base_emit_p_off + d->n_emit_windows * d->base_ctx_radix);
		t_qual_emit_exp += omp_get_thread_num() *
			(sz_qual_emit_p_off + d->n_emit_windows * d->qual_ctx_radix);
#endif

		/* --- update expected values (super redudant comment) --- */
		hmmsbg_update_expected_vals(d, curr_nbhd, t,
				t_exp_counts, t_exp_trans,
				t_base_emit_exp, t_qual_emit_exp,
				kmer_len, ctx_shift, qmax, rseq, 
				conv_q, inv_sum_exp_count);
	}

destroy:
	mempool_destroy(mp);
}

void hmmsbg_viterbi(data *d, read_t *read, void *fdata)
{
	int actg_trans_count[4];
	int cum_actg_trans_count[4];

	const int kmer_len = d->opt->kmer_length;
	const int num_kmers = read->length - kmer_len + 1;
	const kbits_t klen_flag = (1UL << (kmer_len << 1)) - 1;
	const int shift = (kmer_len - 1) << 1;
	const int ctx_shift = (kmer_len - d->opt->base_emit_context) << 1;
	const int qmax = d->opt->max_qual_score + 1;
	const int qmax2 = qmax << 1;
	const int bwidth = d->opt->qscore_bin_width;
	double _uniq_tp = 0.0;

	intptr_t revcomp = (intptr_t) fdata;

	// in reverse complementary decoding mode, must use the first kmer (which
	// is the last kmer on the original read).
	int tmin = d->read_start_positions[read->id];
	if (revcomp && tmin > 0) {
		tmin = read->length - 1 - (d->read_start_positions[read->id] + kmer_len - 1);
	}
	int tmax = read->length - kmer_len;
	if (tmin >= tmax) {
		tmin = 0;	// try starting from the 5' end (since we have larger dmax)
		//return;
	}

	mempool_t *mp = mempool_create(MEMPOOL_DEFAULT_SIZE);

	dict *kdict = d->kmer_dict;

	/* sequence */
	char *rseq = read->sequence;
	char *conv_q = convert_quality_scores(read->qscore, read->length, 
			d->opt->qual_score_offset, bwidth, mp);

	state_nbhd_t *T_nbhds = (state_nbhd_t *) mempool_nalloc(mp,
			sizeof(*T_nbhds) * num_kmers, 16);

	read_t argmax_read = { read->length, 0, 0, NULL, NULL, NULL };
	argmax_read.sequence = (char *) mempool_alloc(mp, 
			read->length + 1);

	/* flag for N (not determined) base, first kmer only */
	kbits_t obs_n_flag = kmer_n_base_flag(rseq + tmin, kmer_len);
	kbits_t next_obs_n_flag = kmer_n_base_flag(rseq + tmin + 1, kmer_len);

	double tnext_max_delta = NAN;
	int tnext_argmax_delta = -1;

	/* convert string literal to numeric representation */
	kbits_t observed_kmers[BITS_TO_U64(read->length << 1)];
	read_seq_to_numeric(rseq, observed_kmers, read->length, kmer_len);

	/* set up first kmer neighborhood */
	kbits_t obs_kid = bitarray_get(observed_kmers, tmin << 1, kmer_len << 1);
	ksubstr_t obs_base = kmer_effective_base(rseq[tmin+kmer_len-1]);
	ksubstr_t obs_next_base = kmer_effective_base(rseq[tmin+kmer_len]);

	T_nbhds[tmin].size = hmmsbg_load_preconstructed_nbhd(
			d, read, obs_kid, 
			obs_n_flag, conv_q, &T_nbhds[tmin], kmer_len, 1, kdict, mp, 
			const_cast<double *>(&_uniq_tp));

	int wsize = d->opt->hd_window_size;
	uint8_t *win_dmax = &d->read_dmax[read->id * d->n_hd_windows];

	/* initial state distribution */
	size_t sz_base_emit_p_off = kmer_len * 20;
	size_t sz_qual_emit_p_off = kmer_len * qmax2;
	uint8_t dmax = 0;
	int t = tmin;
	int early_decode = 0;
	for (; t < tmax; t++) {
		int tnext_kp = t + kmer_len;
		int w_emit = t / d->opt->emit_window_size;

#if !PMR_USE_PRIMITIVE_EMISSION
		double *tnext_base_emit_p = &d->base_emission_p[w_emit * 
			d->base_ctx_radix + sz_base_emit_p_off];
		double *tnext_qual_emit_p = &d->qual_emission_p[w_emit *
			d->qual_ctx_radix + sz_qual_emit_p_off];
#endif

		obs_next_base = kmer_effective_base(rseq[t + kmer_len]);
		obs_kid = bitarray_get(observed_kmers, t << 1, kmer_len << 1);
		kbits_t obs_next_kid = bitarray_get(observed_kmers, (t+1) << 1, 
				kmer_len << 1);

		kmer_t *obs_pk = (kmer_t *) &(dict_find(kdict,
					kbits_cast_to_ptr(obs_kid))->value);

		register kbits_t _new_nf = ((obs_next_base >> 2) << 1) | 
			(obs_next_base >> 2);
		next_obs_n_flag = (next_obs_n_flag >> 2) | (_new_nf << shift);

		state_nbhd_t *curr_nbhd = &T_nbhds[t],
					 *next_nbhd = &T_nbhds[t+1];

		dmax = win_dmax[t / wsize]; 

		int t_nbhd_size = curr_nbhd->size;

		kbits_t *t_states_sorted_sfx = T_nbhds[t].states_sorted_sfx;
		double *t_states_deltas = T_nbhds[t].alphas;
		double **t_kmer_trans_p = T_nbhds[t].kmer_trans_prob;

		/* ----- II. compute one iteration of Viterbi algorithm ----- */

		/* --- II a. set up kmer neighborhood at next position --- */
		memset(actg_trans_count, 0, sizeof(int) << 2);

		if (curr_nbhd->size == 1 && 
				curr_nbhd->kmer_ptrs[0]->last_kmer &&
				curr_nbhd->kmer_ptrs[0]->trans_flag == 0) {
			early_decode = 1;
			break;
		}

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

		/* allocate memory given tnext_nbhd_size */
		char *pmem = (char *) mempool_nalloc(mp,
				hmmsbg_nbhd_alloc_size(tnext_nbhd_size), 16);

		next_nbhd->size = tnext_nbhd_size;
		hmmsbg_setup_nbhd_ptrs(next_nbhd, tnext_nbhd_size, pmem);

		if (tnext_nbhd_size == 0) {
			/* if somehow all pathways terminated up to this point, perform
			 * backtrack immediately */
			break;
		}

		int kmer_idx = 0, index_pfx = 0;

		tnext_max_delta = NAN;
		tnext_argmax_delta = -1;

		register kbits_t mismatch_flag = ~(1UL << obs_next_base) |
			~((obs_next_base >> 2) - 1UL);

		for (int i = 0; i < n_distinct_suffixes; i++) {
			int n_common_sfx_kmers = bitarray_get_pwr2(
					curr_nbhd->ba_distinct_sfx, i, 1) + 1;

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

#if PMR_USE_PRIMITIVE_EMISSION
				double emit_dens = dns_base == obs_next_base ? 
					LOG(1.0-d->base_error_rate) : LOG(d->base_error_rate/3);
#else
				/* if already dmax errors in the (k-1)-prefix, can only emit
				 * the observed base. */
				kbits_t emit_renorm_flag = (sfx_hd == dmax) ? 
					((1U << obs_next_base) & 0xF) : 15; //obs_next_kmer->emit_norm_flag;

				/* compute the emission probability */
				int bctx = calc_bctx(tnext_kid, sfx_hd, shift, tnext_kmer->single_stranded);
				int qctx = (tnext_kid >> shift) == (obs_next_kid >> shift);
				//int qctx = calc_qctx(tnext_kid, obs_next_kid, conv_q + t + 1, kmer_len);

				double emit_dens = hmm_compute_emit_dens(
						dns_base, obs_next_base, 
						(int) conv_q[tnext_kp], 
						emit_renorm_flag,
						tnext_base_emit_p + bctx * 5,
						tnext_qual_emit_p + qctx * qmax);

				/*
				if (dns_base != obs_next_base) {
					if (sfx_hd > 0) emit_dens += 1.5;
					if (sfx_hd > 3) emit_dens += 0.7;
				}
				*/

				/*
				emit_dens = PRODUCT(emit_dens, 
						d->avg_initial_p[kmer_initial_index(obs_pk,
							d->uniq_kmer_init_off)]);
				*/
#endif

				/* we will be use betas as the psi's following Rabiner's
				 * notation. */
				double delta;
				int64_t countdown = 0;

				int argmax_delta = hmmsbg_compute_delta(d, 
						ups_trans_packed, kmer_idx,
						curr_nbhd, &countdown,
						emit_dens, dns_base, &delta); 

				/* if we are at the end of the read, 
				 * check if the kmer is specially flagged as a "last kmer",
				 * i.e. a kmer that cannot be further extended. 
				 * all kmer designated as the "last kmer" will be assigned zero
				 * likelihood */
				// if (t == (tmax - 1) && tnext_kmer->last_kmer) delta = NAN;
				
#if PMR_EXTEND_READS
				/*
				if (t == tmax - 1) {
					// for the last kmer, incorporate the probability of read
					// extension.
					double pext = d->prob_kmer_extension[
							hmmsbg_initial_index(d, tnext_kmer)];
					delta = PRODUCT(delta, pext);
				}
				*/
#endif

				update_max_quantity(delta, index_sfx,
						tnext_max_delta, tnext_argmax_delta);

				next_nbhd->kmer_ptrs[index_sfx] = tnext_kmer;
				next_nbhd->kmer_trans_prob[index_sfx] = 
					hmmsbg_locate_trans_p(d, tnext_kmer, &_uniq_tp);

#if PMR_JUMP
				next_nbhd->jump_ahead_countdown[index_sfx] = countdown;

				if (countdown > 0) {
					kbits_t jump_step = kmer_effective_base(rseq[t + kmer_len + 1]);
					// set the transition probability accordingly
					double *ctx_tp = (double *) mempool_alloc(mp, 
							tnext_kmer->n_trans * sizeof(double));

					kbits_t tpacked = trans_flag2packed(tnext_kmer->trans_flag) >> 3;
					for (int b = 0; b < tnext_kmer->n_trans; ++b, tpacked >>= 2) {
						kbits_t base = tpacked & 3;
						if (base == jump_step) ctx_tp[b] = 0.0;
						else ctx_tp[b] = NAN;
					}

					next_nbhd->kmer_trans_prob[index_sfx] = ctx_tp;
				}
				else if (tnext_kmer->has_multi_copies) {
#else
				if (tnext_kmer->has_multi_copies) {
#endif
					if (tnext_kmer->label == 0 || tnext_kmer->n_trans == 1) {
						// perform backtrack to find the argmax sequence, which
						// is utilized to determine the upstream context.
						kmer_ext_t *kext = &d->nonuniq_kmers[tnext_kmer->index];
						int cdist = kext->ctx_distance;

						uint64_t kups_ctx = kext->n_ctx; 
						if (t + 1 >= cdist) {
							memcpy(argmax_read.sequence, read->sequence, read->length);

							int optimal_st_idx = argmax_delta;
							for (int tb = t; tb > t - cdist; --tb) {
								state_nbhd_t *nbhd = &T_nbhds[tb];

								hmm_correct_one_error(d, &argmax_read, tb,
										kbits_first_base(nbhd->states_sorted_sfx[optimal_st_idx]),
										kmer_effective_base(rseq[tb]));

								optimal_st_idx = (int) nbhd->betas[optimal_st_idx];
							}
							kups_ctx = hmmsbg_determine_ctx(d, &argmax_read, tnext_kmer, t+1);

						}
						//kups_ctx = hmmsbg_determine_ctx(d, read, tnext_kmer, t+1);

						// if kups_ctx == n_ctx, it suggests ambiguous context
						next_nbhd->kmer_ups_ctx[index_sfx] = kups_ctx;
						next_nbhd->kmer_trans_prob[index_sfx] += 
							tnext_kmer->n_trans * kups_ctx;

						if (tnext_kmer->n_trans == 1) // && 
								//!IS_ZERO(next_nbhd->kmer_trans_prob[index_sfx][0]))
							next_nbhd->kmer_trans_prob[index_sfx] = &_uniq_tp;
					}
					else {
						uint64_t kups_ctx = 0;
						// label-2 (implied by else) and non-unique transitions
						double *ctx_tp = (double *) mempool_alloc(mp, 
								tnext_kmer->n_trans * sizeof(double));

						// kups_ctx is a bit-field consisting of 4-bit
						// unsigned integers.
						kups_ctx = hmmsbg_load_ctx_trans_p(d, read, tnext_kmer, t+1, ctx_tp);
						renormalize_trans_p(tnext_kmer->n_trans, ctx_tp);
						next_nbhd->kmer_ups_ctx[index_sfx] = kups_ctx;
						next_nbhd->kmer_trans_prob[index_sfx] = ctx_tp;
					}
				}

				/* set up transition flag */
				bitarray_set_pwr2(next_nbhd->ba_kmer_trans_flag,
						index_sfx, tnext_kmer->trans_flag, 2);
				/* set up Hamming distance */
				bitarray_set_pwr2(next_nbhd->ba_hamming_dist, 
						index_sfx, 
						sfx_hd + ((mismatch_flag >> dns_base) & 1), 3);

				next_nbhd->states_sorted_sfx[index_sfx] = tnext_kid;
				next_nbhd->alphas[index_sfx] = delta;
				next_nbhd->betas[index_sfx] = (double) argmax_delta;
			}

			kmer_idx += n_common_sfx_kmers;
		}

		obs_base = obs_next_base;
		obs_n_flag = next_obs_n_flag;
	}

	if (((early_decode == 0) && (d->opt->partial_decoding_lvl == 0)) && t < tmax) {
		mempool_destroy(mp);

		FILE *stream = d->errors_fp ? d->errors_fp : stdout;
		fprintf(stream, "@%s\n%.*s\n+%s\n%.*s\n",
				read->identifier, 
				read->length, read->sequence,
				read->identifier, read->length, read->qscore
			);
		return;
	}

	int end_of_read = (tmax == t);

	/* identify the most likely pathway at the last position */
	/* and perform partial decoding */
	kbits_t tfin_obs_kid = (t < tmax) ? obs_kid : 
		(obs_kid >> 2) | (obs_next_base << shift);

	tmax = t;

	if (!IS_ZERO(tnext_max_delta)) {
		state_nbhd_t *nbhd = &T_nbhds[tmax];

		if (end_of_read || early_decode) { //  d->opt->partial_decoding_lvl == PDEC_AGGRESSIVE) {
			/* for conservative partial decoding, do not correct errors within
			 * the last kmer to avoid excessive false positives. */
			hmm_correct_errors(d, read, tmax,
					nbhd->states_sorted_sfx[tnext_argmax_delta],
					tfin_obs_kid, obs_n_flag,
					kmer_len);
		}

		int optimal_st_idx = (int) nbhd->betas[tnext_argmax_delta];
		//if (!end_of_read) 
			//fprintf(stream, "# ---- [partial decoding: %d] ----\n", read->id);

		for (int t = tmax - 1; t >= tmin; t--) {
			nbhd = &T_nbhds[t];

			hmm_correct_one_error(d, read, t,
					kbits_first_base(nbhd->states_sorted_sfx[optimal_st_idx]),
					kmer_effective_base(rseq[t]));

			optimal_st_idx = (int) nbhd->betas[optimal_st_idx];
		}
	}

	FILE *stream = d->errors_fp ? d->errors_fp : stdout;
	fprintf(stream, "@%s\n%.*s\n+%s\n%.*s\n",
			read->identifier, 
			read->length, read->sequence,
			read->identifier, read->length, read->qscore
		);

	mempool_destroy(mp);
}

// this function only applies to `pk` which has:
//	* pk->label == 1
//	* pk->has_multi_copies == 1
//	* pk->n_trans > 1
static uint64_t hmmsbg_load_ctx_trans_p(data *d, read_t *read, kmer_t *pk, int t, 
		double *tp)
{
	int k = d->opt->kmer_length;
	int shift = (k - 1) << 1;

	kmer_ext_t *kext = &d->nonuniq_kmers[pk->index];
	double *pk_tp = d->transition_p + kext->poff;
	kbits_t kid = kbits_id_of(pk);

	uint64_t ctx = 0;
	int tp_par_off = 0;

	kbits_t tpacked = trans_flag2packed(pk->trans_flag);
	int n_trans = tpacked & 7, i = 0;
	for (tpacked >>= 3; i < n_trans; ++i, tpacked >>=2) {
		kbits_t base = tpacked & 3;
		kbits_t kid_lbl1 = kmer_reverse_complement(
				(kid >> 2) | (base << shift), k);

		uint32_t n_ctx_i = (kext->n_ctx >> (i * PMR_MAX_CTX_BITS)) & PMR_MAX_NCTX;
		uint64_t ctx_i = 0;

		if (n_ctx_i > 0) {
			kmer_t *pk_lbl1 = (kmer_t *) &(dict_find(d->kmer_dict, 
						kbits_cast_to_ptr(kid_lbl1))->value);
			ctx_i = hmmsbg_determine_ctx(d, read, pk_lbl1, t, true);
		}

		tp[i] = pk_tp[tp_par_off + ctx_i];
		ctx |= (ctx_i << (i * PMR_MAX_CTX_BITS));

		tp_par_off += (n_ctx_i + 1); 
	}
	
	return ctx;
}

/* Convert the upstream/downstream sequence into context indices */ 
static uint64_t hmmsbg_determine_ctx(data *d, read_t *read, kmer_t *pk, int t, 
		bool rev_strand)
{
	int k = d->opt->kmer_length;

	kmer_ext_t kext = d->nonuniq_kmers[pk->index];
	kbits_t *ctxs = &d->kmer_ups_contexts[kext.ctxoff];

	register int nctx = kext.n_ctx;
	kbits_t mask = 0;
	int ctxlen = 0;
	if (kext.ctx_distance >= k) {
		ctxlen = k;
		mask = UINT64_MAX;
	}
	else {
		ctxlen = kext.ctx_distance;
		mask = (1UL << (ctxlen << 1)) - 1;
	}

	kbits_t needle = 0;
	if (pk->label == 0 && rev_strand == false) {
		if (t < kext.ctx_distance) return nctx;
		needle = kmer_seq_to_id(&read->sequence[t - kext.ctx_distance], ctxlen);
	}
	else {
		++kext.ctx_distance;
		// if (rev_strand) ++kext.ctx_distance;

		int toff = t + kext.ctx_distance;
		if (toff + k - 1  >= read->length) return nctx;
		needle = kmer_reverse_complement(
				kmer_seq_to_id(&read->sequence[toff], k), k) & mask;
	}

	auto hit = std::find(ctxs, ctxs + nctx, needle);

	return (int) (hit - ctxs);
}

void hmmsbg_coalesce_lbl2_trans(data *d)
{
	double _uniq_tp = 0.0;
	int kmer_len = d->opt->kmer_length;
	int shift = (kmer_len - 1) << 1;

	int64_t n_kmers = d->kmer_dict->sizemask + 1;

	for (int64_t i = 0; i < n_kmers; ++i) {
		dictentry *de = &d->kmer_dict->table[i];
		if (de->value == NULL) continue;

		kmer_t *pk = (kmer_t *) &(de->value);
		kbits_t kid = kbits_id_of(pk);

		if (pk->has_multi_copies == 0 || pk->label == 0) continue;

		kmer_ext_t *pkext = &d->nonuniq_kmers[pk->index];
		double *tp = hmmsbg_locate_trans_p(d, pk, &_uniq_tp);

		uint32_t n_ctx = pkext->n_ctx;
		int par_off = 0;
		kbits_t tpacked = trans_flag2packed(pk->trans_flag) >> 3;

		// use the sum for the ambiguous contexts
		for (int b = 0; b < pk->n_trans; ++b, tpacked >>= 2) {
			uint32_t n_ctx_b = (n_ctx >> (b * PMR_MAX_CTX_BITS)) & PMR_MAX_NCTX;
			uint32_t n_params = (n_ctx_b == 0) ? 1 : n_ctx_b + 1;

			kbits_t base = tpacked & 3;
			kbits_t k1 = kmer_reverse_complement(
					(kid >> 2) | (base << shift), kmer_len);
			kmer_t *pk1 = (kmer_t *) &(dict_find(
						d->kmer_dict, kbits_cast_to_ptr(k1))->value);

			if (pk1->has_multi_copies == 0) {
				par_off += n_params;
				continue;
			}

			for (int c = 0; c < n_ctx_b; ++c) {
				//if (!IS_ZERO(tp[par_off + c]))
					tp[par_off + c] = tp[par_off + n_ctx_b];  // use the sum if the context is not resolvable. 
			}

			par_off += n_params;
		}
	}
}
void hmmsbg_m_step(data *d)
{
	int kmer_len = d->opt->kmer_length;
	int shift = (kmer_len - 1) << 1;

	double jnt_exp_trans[4] = {0.0, 0.0, 0.0, 0.0};
	// double accmu_exp_trans[4] = {0.0, 0.0, 0.0, 0.0};
	double penalized_tp[4] = {0.0, 0.0, 0.0, 0.0};
	double _etas[4] = {0.0};
	double _uniq_tp = 0.0;

	double thresh = d->opt->penalty_eta / log1p(1/d->opt->penalty_gamma);
	double loglik_penalty = 0.0;
	double sum_k_cexp_cnt = 0.0;

	const double rho_one = 1.1 * log1p(1.0/d->opt->penalty_gamma);

	mstep_data_t md = {
		d->opt->penalty_gamma, 
		_etas, nullptr, 0, 0
	};

	int64_t n_kmers = d->kmer_dict->sizemask + 1;

	double *_init_p = (double *) calloc(d->n_nonzero_kmers + 1,
			sizeof(*_init_p));
	std::vector<double> _multictx_init_p;
	//_multictx_init_p.reserve();

	// ---- label 1 kmer ----
	INFO_TEE(d->log_fp, "Estimating transition probabilities for label-1 kmers...\n");
	for (int64_t i = 0; i < n_kmers; ++i) {
		dictentry *de = &d->kmer_dict->table[i];
		if (de->value == nullptr) continue;

		kmer_t *pk = (kmer_t *) &(de->value);

		if (pk->label == 1) continue;

		int n = pk->n_trans;
		kbits_t kid = kbits_id_of(pk);
		//int n_ctx = pk->has_multi_copies ? d->nonuniq_kmers[pk->index].n_ctx : 1;
		int n_ctx = pk->has_multi_copies ? d->nonuniq_kmers[pk->index].n_ctx : 0;
		double accmu_k_cexp_cnt = 0.0;

		// memset(accmu_exp_trans, 0, 4 * sizeof(*accmu_exp_trans));

		// tweak \eta on the fly
		double _eta = //pk->has_multi_copies ? d->opt->penalty_eta / n_ctx : 
			d->opt->penalty_eta;
		_eta = std::max(_eta, rho_one);
		std::fill(_etas, _etas + 4, _eta);

		thresh = _eta / log1p(1/d->opt->penalty_gamma);

		for (int c = 0; c <= n_ctx; ++c) {
		//for (int c = 0; c < n_ctx; ++c) {
			double *exp_trans_p = hmmsbg_locate_exp_trans(d, pk) + pk->n_trans * c; 
			double *cur_tp = hmmsbg_locate_trans_p(d, pk, &_uniq_tp) + pk->n_trans *c;

			hmmsbg_multictx_lbl2_exp_trans(d, pk, c, &jnt_exp_trans[0]);
			for (int j = 0; j < n; ++j) {
				jnt_exp_trans[j] += exp_trans_p[j];
				//accmu_exp_trans[j] += jnt_exp_trans[j];
			}

			memcpy(&penalized_tp[0], cur_tp, n * sizeof(double));
			
			tflag_t new_tf = {pk->trans_flag, n};
			double k_cexp_cnt = 0.0;
			loglik_penalty += cmstep_estimate_trans_p(d, pk, &penalized_tp[0], 
					&jnt_exp_trans[0],
					// (c == n_ctx) ? &accmu_exp_trans[0] : &jnt_exp_trans[0],
					&k_cexp_cnt, pk->trans_flag, thresh, &md, &new_tf); 

			if (pk->has_multi_copies || pk->n_trans > 1) {
				memcpy(cur_tp, penalized_tp, n * sizeof(double));
			}

			sum_k_cexp_cnt += k_cexp_cnt;
			accmu_k_cexp_cnt += k_cexp_cnt;
			memcpy(exp_trans_p, jnt_exp_trans, n * sizeof(double));

			if (c < n_ctx) {
				if (pk->has_multi_copies)
					_multictx_init_p.push_back(k_cexp_cnt);

				// update the old transition probabilities in situ
				// memcpy(exp_trans_p, jnt_exp_trans, n * sizeof(double));
			}
			else {
				if (pk->has_multi_copies) {
					//_multictx_init_p.push_back(accmu_k_cexp_cnt);
					_multictx_init_p.push_back(k_cexp_cnt);
					// total counts
					_multictx_init_p.push_back(accmu_k_cexp_cnt);
				}
				else
					_init_p[hmmsbg_initial_index(d, pk)] = k_cexp_cnt;

				// update the old transition probabilities in situ
				// memcpy(exp_trans_p, accmu_exp_trans, n * sizeof(double));
			}

			/*
			if (c == n_ctx - 1) {
				if (pk->has_multi_copies) {
					_multictx_init_p.push_back(accmu_k_cexp_cnt);
				}
				else
					_init_p[hmmsbg_initial_index(d, pk)] = k_cexp_cnt;
			}
			*/


			// new transition flags will be ignored.
			// if a kmer with multiple exit transitions now has a single
			// non-zero transition, we will still keep all the flags, although
			// some of the probabilities will be zero.

		}

		pk->dirty = 1;

	}

	INFO_TEE(d->log_fp, "Estimating transition probabilities for label-2 kmers...\n");

	// ---- label 2 kmer ----
	for (int i = 0; i < n_kmers; ++i) {
		dictentry *de = &d->kmer_dict->table[i];
		if (de->value == NULL) continue;

		kmer_t *pk = (kmer_t *) &(de->value);

		if (pk->label == 0) continue;
		int n = pk->n_trans;
		kbits_t kid = kbits_id_of(pk);

		double k_cexp_cnt = 0.0;
		k_cexp_cnt = hmmsbg_update_lbl2_trans_p(d, pk, _init_p);

		_init_p[hmmsbg_initial_index(d, pk)] = k_cexp_cnt;
		sum_k_cexp_cnt += k_cexp_cnt;

		pk->dirty = 1;
	}

	// compute the initial state distribution parameters.
	for (int i = 0; i < n_kmers; ++i) {
		dictentry *de = &d->kmer_dict->table[i];
		if (de->value == NULL) continue;

		kmer_t *pk = (kmer_t *) &(de->value);

		if (pk->dirty == 0) continue;

		kbits_t kid = kbits_id_of(pk);
		kbits_t rev_kid = kmer_reverse_complement(kid, kmer_len);
		kmer_t *rev_pk = (kmer_t *) &(dict_find(
					d->kmer_dict, kbits_cast_to_ptr(rev_kid))->value);

		int _idx = hmmsbg_initial_index(d, pk);
		int _ridx = hmmsbg_initial_index(d, rev_pk);
		if (pk->has_multi_copies == 0 &&
				rev_pk->has_multi_copies == 0) {

			/* EXPERIMENTAL: FIXME!!! */
			double avg_kcexp_cnt = (_init_p[_idx] + _init_p[_ridx]) / 2.0;
			/*
			if (avg_kcexp_cnt < thresh &&
					(pk->trans_flag == 0 || rev_pk->trans_flag == 0)) {
				// if a kmer is below the threshold, and appears to be the
				// first kmer, push its initial state prob. to zero, 
				// mimicking the penalty on the initial state distribution.
				avg_kcexp_cnt = d->opt->penalty_gamma;
			}
			*/
			/* END OF EXPERIMENTAL CODE: FIXME!! */
			double mean_init_p = LOG(avg_kcexp_cnt / sum_k_cexp_cnt);

			d->avg_initial_p[_idx] = mean_init_p;
			d->avg_initial_p[_ridx] = mean_init_p;
		}
		else {
			if (pk->has_multi_copies && pk->label == 0) {
				uint32_t n_ctx = d->nonuniq_kmers[pk->index].n_ctx;
				int _cidx = (int) d->avg_initial_p[_idx];	
				for (int c = 0; c <= n_ctx + 1; ++c) {
					d->multictx_initial_p[_cidx + c] = LOG(_multictx_init_p[_cidx + c] /
						sum_k_cexp_cnt);
				}
			}
			else {
				d->avg_initial_p[_idx] = LOG(_init_p[_idx] / sum_k_cexp_cnt);
			}

			if (rev_pk->has_multi_copies && rev_pk->label == 0) {
				uint32_t n_ctx = d->nonuniq_kmers[rev_pk->index].n_ctx;
				int _cidx = (int) d->avg_initial_p[_ridx];	
				for (int c = 0; c <= n_ctx + 1; ++c) {
					d->multictx_initial_p[_cidx + c] = LOG(_multictx_init_p[_cidx + c] /
						sum_k_cexp_cnt);
				}
			}
			else {
				d->avg_initial_p[_ridx] = LOG(_init_p[_ridx] / sum_k_cexp_cnt);
			}
		}

		pk->dirty = 0;
		rev_pk->dirty = 0;
	}

	// compute the weighted average probabilities for ambiguous contexts
	for (int i = 0; i < n_kmers; ++i) {
		dictentry *de = &d->kmer_dict->table[i];
		if (de->value == NULL) continue;

		kmer_t *pk = (kmer_t *) &(de->value);
		kbits_t kid = kbits_id_of(pk);

		if (pk->has_multi_copies == 0) continue;

		kmer_ext_t *pkext = &d->nonuniq_kmers[pk->index];
		double *tp = hmmsbg_locate_trans_p(d, pk, &_uniq_tp);
		double *initp = hmmsbg_initial_p(d, pk);
		if (pk->label == 0) {
			for (int b = 0; b < pk->n_trans; ++b) {
				double avg_tp = 0.0;
				double denom = 0.0;
				for (int c = 0; c <= pkext->n_ctx; ++c) {
					denom += EXP(initp[c]);
					avg_tp += EXP(PRODUCT(initp[c], 
								tp[c * pk->n_trans + b]));
				}

				// overwrite the last set of trans. prob.
				tp[pkext->n_ctx * pk->n_trans + b] = 
					LOG(avg_tp / denom);
			}
		}
		else {
			uint32_t n_ctx = pkext->n_ctx;
			int par_off = 0;
			kbits_t tpacked = trans_flag2packed(pk->trans_flag) >> 3;

			// use the sum for the ambiguous contexts
			for (int b = 0; b < pk->n_trans; ++b, tpacked >>= 2) {
				uint32_t n_ctx_b = (n_ctx >> (b * PMR_MAX_CTX_BITS)) & PMR_MAX_NCTX;
				// uint32_t n_params = (n_ctx_b == 0) ? 1 : n_ctx_b + 1;

				kbits_t base = tpacked & 3;
				kbits_t k1 = kmer_reverse_complement(
						(kid >> 2) | (base << shift), kmer_len);
				kmer_t *pk1 = (kmer_t *) &(dict_find(
							d->kmer_dict, kbits_cast_to_ptr(k1))->value);

				if (pk1->has_multi_copies == 0) {
					par_off += (n_ctx_b + 1);
					//par_off += n_params;
					continue;
				}

				double *initp_pk1 = hmmsbg_initial_p(d, pk1);
				// double avg_tp = 0.0, denom = 0.0;
				double sum_tp = 0.0;
				for (int c = 0; c <= n_ctx_b; ++c) {
					// avg_tp += EXP(PRODUCT(tp[par_off + c], initp_pk1[c]));
					// denom += EXP(initp_pk1[c]);
					sum_tp += EXP(tp[par_off + c]);
				}

				sum_tp = LOG(sum_tp);
				// EXPERIMENTAL: always coalesce!!!
				for (int c = 0; c <= n_ctx_b; ++c) {
					//if (!IS_ZERO(tp[par_off + c]))
						tp[par_off + c] = sum_tp;  // use the sum if the context is not resolvable. 
				}
				// !!! END OF EXPERIMENTAL FEATURE !!!
				
				//tp[par_off + n_ctx_b] = LOG(avg_tp / denom);  // weighted average
				//tp[par_off + n_ctx_b] = LOG(sum_tp);  // sum of all contexts

				par_off += (n_ctx_b + 1);
				//par_off += n_params;
			}
		}
	}

	d->total_kmer_exp_counts = sum_k_cexp_cnt;

	mstep_emission(d);

	free(_init_p);
}

void hmmsbg_update_expected_vals(data *d, state_nbhd_t *curr_nbhd, 
		int t, double *t_exp_counts, double *t_exp_trans, 
		double *t_base_emit_exp, double *t_qual_emit_exp,
		int kmer_len, int ctx_shift, int qmax, 
		const char *rseq, const char *conv_q,
		double inv_sum_exp_count)
{
	double *ref_exp_trans[4];

	int t_kp = t + kmer_len - 1;
	const int shift = (kmer_len - 1) << 1;

	int t_kmer_idx = 0;
	int n_distinct_suffixes = curr_nbhd->n_uniq_sfx;
	double *kmer_exp_trans = t_exp_trans;

	for (int i = 0; i < n_distinct_suffixes; i++) {
		kbits_t _tf = bitarray_get_pwr2(curr_nbhd->ba_pfx_sfx_flags,
				i, 3);
		kbits_t ups_trans_packed = trans_flag2packed(_tf & 15);
		kbits_t dns_flag = _tf >> 4;

		register int j;
		register int n = ups_trans_packed & 7;
		for (j = n; j > 0; --j) {
			kmer_t *ref_k = curr_nbhd->kmer_ptrs[t_kmer_idx];
			kbits_t ref_kid = kbits_id_of(ref_k); 

			/* quick guess */
			double norm_exp_count = 1.0;
			if (n_distinct_suffixes > 1 || n > 1) {
				register double _norm_cnt = PRODUCT(t_exp_counts[t_kmer_idx], 
							inv_sum_exp_count);
				//if (fabs(_norm_cnt) < 1e-4) {
					// linear approximation to avoid __slowexp()
					//norm_exp_count = 1 + _norm_cnt;
				//	_slowexp++;
				//}
				//else {
					norm_exp_count = EXP(_norm_cnt);
				//}
			}

			
			kbits_t k_trans_packed = trans_flag2packed(
					ref_k->trans_flag);

			double *_kref_exp = hmmsbg_locate_exp_trans(d, ref_k);
			if (ref_k->has_multi_copies && ref_k->label == 1 &&
					ref_k->n_trans > 1) {
				// special case, where the memory layout is different
				uint64_t ctx = curr_nbhd->kmer_ups_ctx[t_kmer_idx];
				int par_off = 0;
				uint32_t n_ctx = d->nonuniq_kmers[ref_k->index].n_ctx;
				for (int b = 0; b < ref_k->n_trans; ++b, ctx >>= PMR_MAX_CTX_BITS, n_ctx >>= PMR_MAX_CTX_BITS) {
					uint64_t ctx_b = ctx & PMR_MAX_NCTX;
					ref_exp_trans[b] = &_kref_exp[par_off + ctx_b];
					par_off += ((n_ctx & PMR_MAX_NCTX) + 1);
				}
			}
			else {
				if (ref_k->has_multi_copies) {
					_kref_exp += ref_k->n_trans * curr_nbhd->kmer_ups_ctx[t_kmer_idx];
				}

				for (int b = 0; b < ref_k->n_trans; ++b, ++_kref_exp) {
					ref_exp_trans[b] = _kref_exp;
				}
			}

			int l = k_trans_packed & 7;
			for (int b = 0; b < l ; b++) {
				register double _norm_cnt = PRODUCT(kmer_exp_trans[b], 
							inv_sum_exp_count);
				/*
				if (fabs(_norm_cnt) < 1e-4) {
					//atomic_add_dbl(&ref_exp_trans_p[b], 1 + _norm_cnt);
					//ref_exp_trans_p[b] += (1 + _norm_cnt);

					_slowexp++;
				}
				*/
				//else {
					atomic_add_dbl(ref_exp_trans[b], EXP(_norm_cnt));
					//ref_exp_trans_p[b] += EXP(_norm_cnt);
				//}
			}

			/* -- update expected emissions -- */
			if (t > 0) {
				kbits_t st = ref_kid >> shift;
				kbits_t xt = kmer_effective_base(rseq[t_kp]);
				int yt = (int) conv_q[t_kp];

#if PMR_USE_PRIMITIVE_EMISSION
				atomic_add_dbl(&d->error_rate_emit_exp[(xt != st)],
						norm_exp_count);
#else
				int bctx = curr_nbhd->base_emit_ctx[t_kmer_idx];
				int qctx = curr_nbhd->qual_emit_ctx[t_kmer_idx];

				t_base_emit_exp[bctx * 5 + xt] += norm_exp_count;
				t_qual_emit_exp[qctx * qmax + yt] += norm_exp_count;

				/*
				atomic_add_dbl(&t_base_emit_exp[bctx * 5 + xt], norm_exp_count);
				atomic_add_dbl(&t_qual_emit_exp[qctx * qmax + yt], norm_exp_count);
				bctx_sum_[bctx * 5 + xt] += norm_exp_count;
				qctx_sum_[qctx * qmax + yt] += norm_exp_count;
				*/
#endif
			}
			else {
				for (int l = 0; l < kmer_len; l++, ref_kid >>= 2) {
					kbits_t si = ref_kid & 3;
					kbits_t xi = kmer_effective_base(rseq[l]);
					int yi = (int) conv_q[l];

#if PMR_USE_PRIMITIVE_EMISSION
					atomic_add_dbl(&d->error_rate_emit_exp[(xi != si)],
							norm_exp_count);
#else
					/*
					atomic_add_dbl(&t_base_emit_exp[l * 20 + si * 5 + xi], 
							norm_exp_count);
					
					atomic_add_dbl(&t_qual_emit_exp[((l << 1) + (si == xi)) * qmax + yi],
							norm_exp_count);
					*/

					t_base_emit_exp[l * 20 + si * 5 + xi] += norm_exp_count;
					t_qual_emit_exp[((l << 1) + (si == xi)) * qmax + yi] +=
						norm_exp_count;
#endif
				}
			}

			t_kmer_idx++;
			kmer_exp_trans += ref_k->n_trans;
		}
	}

	/*
	for (int i = 0; i < 20; ++i) {
		if (bctx_sum_[i] > 0.0)
			atomic_add_dbl(&t_base_emit_exp[i], bctx_sum_[i]); 
	}

	for (int i = 0; i < (qmax << 1); ++i) {
		if (qctx_sum_[i] > 0.0)
			atomic_add_dbl(&t_qual_emit_exp[i], qctx_sum_[i]);
	}

	delete[] bctx_sum_;
	delete[] qctx_sum_;
	*/
}

static void hmmsbg_debug_console(data *d)
{
	char cmd_buf[256];
	printf("PREMIER M-step debug console.\nUse '?' for supported commands.\n");
	while (1) {
		printf(">>> ");
		fgets(cmd_buf, 255, stdin);

		if (cmd_buf[0] == 'q') {
			printf("Exited from M step debug console.");
			fflush(stdout);
			break;
		}
		/* print out kmer information */

		if (cmd_buf[0] == 'k' || strncmp(cmd_buf, "kmer", 4) == 0) {
			kbits_t kid = 0;
			if (isalpha(cmd_buf[2])) {
				kid = kmer_seq_to_id(&cmd_buf[2], d->opt->kmer_length);
			}
			else {
				kid = strtol(&cmd_buf[2], NULL, 0);
			}

			hmmsbg_debug_kmer(d, kid);
		}

		if (cmd_buf[0] == 'r' || strncmp(cmd_buf, "rkmer", 5) == 0) {
			kbits_t kid = 0;
			if (isalpha(cmd_buf[2])) {
				kid = kmer_seq_to_id(&cmd_buf[2], d->opt->kmer_length);
			}
			else {
				kid = strtol(&cmd_buf[2], NULL, 0);
			}
			hmmsbg_debug_kmer(d, 
					kmer_reverse_complement(kid, d->opt->kmer_length));
		}

		// e 2222
		if (cmd_buf[0] == 'e' || strncmp(cmd_buf, "emit", 4) == 0) {

		}

		if (cmd_buf[0] == 'h' || cmd_buf[0] == '?' || strncmp(cmd_buf, "help", 4) == 0) {
			printf("Command list:\n\n"
					" h/?                Print out this help information.\n"
					" k <kmer|kmerid>    Query kmer information.\n"
					" r <kmer|kmerid>    Query a reverse complementary kmer.\n"
					" B <kmer|kmerid>    Perform a backtrack from the kmer.\n"
					" q                  Quit the debug console.\n");
		}

	}
}

static void hmmsbg_debug_kmer(data *d, kbits_t kid)
{
	const char *BASES = "ACTG";
	char kseq[32];
	double jnt_exp_trans[4] = {0.0, 0.0, 0.0, 0.0};
	double uniq_tp = 0.0;

	const int kmer_len = d->opt->kmer_length;
	const int shift = (kmer_len - 1) << 1; 


	dictentry *de = dict_find(d->kmer_dict,
			kbits_cast_to_ptr(kid));

	if (de == NULL) {
		printf("[ERROR] Kmer %lu does not exist.\n", kid);
		return;
	}

	kmer_t *pk = (kmer_t *) &(de->value);
	printf("[Kmer %lu]\n"
			" * Multi-context: %u\n"
			" * Label: %u\n"
			" * Index: %u\n",
			kid, pk->has_multi_copies,
			pk->label, 
			pk->index);

	if (pk->label == 0 && pk->has_multi_copies) {
		kmer_ext_t *pkext = &d->nonuniq_kmers[pk->index];
		printf(" * Expected counts: ");
		double *initp = hmmsbg_initial_p(d, pk);
		for (int c = 0; c <= pkext->n_ctx; ++c) {
			printf("%.2g ", d->total_kmer_exp_counts * EXP(initp[c]));
		}
		printf("\n");
	}
	else {
		printf(" * Expected counts: %g\n",
				d->total_kmer_exp_counts * EXP(*hmmsbg_initial_p(d, pk)));
	}

	if (pk->has_multi_copies) {
		kmer_ext_t *pkext = &d->nonuniq_kmers[pk->index];
		int n_ctx = pkext->n_ctx;
		printf(" * Contexts:\n");
		if (pk->label == 0 || pk->n_trans == 1) {
			kbits_t *ctxs = &d->kmer_ups_contexts[pkext->ctxoff];
			for (int c = 0; c < pkext->n_ctx; ++c) {
				int ctxlen = pkext->ctx_distance > kmer_len ?
					kmer_len : pkext->ctx_distance;
				kmer_id_to_seq(ctxs[c], ctxlen, kseq);
				kseq[ctxlen] = '\0';
				printf("  >> [%d] %s / %d\n", c, kseq, pkext->ctx_distance);
			}
		}
		else {
			kbits_t tpacked = trans_flag2packed(pk->trans_flag) >> 3;
			for (int i = 0; i < pk->n_trans; ++i, tpacked >>= 2) {
				kbits_t base = tpacked & 3;
				kbits_t k1 = kmer_reverse_complement(
						(kid >> 2) | (base << shift), kmer_len);
				kmer_t *pk1 = kptr_locate(d, k1);

				if (pk1->has_multi_copies == 0) continue;

				kmer_ext_t *pkext1 = &d->nonuniq_kmers[pk1->index];
				int n_ctx = pkext1->n_ctx;

				kbits_t *ctxs = &d->kmer_ups_contexts[pkext1->ctxoff];

				for (int c = 0; c < n_ctx; ++c) {
					int ctxlen = pkext1->ctx_distance > kmer_len ?
						kmer_len : pkext1->ctx_distance;
					kmer_id_to_seq(ctxs[c], ctxlen, kseq);
					kseq[ctxlen] = '\0';
					printf("  >> [%c] %s / %d\n", BASES[base], kseq, pkext1->ctx_distance);
				}
			}
		}
	}

	if (pk->trans_flag > 0) {
		double *exp_trans_p = hmmsbg_locate_exp_trans(d, pk);
		double *cur_tp = hmmsbg_locate_trans_p(d, pk, &uniq_tp);

		int n_ctx = pk->has_multi_copies ? 
			d->nonuniq_kmers[pk->index].n_ctx : 0;

		if (pk->label == 0) {
			printf(" * Trans. prob. (prev. iter): \n");
			for (int c = 0; c <= n_ctx; ++c) {
				printf("  >> ");
				kbits_t tpacked = trans_flag2packed(pk->trans_flag) >> 3;
				for (int i = 0; i < pk->n_trans; ++i, tpacked >>= 2) {
					printf("%c (%6.3e) ", BASES[tpacked & 3], 
							EXP(cur_tp[c * pk->n_trans + i]));
				}
				printf("\n");
			}

			printf("\n * Expected trans: \n");
			for (int c = 0; c <= n_ctx; ++c) {
				hmmsbg_multictx_lbl2_exp_trans(d, pk, c, &jnt_exp_trans[0]);
				printf("  >> ");
				kbits_t tpacked = trans_flag2packed(pk->trans_flag) >> 3;
				for (int i = 0; i < pk->n_trans; ++i, tpacked >>= 2) {
					printf("%c <%g / %g> ", BASES[tpacked & 3], 
							exp_trans_p[c * pk->n_trans + i], 
							exp_trans_p[c * pk->n_trans + i] + jnt_exp_trans[i]);
				}
				printf("\n");
			}
			printf("\n");
		}
		else {
			kbits_t tpacked = trans_flag2packed(pk->trans_flag) >> 3;
			uint32_t n_ctx = pk->has_multi_copies ? 
				d->nonuniq_kmers[pk->index].n_ctx : 0;

			int par_off = 0;
			printf("\n * Trans. prob. (prev. iter): \n");
			for (int i = 0; i < pk->n_trans; ++i, tpacked >>= 2) {
				uint32_t n_ctx_i = (n_ctx >> (i * PMR_MAX_CTX_BITS)) & PMR_MAX_NCTX;	
				for (int c = 0; c <= n_ctx_i; ++c) {
					printf("  >> %c (%6.3e) \n", BASES[tpacked & 3], 
							EXP(cur_tp[par_off + c]));
				}

				par_off += (n_ctx_i + 1);
			}

			printf("\n * Expected trans: \n");
			par_off = 0;
			tpacked = trans_flag2packed(pk->trans_flag) >> 3;
			for (int i = 0; i < pk->n_trans; ++i, tpacked >>= 2) {
				uint32_t n_ctx_i = (n_ctx >> (i * PMR_MAX_CTX_BITS)) & PMR_MAX_NCTX;	

				kbits_t k1 = kmer_reverse_complement(
						(kid >> 2) | ((tpacked & 3) << shift), kmer_len);
				kmer_t *pk1 = kptr_locate(d, k1);
				int tidx = base_to_transidx(pk1->trans_flag,
						complementary_base(kid & 3));

				for (int c = 0; c <= n_ctx_i; ++c) {
					double exp_trans_k1 = *(hmmsbg_locate_exp_trans(d, pk1) + 
							(c * pk1->n_trans) + tidx);

					printf(" >> %c <%g / %g> \n", BASES[tpacked & 3], 
							exp_trans_p[par_off + c],
							exp_trans_p[par_off + c] + exp_trans_k1);
				}

				par_off += (n_ctx_i + 1);
			}
		}
	}
	else {
		printf(" x No outgoing transitions.\n");
	}
}
