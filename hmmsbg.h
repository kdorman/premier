#ifndef __PMR_HMMBSG_H__
#define __PMR_HMMBSG_H__

extern "C" {
#include "em.h"
#include "kmer.h"
#include "premier.h"
#include "mstep.h"
}

#ifdef __cplusplus
#define EXTERNC extern "C"
#else
#define EXTERNC
#endif

static inline size_t hmmsbg_nbhd_alloc_size(int n) 
{
	return (((n * (sizeof(double) + sizeof(kbits_t))) << 1) +
			(n * sizeof(intptr_t) * 3) +
			(((BITS_TO_U64(n << 1)) + 2 * (BITS_TO_U64(n << 3)) + 
			  (BITS_TO_U64(n << 3))) << 3) + (5 * n * sizeof(int64_t)));
}

static inline double *hmmsbg_locate_trans_p(data *d, kmer_t *pk, double *uniq_tp)
{
	return (pk->unique_trans && pk->has_multi_copies == 0) ? uniq_tp : 
			d->transition_p + d->nonuniq_kmers[pk->index].poff;
}


static inline double *hmmsbg_locate_exp_trans(data *d, kmer_t *pk)
{
	if (pk->has_multi_copies || pk->n_trans > 1)
		return d->exp_trans_p + d->nonuniq_kmers[pk->index].poff;
	else
		return d->exp_trans_p + (pk->index + d->uniq_kmer_param_off);
}

static inline uint32_t hmmsbg_initial_index(data *d, kmer_t *pk)
{
	uint32_t idx = pk->index;

	if (pk->has_multi_copies) {
		idx += d->multictx_kmer_init_off;
	}
	else if (pk->n_trans == 1) {
		idx += d->uniq_kmer_init_off;
	}

	return idx + (pk->n_trans > 0);
}

static inline double* hmmsbg_initial_p(data *d, kmer_t *pk)
{
	uint32_t index = hmmsbg_initial_index(d, pk);
	if (pk->has_multi_copies && pk->label == 0) {
		return &d->multictx_initial_p[(uint32_t) d->avg_initial_p[index]];
	}
	return &d->avg_initial_p[index];
}

void hmmsbg_e_step(data *d, read_t *read, void *fdata);
void hmmsbg_viterbi(data *d, read_t *read, void *fdata);
void hmmsbg_m_step(data *d);
void hmmsbg_update_expected_vals(data *d, state_nbhd_t *curr_nbhd, 
		int t, double *t_exp_counts, double *t_exp_trans, 
		double *t_base_emit_exp, double *t_qual_emit_exp,
		int kmer_len, int ctx_shift, int qmax, 
		const char *rseq, const char *conv_q,
		double inv_sum_exp_count);
uint32_t hmmsbg_backtrack(data *d, kbits_t kid, int max_paths, int max_len,
		std::vector<kbits_t> &vctx, 
		std::unordered_set<std::string> &k_ups_ctxs);

EXTERNC void HMMSBG(data *d);

#undef EXTERNC

#endif
