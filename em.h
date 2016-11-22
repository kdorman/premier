#ifndef __PREMIER_EM_H__
#define __PREMIER_EM_H__

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <sys/time.h>
#include <signal.h>

#include "bitarray.h"
#include "premier.h"

#define MSTEP_FLIP_LABEL 0x1
#define MSTEP_NO_EMISSION 0x2

/* FIXME: consider using lookup table instead */
#define nbhd_alloc_size(n)	\
		(((n * (sizeof(double) + sizeof(kbits_t))) << 1) +  \
		 (n * sizeof(intptr_t) * 3) +                          \
		 (((BITS_TO_U64(n << 1)) + 2 * (BITS_TO_U64(n << 3)) + \
		  (BITS_TO_U64(n << 3))) << 3) + (3 * n * sizeof(int64_t)))

void EM(data *d, struct timeval *ts);
void hmm_e_step_wrapper(data *d, read_t *read, void *fdata);
void hmm_e_step(data *d, read_t *read, void *fdata, read_t *argmax_read);
void hmm_probe_nbhd_size(data *d, read_t *read, void *fdata);

double hmm_compute_alpha(kbits_t ups_trans_packed, int ups_kmer_idx,
		kmer_t **t_kmer_ptrs, double *t_state_alphas, double **t_kmer_trans_p, 
		double emit_dens, int dns_base);
double hmm_compute_beta(kbits_t dns_trans_packed,
		int dns_kmer_idx, kbits_t _tf, double *tnext_state_betas, 
		double *kmer_trans_p, double *kmer_exp_trans, 
		int64_t *tnext_sfx_order, double *emit_dens, double alpha);
double hmm_compute_delta(kbits_t ups_trans_packed, 
		int ups_kmer_idx, kmer_t **t_kmer_ptrs, double *t_state_deltas, 
		double **t_kmer_trans_p, double emit_dens, int dns_base, double *delta);
double hmm_compute_emit_dens(kbits_t st, kbits_t xt, int yt, kbits_t enf,
		double *base_emit_p, double *qual_emit_p);
double hmm_compute_kmer_emission(data *d, kbits_t obs_kid, 
		kbits_t obs_n_flag, kmer_t *pk, int kmer_len, const char *conv_q);

void hmm_correct_errors(data *d, read_t *read, int tmax,
		kbits_t st, kbits_t xt, kbits_t xt_n_flag, int kmer_len);
void hmm_correct_one_error(data *d, read_t *read, int t,
		kbits_t st, kbits_t xt);

char *convert_quality_scores(const char *qscore, int len, int offset,
		int bwidth, mempool_t *mp);

void hmm_read_starting_pos(data *d, read_t *read, void *fdata);
void hmm_dump_strandness(data *d, read_t *read, void *fdata);

void hmm_kcov_quantile(data *d, read_t *read, void *fdata);

void hmm_setup_nbhd_ptrs(state_nbhd_t *nbhd, int n, 
		char *pmem);

void hmm_init_e_step(data *d);
void hmm_init_model_params(data *d);

int hmm_load_preconstructed_nbhd(data *d, int rid,
		kbits_t obs_1st_kid, kbits_t obs_n_flag, const char *conv_q, 
		state_nbhd_t *first_nbhd, int kmer_len, int compute_emit_dens, 
		dict *kdict, mempool_t *mp, const double *uniq_tp);

int hmm_count_distinct_suffixes(state_nbhd_t *nbhd, 
		int *actg_trans_count, kbits_t t_dmax, kbits_t obs_kid,
		kbits_t obs_n_flag, kbits_t obs_next_base, double *unique_trans_p);

static inline int calc_bctx(kbits_t st, kbits_t hd, int shift, int single_stranded)
{
	// clzl: number of leading zeros before the most significant bit
#if PMR_NO_DEP_EMISSION
	return (st >> shift);
	//return (st >> shift) + 4 * single_stranded;
#else
	if (hd == 0) return (st >> shift);

	int _msig = 64 - __builtin_clzl(hd);
	if (_msig > 5) raise(SIGINT);
	return 4 * _msig + (st >> shift);
#endif
}

static inline int calc_qctx(kbits_t st, kbits_t xt, char *yt_pfx, int k)
{
	uint32_t n_q2 = 0;
	int shift = (k - 1) << 1;
	// FIXME: this for loop can be optimizied, and only compute once.
	for (int p = 0; p < k - 1; ++p) n_q2 += (yt_pfx[p] == 2);
	if (n_q2 == 0) return 6 * ((st >> shift) == (xt >> shift));

	int _msig = 32 - __builtin_clz(n_q2);
	// 6 is the number of distinct values of n_q2:
	// 0, 1, 2, 4, 8, 16
	return 6 * ((st >> shift) == (xt >> shift)) + _msig;
}

#endif
