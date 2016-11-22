#ifndef __PREMIER_MSTEP_H__
#define __PREMIER_MSTEP_H__

#include <stdint.h>
#include <math.h>
#include "premier.h"
#include "numeric.h"

#ifdef __POPCNT__
#include <nmmintrin.h>
#endif

typedef struct tflag_s tflag_t;

double cmstep_estimate_trans_p(data *d, kmer_t *pk,
		double *tp, double *exp_trans_p, 
		double *initial_p, uint64_t trans_flag, double thresh, 
		mstep_data_t *md, tflag_t *tf);

void cmstep_reciprocal_trans_p(data *d, kmer_t *pk, double mu, double *recip_tp);
uint32_t mstep_reciprocal_exp_trans(data *d, kmer_t *pk, double *recip_exp_trans);

void cmstep_penalized_transition(data *d, int pen2, int flip_label);
void mstep_emission(data *d);
void mstep_alt_emission(data *d);
void mstep_initial(data *d);

void mstep_strand_specific_errors(data *d);


/* two penalty version */
double mstep_pen1_lambda_support(void *fdata);
double mstep_pen1_lagrange_cstr_func(double lambda, void *fdata);
double mstep_pen1_lagrange_cstr_derv(double lambda, void *fdata);
int32_t mstep_pen1_determine_signs(double lamb_lb, void *fdata);

double mstep_pen2_lambda_support(void *fdata);
double mstep_pen2_lagrange_cstr_func(double lambda, void *fdata);
double mstep_pen2_lagrange_cstr_derv(double lambda, void *fdata);
int32_t mstep_pen2_determine_signs(double lamb_lb, void *fdata);

typedef double (*msp_lcstr_func)(double lambda, void *fdata);
typedef double (*msp_lcstr_derv)(double lambda, void *fdata);

double mstep_newton(double (*fx)(double x, void *data), 
		double (*fderv)(double x, void *data), double x0, double eps,
		int maxiter, void *fdata);
double mstep_bisection(double (*fx)(double x, void *data), double a, double b, 
		double eps, void *fdata);
double cmstep_try_bisection(mstep_data_t *d2, double lamb_lb,
		int32_t pos_signs, int *signs, 
		msp_lcstr_func plfunc, msp_lcstr_derv plderv);

void mstep_extend_all_kmers(data *d);

struct tflag_s {
	uint32_t trans_flag;
	int n_nonzero_trans;
};

static void renormalize_trans_p(int n, double *tp)
{
	double sum_tp = 0.0;
	for (int j = 0; j < n; ++j) {
		sum_tp += EXP(tp[j]);
	}

	for (int j = 0; j < n; ++j) {
		tp[j] -= LOG(sum_tp);
	}
}

#endif
