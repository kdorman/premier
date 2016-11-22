#ifndef __DATA_H
#define __DATA_H

#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <time.h>

#include "premier.h"
#include "kmer.h"
#include "numeric.h"

#define PREALLOC_NBHD_SIZE 16 
#define LINE_BUFFER_SIZE 1024

#define LOG_ONE_FOURTH -1.386294361119891

enum {NU_A = 0, NU_C, NU_T, NU_G};

/* function declarations */
void allocate_trans_params(data *d, int halve_threshold);
int compute_kmer_hamming_nbhds(data *d);

void kmer_id2seq(kbits_t kid, int size);

/* --- precomputed tables --- */
int compute_error_positions(data *d);
int compute_binomial_coefficients(data *d);
int compute_all_errors(data *d);
int compute_phred_to_probability(data *d);

void adjust_kmer_strandness(data *d);
void adjust_kmer_labeling(data *d);
void readjust_kmer_labeling(data *d);

typedef struct {
	kbits_t pmask;
	int mask_shift;
} mask_t;

#endif
