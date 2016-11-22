#ifndef __DATA_H
#define __DATA_H

#include <stdlib.h>
#include <stdint.h>
#include <pthread.h>
#include <time.h>

#include "dict.h"
#include "iodata.h"

#define LINE_BUFFER_SIZE 128

enum {FASTQ_INIT, FASTQ_ID, FASTQ_SEQ, FASTQ_QID, FASTQ_Q};
enum {NU_A = 0, NU_C, NU_T, NU_G};

typedef pthread_mutex_t mutex_t;
typedef double DOUBLE;

typedef struct {
	int kmer_length;
	int use_log_calculations;
	int n_processor;

	char *fastq_file;
} options;

typedef struct {
	dict *kmer_dict;
	DOUBLE ll, pll;

	dictsize_t num_unique_kmer;	/* number of unique oligomers observed */
	options *opt;
} data;

typedef struct {
	uint32_t id;	/* k-mer unique id. each two bits are used to represent a nucleotide, and 16mer at maximum */
	uint32_t count;
	mutex_t lock;	/* read/write lock */

	DOUBLE stationary_p;	/* stationary probability for given kmer */
	DOUBLE transition_p[4];	/* probility of current kmer transitioning into ACGT */

	DOUBLE *gamma;  /* used for Baum-Welch algorithm */
	DOUBLE **xi;
} kmer_t;

typedef struct {
	uint32_t id;	/* same as kmer id */
	char *qscore;	/* quality score associated with kmer; might be needed by emission distribution */

	/* NOTE: since uint32_t can hold at maximum a 16-mer, we can use uint64_t
	 * to hold flags for all 16*4 neighbors that have distance up to 1 to a
	 * given 16-mer. 
	 *
	 * Therefore, the code works for kmers at most 16bp long, and only works
	 * for neighborhood with at most 1 hamming distance.
	 *
	 * More complex error patterns require advanced data structures.
	 */
	uint64_t nbhd_mask;	/* each bit represents if a corresponding kmer exists */
	kmer_t **nbhd_kmers;	/* pointers to correponding kmer struct */
	DOUBLE *alpha;	/* forward algorithm */
	DOUBLE *beta;	/* backward algorithm */
} readnode_t;

/* NOTE: do we need to store every read in the memory?
 * One immediate problem that we will encounter by storing all reads in memory
 * is that for large genomes and high coverages, the data simply can't scale
 * into the memory.
 *
 * And hence we'd better read the fastq file every time, or use some other
 * compressed format to (temporarily) store the sequences.
 */

typedef struct {
	int length;
	char *sequence;
	char *qscore;	/* quality score associated with a read */

	int n_oligomer;
	readnode_t **oligomers;
} read_t;

/* function declarations */
read_t *read_create();
void read_destroy(read_t *tmp_read);
int read_parse(data *d, read_t *tmp_read);

int kmer_incr_count_by_id(data *d, uint32_t kid);
kmer_t *kmer_create(uint32_t kid);
void kmer_destroy(kmer_t *kmer);

uint32_t kmer_seq2id(const char *s, int size);


#endif
