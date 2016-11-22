#ifndef __PREMIER_H__
#define __PREMIER_H__

#include <stdint.h>
#include <math.h>

#include "kmerdict.h"
#include "mempool.h"


/* 20 = 4 (true base) * 5 (observed base, ACTGN) */
#define PMR_BASE_EMIT_COMBS 20
#define PMR_UNIF_BASE_EMIT_P -1.6094379124341002818 
#define PMR_LINE_BUFFER_SIZE 1024
#define PMR_READ_BUFFER_SIZE 1073741824 //(1<<30) 

#define PMR_NO_DEP_EMISSION 1
#define PMR_READ_EXTENSION_TWO_STRANDS 0

#define PMR_MAX_CTX_BITS 8
#define PMR_MAX_NCTX ((1U << PMR_MAX_CTX_BITS) - 1)

#define PMR_HMMSBG 0
#define PMR_EXTEND_READS 0

#define PMR_MODE_ORACLE 1
#define PMR_MODE_ERROR_CORRECTION 2
#define PMR_MODE_PERFECT_ORACLE 3

/* --- some precompiled switches --- */

/* faster, but less numerically stable */
#define PMR_USE_LOG_ADD 0
#define PMR_ENABLE_DIAGNOSTICS 0
#define PMR_JUMP 0

/* accept phred quality scores as probabilities of error,
 * see --phred for handling of quality score 2 */
//#define PMR_ALTERNATIVE_EMISSION 1

// supplement transitions using the true de Bruijn graph in the Oracle mode.
#define PMR_ORACLE_SUPPLEMENT_TRANS 0

#define PMR_BATCH_GROUP_SIZE 10000

//#define PMR_DUMP_NBHD_SIZE 0
//#define PMR_APPROX_HMM 1

/* use the most basic form of emission distribution:
 * no quality score, and a simple error rate parameter for the base emission */
//#define PMR_USE_PRIMITIVE_EMISSION 1

/* use standard penalty function instead of the penalty that combines coverage
 * of both strands.
 */
/* combine coverage from original reads and their reverse complements
 * (artificially created) to initialize the transition parameters */
//#define PMR_USE_JOINT_INITIALIZATION 1 
#define PMR_COMBINE_TRANS_FLAG 1

/* The length, measured by units of uint64_t type, to store kmers in
 * 2-bit numeric representation.
 * By default, a uint64_t type is used which permits the maximum of 
 * kmer length to be 31 (highest two bits reserved).
 * If a kmer size beyond 31 is needed, adjust this value by specifying
 * -DPMR_KMER_WIDTH=d option to the compiler.
 * */
#ifndef PMR_KMER_WIDTH
#define PMR_USE_SCALAR_KMER_TYPE 1
#define PMR_KMER_WIDTH 1
#endif

// log(.Machine$double.eps)
#define PMR_NBHD_THIN_THRESHOLD 36.04365

/* NOTE: not in regular ACGT order for computation convenience */
enum bases {BASE_A, BASE_C, BASE_T, BASE_G, BASE_N};
enum partial_decoding_lvls {PDEC_DISABLED, PDEC_CONSERVATIVE, PDEC_AGGRESSIVE};

#if PMR_KMER_WIDTH == 1
typedef uint64_t kbits_t;
#else
typedef struct kbits_s {
	uint64_t u64[PMR_KMER_WIDTH];
}
#endif

/* used for a substring shorter than 32 bases.
 * when PMR_KMER_WIDTH == 1, this is equivalent to kbits_t.
 * the equivalence disqualifies, however, when longer kmer width
 * is chosen by user. */
typedef uint64_t ksubstr_t;

typedef struct kmer_s kmer_t;
typedef struct kmer_nbhd_s kmer_nbhd_t;
typedef struct kmer_ext_s kmer_ext_t;
typedef struct kmer_counter_s kmer_counter_t;
typedef struct state_nbhd_s state_nbhd_t;
typedef struct _options options;

typedef union {
	double dbl;
	uint64_t u64;
} dbl_t;

typedef struct {
	int size;
	int used;
	kmer_t **kmer_refs; /*< pointer to these reference kmers */
	/* hash_t *raw_hashes; */
} nbhd_t;

typedef struct {
	kbits_t masked_id;
	kbits_t pmask;	/* position_mask */
} nbhd_key_t;

typedef struct {
	int length;	/*< short read length */
	size_t id;	/* read identifier */
	size_t batch_id;

	char *identifier;
	char *sequence;	/*< literal read sequence */
	char *qscore;	/*< quality score associated with a read */
} read_t;

#if 0
struct kmer_s {
#ifdef PMR_ENABLE_DIAGNOSTICS
	/* expected count computed in previous E-step */
	double prev_exp_count;
#endif
	double count;
	double exp_count;
	double init_p;
	double *transition_p;
	double exp_trans_p[4];

	kbits_t id;
	uint32_t trans_flag : 4;
	uint32_t unique_trans : 1;
	uint32_t alloc_nbhd_size : 4;
	uint32_t emit_norm_flag : 4;
	uint32_t n_neighbors : 19;

	void *data;
};
#endif

/* a packed kmer struct fit in a 64bit pointer */
struct kmer_s {
	uint32_t trans_flag: 4;
	uint32_t unique_trans: 1;
	uint32_t n_trans: 3;
	/* partition transition probabilities into two sets, which are reverse
	 * complements to each other. The set is labeled by this bit. */
	uint32_t label: 1;
	uint32_t prev_n_trans: 3;
	uint32_t dirty: 1; /* has been partially updated in CM-step */
	/* trans. flag from the reverse complementary strand */
	uint32_t alt_trans_flag: 4; 
	/* valid in the sense of having coverage above threshold */
	// uint32_t valid : 1;
	uint32_t coverage: 10;
	/* is this kmer *ONLY* observed at the 3' end of a read ? */
	uint32_t last_kmer : 1; 
	uint32_t ambig_labeled : 1;
	uint32_t has_multi_copies : 1;
	uint32_t single_stranded: 1;
	uint32_t reserved: 1;
	/* for a kmer with unique_trans == 0, `index` corresponds to both the 
	 * index of element in array `data.nonuniq_kmers`;
	 * if unique_trans == 1, it locates the conditional expected transition
	 * for the unique transition in `data.exp_trans_p`.
	 */
	uint32_t index;
};

struct kmer_ext_s {
	/* offset to find the corresponding parameters */
	uint32_t poff;  
	uint32_t ctxoff;
	// for label-1 kmers, n_ctx is a single number
	// for label-2 kmers, n_ctx consists of 4, 8-bit numbers (255 contexts max)
	uint32_t n_ctx;
	// for label-2 kmers, ctx_distance consists of 4, 8-bit numbers 
	// (max. dist=255) 
	uint32_t ctx_distance;
	//uint32_t reserved : 6;
};

/* FIXME: consider rearrange memory to improve locality / optimize
 * cache usage */
struct state_nbhd_s {
	int size;
	int n_uniq_sfx;
	int n_trans;
	int max_hd;
	kmer_t **kmer_ptrs;
	double **kmer_trans_prob;

	uint64_t *ba_kmer_trans_flag;
	uint64_t *ba_distinct_sfx;
	uint64_t *ba_hamming_dist;
	uint64_t *ba_pfx_sfx_flags;

	int64_t *suffixes_order;
	kbits_t *states_sorted_pfx;
	kbits_t *states_sorted_sfx;

	double *alphas;
	double *betas;
	double *emit_dens;
	int64_t *base_emit_ctx;
	int64_t *qual_emit_ctx;
	uint64_t *kmer_ups_ctx;
	int64_t *jump_ahead_countdown;
};

typedef struct {
	dict *kmer_dict;
	dict *nbhd_dict;

	options *opt;

	int max_nbhd_size;
	int max_read_length;
	int n_hd_windows;
	int n_emit_windows;
	uint64_t n_reads;

	size_t n_qual_context;
	size_t qual_ctx_radix;
	size_t n_base_context;
	size_t base_ctx_radix;
	size_t init_ctx_radix;

	size_t n_uniq_trans_kmer;
	size_t n_zeroed_out_kmer;
	size_t n_composite_trans;
	size_t n_nonzero_kmers;

	/* --- Baum-Welch --- */
	double loglik;
	double total_kmer_exp_counts;

	/* - transition - */
	double *transition_p;
	double *exp_trans_p;
	uint64_t n_total_trans;
	uint64_t n_nonuniq_trans;
	uint64_t n_nonuniq_kmers;

	size_t n_first_kmers;
	uint32_t model_complexity;
	uint32_t model_obs_complexity; // complexity under strict uniqueness

	double prop_nonuniq_kmers;

	kmer_ext_t *nonuniq_kmers;

	/* - emission - */
#if PMR_USE_PRIMITIVE_EMISSION
	double base_error_rate;
	double error_rate_emit_exp[2];
#else
	double *base_emission_p;
	double *qual_emission_p;
	double *base_emit_exp;
	double *qual_emit_exp;

	double *qual_err_rate;
#endif

	/* initial distribution */
	double *avg_initial_p;
	double *multictx_initial_p;

	/* probability that a kmer cannot extend on the de Bruijn graph,
	 * due to loss of coverage, end of genome, or being an error */
	double *prob_kmer_extension;

	/* two offsets used to compute the actual index in the arrays above */
	int uniq_kmer_init_off;
	int uniq_kmer_param_off;
	int multictx_kmer_init_off;

	// # of kmer transitions removed, when all expected transitions for a
	// certain kmer is below the penalizing threshold (\rho).
	int n_trans_removed_below_rho;

	/* --- Precomputed tables --- */
	char **n_choose_k;	/*< binomial coefficients (n choose k) */
	char **position_perms;	/*< permutation of all {k \choose d} erroneous positions */
	uint64_t *position_masks;	/*< masks for {k \choose d} erroneous positions. */
	uint64_t *error_perms;
	//kbits_t *transition_masks;
	double *phred_to_prob; /*< table that converts Phred quality score to error probability. */ 
	int *read_start_positions;
	double *read_em_time_cov;
	uint8_t *read_dmax;
	uint64_t *kmer_ups_contexts;

	kmer_nbhd_t *preconstructed_nbhds;

	/* --- output files --- */
	FILE *log_fp;
	FILE *errors_fp;
	FILE *params_fp;
	FILE *sdump_fp;
#ifdef PMR_ENABLE_DIAGNOSTICS
	FILE *viterbi_diag_fp;
#endif

#ifdef PMR_DUMP_NBHD_SIZE
	FILE *nbhd_size_fp;
#endif

	/* --- mmap --- */
	char *fastq_mmap;
	size_t fastq_size;
	int fastq_fd;

	mempool_t *mp;
	mempool_t *tmpmp;
} data;

struct _options {
	int kmer_length;
	int n_threads; 

	int preconstructed_nbhd_hd;
	int max_hamming_dist; 

	size_t max_approx_nbhd_size;

	int mode;
	int pagesize;

	int qual_score_offset;
	int max_qual_score;
	/* binning quality scores into more coarse discrete scales */
	int qscore_bin_width;

	int n_max_kmer_contexts;
	int max_backtrack_dist;

	int base_emit_context;
	int qual_emit_context;

	int max_iter;
	double tol;
	double error_rate;

	double kcov_q75;
	double kcov_q90;

	int hd_window_size;
	int emit_window_size;

	uint32_t output_params : 1;
	uint32_t random_param_init : 1; 
	uint32_t disable_penalized_em : 1;
	uint32_t enable_params_output : 1;
	uint32_t filename_append_timeinfo : 1;
	uint32_t filename_append_params : 1;
	uint32_t penalize_init_params : 1;
	uint32_t distinct_q2 : 1;
	uint32_t enable_error_watch : 1;
	uint32_t partial_decoding_lvl : 2;
	uint32_t init_penalty_func: 2;
	uint32_t dump_strandness: 1;
	uint32_t em_penalty_func: 2;
	uint32_t auto_kmer_length : 1;
	uint32_t auto_penalty_eta : 1;
	uint32_t viterbi_revcomp : 1;
	uint32_t is_tty : 1;
	uint32_t mstep_debug_console: 1;
	uint32_t diagnose_mode : 1;
	uint32_t output_ordered : 1;
	uint32_t adaptive_hamming_dist : 1;
	uint32_t update_qual_score : 1;
	uint32_t qscore_as_phred: 2;
	uint32_t dummy : 7; 

	/* penalized EM */
	double penalty_eta;
	double penalty_gamma;

	const char *fastq_file;
	char *output_prefix;
	char *ground_truth_file;
	char *qscore_pmf_file;
	char *params_data;
};

typedef struct {
	double gamma;
	double *eta;
	double *xi;
	int32_t signs;
	int n_xi;
} mstep_data_t ;

typedef struct {
	kbits_t id;
	kbits_t mut_flag;
} mut_kmer_t;

struct kmer_nbhd_s {
	int n;
	kmer_t **nbhd;
};

typedef struct {
	char identifier[8]; /* PMRPARAM */
	uint32_t kmer_length: 6;
	uint32_t max_qscore: 8;
	uint32_t max_read_length: 10; 
	uint32_t max_hamming_dist: 4;
	uint32_t reserved : 4;
	double penalty_eta;
	double penalty_gamma;
	double loglik;
	uint16_t base_ctx_radix;
	uint16_t qual_ctx_radix;
	uint32_t n_uniq_trans_kmer;
	uint32_t n_multi_trans_kmer;
	uint64_t off_trans_p;
	uint64_t off_emit_p;
	uint64_t off_kmers;
} dumpfile_header_t;

typedef struct {
	kbits_t id;
	double exp_count;
} dumpfile_kmer_t ;

void hmm_build_1st_nbhd(data *d, read_t *read, void *fdata);
void hmm_determine_dmax(data *d, read_t *read, void *fdata);
void hmm_adaptive_dmax(data *d, read_t *read, void *fdata);
void HMMSBG(data *d);
void run_oracle(data *d);
void oracle_init_em(data *d);

#endif
