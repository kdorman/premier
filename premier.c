#define _GNU_SOURCE

#include <sys/time.h>
#include <stdio.h>
#include <unistd.h>
#include <limits.h>
#include <getopt.h>
#include <time.h>
#include <stdbool.h>
//#include <signal.h>

#include "premier.h"
#include "read.h"
#include "kmer.h"
#include "data.h"
#include "iodata.h"
#include "em.h"
#include "numeric.h"
#include "nbhd.h"

static int make_options(options **opt, int argc, char **argv);
static int make_data(data **d, options *opt);
static int destroy_data(data *d);
static void sample_probing_reads(data *d, int *rsel, int N, int n);
static void print_usage(const char *basename);
//static void sigint_handler(int signum);

static void resize_kspectrum(data *d, void *fdata) {
	dict_try_resize(d->kmer_dict);
}

//static data *__d;

static char *pmr_usage = "Usage: premier <command> [options] <fastq_file>\n\n"
	"COMMANDS\n"
	"  correct               Performs standard error correction. \n"
	"  oracle                Evaluate performance of a k,d-local Oracle. \n\n"
	"  poracle               Evaluate performance of a k,d-local perfect Oracle. \n\n"
	"OPTIONS\n"
	"  -k K[,auto]           Kmer length used for error correction. \n"
	"                        Append \",auto\" after K to enable automatic\n"
	"                        selection of kmer length.\n"
	"  -d DMAX               Maximal number of errors allowed in each kmer\n"
	"  Model estimation:\n"
	/*
	"  --no-penalized-em     Turn off penalized estimation.\n"
	"  --no-penalized-init   Do not impose penalty on HMM parameter\n"
	"                        initialization.\n"
	*/
	"  -T --thresh=THRESH    [Penalized EM] Threshold for the approximate\n"
	"                        L0-type penalty.\n"
	"  -E --eta=ETA          [Penalized EM] Determing the extent of penalty.\n"
	"                        (It is more intuitive to use -T than this.)\n"
	"  -G --gamma=GAMMA      [Penalized EM] Determing how well the penalty\n"
	"                        approximates an L0-type penalty.\n"
	"  -R --rel-tol=RELTOL   Convergence criterion in relative tolerance\n"
	"  --error-rate=ER       Estimated per-base average error rate, which is\n"
	"                        used to initialize the emission parameters.\n"
#ifdef PMR_ALTERNATIVE_EMISSION
	"  --phred={1, 2}        Interpret all (2) quality scores as error probabilities\n"
	"                        or all but 2 (1) as probabilities: error probability\n"
	"                        assigned to 2 is given by --error-rate argument\n"
#endif
	/*
	"  Emission distribution:\n"
	"  -B CTX_LEN            Length of upstream context that base emissions\n"
	"                        depend on.\n"
	"  -Q CTX_LEN            Context length for quality score emission\n\n"
	"  Model loading/dumping:\n"
	"  -p                    Dump HMM parameters to a local file. Output is \n"
	"                        suppressed if `-o` option is not specified.\n"
	"  -P PARAM_FILE         Load previously dumped parameters in PARAM_FILE\n\n"
	*/
	"  Results & output:\n"
	"  -o OUT_PREFIX         Prefix for all output results, logs, diagostics\n"
	"                        and HMM parameters\n"
#ifdef PMR_ENABLE_DIAGNOSTICS
	/*
	"  --diagnose={0, 1}     Generate a full (1) memory snapshot of Viterbi\n"
	"                        decoding or a lite (0) version, with data from\n"
	"                        decoded and observed sequences only.\n"
	*/
#endif
	"";

static int cmp_dbl(const void *a, const void *b)
{
	const register double *pa = a;
	const register double *pb = b; 

	return *pa > *pb;
}


int main(int argc, char **argv)
{
	data *d = NULL;
	options *opt = NULL;

	struct timeval ts, te;
	time_t rt;

	if (make_options(&opt, argc, argv))
		exit(EXIT_FAILURE);

	/* record timestamp */
	gettimeofday(&ts, NULL);

	//signal(SIGINT, sigint_handler);

	if (make_data(&d, opt)){
		/* error processing here */
		exit(EXIT_FAILURE);
	}

	if (d->log_fp != NULL) {
		fprintf(d->log_fp, "# Call: ");
		for (int i = 0; i < argc; ++i) {
			fprintf(d->log_fp, "%s ", argv[i]);
		}
		fprintf(d->log_fp, "\n");
	}

	if (d->opt->mode == PMR_MODE_ORACLE) {
		// build 1st neighborhood
		oracle_init_em(d);
	}
	if (d->opt->mode == PMR_MODE_PERFECT_ORACLE) {
		run_oracle(d);
	}
	else {

		/* data prepared, proceed to em */
		EM(d, &ts);

#if PMR_HMMSBG
		HMMSBG(d);
#endif
	}

	/* timing info */
	time(&rt);
	gettimeofday(&te, NULL);

	if (d->log_fp != NULL){
		fprintf(d->log_fp, "# PREMIER finished at %s\n", ctime(&rt));
		fprintf(d->log_fp, "# Total seconds elapsed: %.3f\n",
				(0.0 + te.tv_sec + (double) te.tv_usec/1e6) - 
				(0.0 + ts.tv_sec + (double) ts.tv_usec/1e6));
		fclose(d->log_fp);
	}

	if (d->errors_fp != NULL)
		fclose(d->errors_fp);

	if (d->params_fp != NULL)
		fclose(d->params_fp);
	/*
	if (d->ll_fp != NULL)
		fclose(d->ll_fp);
	*/

	/* cleanups */
	destroy_data(d);

	return 0;
}

int make_data(data **d, options *opt)
{
	int kmer_len = opt->kmer_length;

	data *pd;
	time_t rt;

	*d = (data *) malloc(sizeof(**d));
	//__d = d;

	if (*d == NULL){
		return message(stderr, __FILE__, __LINE__,
				ERROR_MSG, MEMORY_ALLOCATION, NULL);
	}

	pd = *d;
	pd->opt = opt;

	/* initialize some default values */
	pd->max_read_length = 0;
	pd->total_kmer_exp_counts = 0.0;

	pd->errors_fp = NULL;
	pd->log_fp = NULL;
	pd->params_fp = NULL;
#ifdef PMR_ENABLE_DIAGNOSTICS
	pd->viterbi_diag_fp = NULL;
#endif
	pd->preconstructed_nbhds = NULL;
	pd->read_start_positions = NULL;
	pd->read_dmax = NULL;
	pd->n_hd_windows = 0;
	pd->read_em_time_cov = NULL;
	pd->n_reads = 0;

	pd->transition_p = NULL;
	pd->exp_trans_p = NULL;
	pd->avg_initial_p = NULL;
	pd->prob_kmer_extension = NULL;
	pd->n_total_trans = 0;

	if (pd->opt->output_prefix != NULL) {
		char *out_fname = io_cat_output_filename(pd->opt, "", "fastq");	
		if (out_fname != NULL) {
			pd->errors_fp = fopen(out_fname, "w+");
			free(out_fname);
		}

		out_fname = io_cat_output_filename(pd->opt, "", "log");
		if (out_fname != NULL) {
			pd->log_fp = fopen(out_fname, "w+");
			/* set the log file to line buffer mode */
			setvbuf(pd->log_fp, NULL, _IOLBF, 2048);
			free(out_fname);
		}

		char fastq_realpath[PATH_MAX + 1];
		realpath(pd->opt->fastq_file, fastq_realpath);

		/* generate some useful information in the output */
		time(&rt);
		fprintf(pd->log_fp, "# HMMEC starts at %s\n", ctime(&rt));
		fprintf(pd->log_fp, "# Sequencing FASTQ datafile: %s\n", 
				fastq_realpath);	

		if (pd->opt->auto_kmer_length)
			fprintf(pd->log_fp, "# Kmer length: %d\n", kmer_len);
		else
			fprintf(pd->log_fp, "# Kmer length: [auto]\n");

		fprintf(pd->log_fp, "# Max substitution errors: %d\n",
				pd->opt->max_hamming_dist);

		if (pd->opt->output_params) {
			out_fname = io_cat_output_filename(pd->opt, "", "model");
			if (out_fname != NULL) {
				pd->params_fp = fopen(out_fname, "w+");
				free(out_fname);
			}
		}
	}

	INFO("Creating memory pools...\n");
	pd->mp = mempool_create(MEMPOOL_DEFAULT_SIZE);

	int k_min = pd->opt->kmer_length - 4;
	int k_max = 30; 
	int *read_selector = NULL;

	dictsize_t dict_init_size = (dictsize_t) pow(4, kmer_len < 22 ? 11 : kmer_len - (kmer_len >> 1));

	if (pd->opt->auto_kmer_length) {
		INFO_TEE(pd->log_fp, "Parameter tuning: finding the optimal k...\n");

		double prev_nonuniq_prop = 0.0;
		uint32_t prev_model_complexity = 0;

		for (int kmer_len = k_min; kmer_len <= k_max; kmer_len += 2) {
			pd->tmpmp = mempool_create(MEMPOOL_DEFAULT_SIZE);

			pd->opt->kmer_length = kmer_len;
			pd->transition_p = NULL;
			pd->exp_trans_p = NULL;
			pd->avg_initial_p = NULL;
			pd->n_total_trans = 0;

			pd->kmer_dict = dict_create(dict_init_size);

			INFO_TEE(pd->log_fp, "Constructing k-spectrum for k=%d...\n", pd->opt->kmer_length);
			io_iterate_reads_fastq(pd, pd->opt->n_threads, false, read_parse, resize_kspectrum, NULL);

			dict_init_size = pd->kmer_dict->sizemask + 1;

			if (kmer_len == k_min) {
				read_selector = calloc(pd->n_reads, sizeof(*read_selector));

				if (pd->opt->qual_score_offset <= 64) 
					pd->opt->qual_score_offset = 33;
				else
					pd->opt->qual_score_offset = 64;

				pd->opt->max_qual_score -= pd->opt->qual_score_offset;

				pd->preconstructed_nbhds = calloc(pd->n_reads, sizeof(kmer_nbhd_t));
				pd->read_em_time_cov = calloc(pd->n_reads << 1, sizeof(double));
				int _tmax = (pd->max_read_length - pd->opt->kmer_length + 1);
				pd->n_hd_windows = _tmax / pd->opt->hd_window_size;
				if (pd->n_hd_windows * pd->opt->hd_window_size < pd->max_read_length) ++pd->n_hd_windows;
				pd->read_dmax = calloc(pd->n_reads * pd->n_hd_windows, sizeof(uint8_t));
				pd->read_start_positions = calloc(pd->n_reads, sizeof(int));

				INFO_TEE(pd->log_fp, "Quality score format: Phred+%d, range : (0, %d)\n", 
						pd->opt->qual_score_offset, pd->opt->max_qual_score);

				io_iterate_reads_fastq(pd, pd->opt->n_threads, pd->opt->is_tty,
						read_qual_select, NULL, read_selector);

				int n_sel_reads = 0;
				for (int i = 0; i < pd->n_reads; ++i) 
					n_sel_reads += read_selector[i];

				// randomly sample 1% of the reads for probing model complexity
				sample_probing_reads(pd, read_selector, n_sel_reads, (int) ((double) n_sel_reads * 0.01));

			}

			INFO_TEE(pd->log_fp, "K-spectrum (k=%d) size: %u\n", kmer_len, pd->kmer_dict->used);

			INFO_TEE(pd->log_fp, "Adjusting kmer labeling (to leverage double-stranded dependence)...\n");
			adjust_kmer_labeling(pd);

			INFO_TEE(pd->log_fp, "Allocating memory for transition parameters...\n");
			allocate_trans_params(pd, 1);

			INFO_TEE(pd->log_fp, "Determining read-specific max. Hamming distances...\n");
			io_iterate_reads_fastq(pd, pd->opt->n_threads, pd->opt->is_tty, hmm_kcov_quantile, NULL, NULL);
			qsort(pd->read_em_time_cov, pd->n_reads, sizeof(double), cmp_dbl);

			pd->opt->kcov_q75 = pd->read_em_time_cov[(int) (0.75 * pd->n_reads)];
			pd->opt->kcov_q90 = pd->read_em_time_cov[(int) (0.90 * pd->n_reads)];
			INFO_TEE(pd->log_fp, "Avg. kmer coverage Q=75/90 quantile: %lg/%lg\n",
					pd->opt->kcov_q75, pd->opt->kcov_q90);
			io_iterate_reads_fastq(pd, pd->opt->n_threads, pd->opt->is_tty, hmm_determine_dmax, NULL, read_selector);

			INFO_TEE(pd->log_fp, "Building first neighborhoods for selected (probing) reads...\n");
			io_iterate_reads_fastq(pd, pd->opt->n_threads, pd->opt->is_tty, hmm_build_1st_nbhd, NULL, read_selector);

			pd->model_complexity = 0;
			pd->model_obs_complexity = 0;

			INFO_TEE(pd->log_fp, "Probing neighborhood size (model complexity)...\n");
			io_iterate_reads_fastq(pd, pd->opt->n_threads, pd->opt->is_tty, hmm_probe_nbhd_size, NULL,
					read_selector);

			double Ck = (double) prev_model_complexity / pd->model_complexity;
			double Rk = (double) prev_nonuniq_prop / pd->prop_nonuniq_kmers;

			INFO_TEE(pd->log_fp, "Model complexity (k=%d): %u/%u, ratio=%lg. \n"
					"Prop. of non-unique (solid) kmers: %lg/%lg, ratio=%lg\n",
					kmer_len, pd->model_complexity, prev_model_complexity, Ck,
					pd->prop_nonuniq_kmers, prev_nonuniq_prop, Rk);

			prev_model_complexity = pd->model_complexity;
			prev_nonuniq_prop = pd->prop_nonuniq_kmers;

			// data cleaning
			free(pd->transition_p);
			free(pd->avg_initial_p);
			free(pd->exp_trans_p);

			dict_destroy(pd->kmer_dict);
			mempool_destroy(pd->tmpmp);

			memset(pd->preconstructed_nbhds, 0, pd->n_reads * sizeof(kmer_nbhd_t));
			memset(pd->read_em_time_cov, 0, pd->n_reads * 2 *  sizeof(double));
			memset(pd->read_dmax, 0, pd->n_reads * pd->n_hd_windows * sizeof(uint8_t));

			if (Ck > 0.0 && Ck <= 1.25 && Rk > 0.0 && Rk <= 1.25) {
				pd->opt->kmer_length = kmer_len - 2;
				INFO_TEE(pd->log_fp, "Kmer length selected to be: %d\n", pd->opt->kmer_length);
				break;
			}

			if (kmer_len == k_max) {
				pd->opt->kmer_length = kmer_len;
				INFO_TEE(pd->log_fp, "Kmer length selected to be: %d\n", pd->opt->kmer_length);
				break;
			}
		}
	}


	pd->kmer_dict = dict_create(dict_init_size);

	INFO_TEE(pd->log_fp, "Constructing k-spectrum for k=%d...\n", pd->opt->kmer_length);
	io_iterate_reads_fastq(pd, pd->opt->n_threads, false, read_parse, resize_kspectrum, NULL);

	/* Iterate over all available reads provided by the FASTQ data file,
	 * with callback function 'read_parse'.
	 *
	 * After this step, we will fill the (*d)->kmer_dict buckets with
	 * observed kmers (k-spectrum).
	 */

	INFO_TEE(pd->log_fp, "Size of k-spectrum: %zu\n", (size_t) pd->kmer_dict->used);

	if (pd->opt->qual_score_offset <= 35) 
		pd->opt->qual_score_offset = 33;
	else if (pd->opt->qual_score_offset > 35 && pd->opt->qual_score_offset <= 66) 
		pd->opt->qual_score_offset = 64;

	pd->opt->max_qual_score -= pd->opt->qual_score_offset;

	INFO_TEE(pd->log_fp, "Quality score format: Phred+%d, range : (0, %d)\n", 
			pd->opt->qual_score_offset, pd->opt->max_qual_score);

	int qmax = pd->opt->max_qual_score + 1;
	int tmax = pd->max_read_length - pd->opt->kmer_length + 1;

	pd->n_emit_windows = (tmax-1) / pd->opt->emit_window_size;
	if (pd->n_emit_windows * pd->opt->emit_window_size < (tmax-1)) ++pd->n_emit_windows;

#if PMR_NO_DEP_EMISSION
	pd->n_base_context = 4; // [ACGT] x {0, 1, 2, 4, 8, 16}
	pd->base_ctx_radix = 4 * 5; // 5 possible emission states x 4 contexts
#else
	pd->n_base_context = 4 * 6; // [ACGT] x {0, 1, 2, 4, 8, 16}
	pd->base_ctx_radix = 4 * 6 * 5; // 5 possible emission states x (4*6) contexts
#endif
	pd->n_qual_context = 2; // {0, 1} 
	pd->qual_ctx_radix = 2* qmax; // qmax possible emission states x (2) contexts

#if PMR_USE_PRIMITIVE_EMISSION
	pd->base_error_rate = pd->opt->error_rate;
#else
	INFO_TEE(pd->log_fp, "Allocating memory for HMM parameters...\n");
	/* No context is used for emitting the first kmer, only base by base
	 * emission assuming independence. Therefore 4 (ACTG) * 5 (ACTGN) 
	 * possibilities */
	size_t alloc_base_emit_slots = (4 * 5 * pd->opt->kmer_length +
			pd->n_emit_windows * pd->base_ctx_radix) * 
		sizeof(*pd->base_emission_p);

	/* similarly, within the first kmer, each position can emit 2 x qmax
	 * possible y_i[t] distinct quality scores.
	 */
	size_t alloc_qual_emit_slots = ((qmax << 1) * pd->opt->kmer_length +
			pd->n_emit_windows * pd->qual_ctx_radix) * 
		sizeof(*pd->qual_emission_p);

	pd->base_emission_p = mempool_alloc(pd->mp, alloc_base_emit_slots);
	pd->qual_emission_p = mempool_alloc(pd->mp, alloc_qual_emit_slots);

	pd->base_emit_exp = mempool_alloc(pd->mp, 
			pd->opt->n_threads * alloc_base_emit_slots);
	pd->qual_emit_exp = mempool_alloc(pd->mp, 
			pd->opt->n_threads * alloc_qual_emit_slots);
#endif

//	if (pd->opt->mode == PMR_MODE_ERROR_CORRECTION) {
		INFO_TEE(pd->log_fp, "Adjusting kmer labeling (to leverage double-stranded dependence)...\n");
		adjust_kmer_labeling(pd);
//	}

	INFO_TEE(pd->log_fp, "Allocating memory for transition parameters...\n");
	allocate_trans_params(pd, 1);

	//INFO("Probing neighborhood size (model complexity)...\n");
	//probe_model_complexity(pd);

	if (pd->preconstructed_nbhds == NULL)
		pd->preconstructed_nbhds = calloc(pd->n_reads, sizeof(kmer_nbhd_t));

	if (pd->read_start_positions == NULL)
		pd->read_start_positions = calloc(pd->n_reads, sizeof(int));

	if (pd->read_em_time_cov == NULL)
		pd->read_em_time_cov = calloc(pd->n_reads << 1, sizeof(double));

	if (pd->read_dmax == NULL) {
		if (pd->n_hd_windows == 0) {
			pd->n_hd_windows = tmax / pd->opt->hd_window_size;
			if (pd->n_hd_windows * pd->opt->hd_window_size < tmax) ++pd->n_hd_windows;
		}
		pd->read_dmax = calloc(pd->n_reads * pd->n_hd_windows, sizeof(uint8_t));
	}
	
	

	/*
	INFO_TEE(pd->log_fp, "Determining read-specific max. Hamming distances...\n");
	io_iterate_reads_fastq(pd, pd->opt->n_threads, pd->opt->is_tty, hmm_kcov_quantile, NULL);
	qsort(pd->read_em_time_cov, pd->n_reads, sizeof(double), cmp_dbl);

	pd->opt->kcov_q75 = pd->read_em_time_cov[(int) (0.75 * pd->n_reads)];
	pd->opt->kcov_q90 = pd->read_em_time_cov[(int) (0.90 * pd->n_reads)];
	INFO_TEE(pd->log_fp, "Avg. kmer coverage Q=75/90 quantile: %lg/%lg\n",
			pd->opt->kcov_q75, pd->opt->kcov_q90);
	io_iterate_reads_fastq(pd, pd->opt->n_threads, pd->opt->is_tty, hmm_determine_dmax, NULL);
	*/

	if (pd->opt->mode == PMR_MODE_ERROR_CORRECTION) {
		readjust_kmer_labeling(pd);

		//io_iterate_reads_fastq(pd, hmm_dump_strandness, pd->sdump_fp);
		INFO_TEE(pd->log_fp, "Building first neighborhoods for all reads...\n");
		io_iterate_reads_fastq(pd, pd->opt->n_threads, pd->opt->is_tty, hmm_build_1st_nbhd, NULL, NULL);
		
		INFO_TEE(pd->log_fp, "Determining read-specific max. Hamming distances...\n");
		io_iterate_reads_fastq(pd, pd->opt->n_threads, pd->opt->is_tty, hmm_adaptive_dmax, NULL, NULL);
	}

	return NO_ERROR;
}

int destroy_data(data *d)
{
	/* destroy all dictionaries */
	dict_destroy(d->kmer_dict);

	/* release all memory allocated for memory pools */
	mempool_destroy(d->mp);

	free(d->preconstructed_nbhds);
	free(d->read_start_positions);
	free(d->opt);
	free(d);

	return 0;
}

int make_options(options **opt, int argc, char **argv)
{
	static struct option pmr_options[] = {
		{"eta", required_argument, 0, 'E'},
		{"gamma", required_argument,  0, 'G'},
		{"rel-tol", required_argument, 0, 'R'},
		{"ground-truth-file", required_argument, 0, 'g'},
		{"quality-score-pmf", required_argument, 0, 'Q'},
		{"penalized-em", required_argument, 0, 0},
		{"penalized-init", required_argument, 0, 0},
		{"dump-strandness", no_argument, 0, 0},
		{"error-rate", required_argument, 0, 0},
		{"diagnose", optional_argument, 0, 0},
		{"qscore-bin-width", required_argument, 0, 0},
		{"error-watch", required_argument, 0, 0},
		{"partial-decoding", required_argument, 0, 0},
		// --phred [1|2]
		// 1: treat quality scores as true sequencing error rates, except for
		// quality score 2, for which the error rate is provided by -E option.
		// 2: treat ALL quality scores as true sequencing error rates,
		// INCLUDING quality score 2, which typically indicates low-quality
		// regions and do not encode error rates.
		{"phred", required_argument, 0, 0},
		{0, 0, 0, 0},
	};

	size_t arg_len;
	int invalid_opts = 0;

	*opt = (options *) malloc(sizeof(**opt));

	options *popt = *opt;

	popt->n_threads = 1;
	popt->kmer_length = 0;
	popt->max_hamming_dist = 0;
	popt->max_iter = -1;
	popt->preconstructed_nbhd_hd = 1;
	popt->pagesize = getpagesize();
	popt->qual_score_offset = 255;
	popt->max_qual_score = 40;
	popt->qscore_bin_width = 1;
	popt->random_param_init = 0;
	popt->filename_append_timeinfo = 1;
	popt->filename_append_params = 1;
	popt->dump_strandness = 0;
	popt->diagnose_mode = 0;
	popt->update_qual_score = 1;
	popt->enable_error_watch = 0;
	popt->qscore_as_phred = 0;
	popt->viterbi_revcomp = 0;
	popt->output_ordered = 1;
	popt->adaptive_hamming_dist = 1;
	popt->partial_decoding_lvl = 0;
	popt->auto_penalty_eta = 0;
	popt->auto_kmer_length = 0;
	popt->qual_emit_context = 1;
	popt->base_emit_context = 3;
	popt->em_penalty_func = 2;
	popt->init_penalty_func = 2;
	popt->max_approx_nbhd_size = 1000;

	popt->is_tty = isatty(fileno(stderr));
	popt->mstep_debug_console = false;

	// backtrack
	popt->n_max_kmer_contexts = 127;
	popt->max_backtrack_dist = 0;

	popt->ground_truth_file = NULL;
	popt->output_prefix = NULL;
	/* default emission distribution */
	popt->output_params = 0;

	popt->disable_penalized_em = 0;
	popt->penalize_init_params = 1;
	popt->error_rate = 0.01;

	popt->penalty_gamma = 1e-20;
	popt->penalty_eta = NAN;
	popt->tol = 1e-4;

	popt->mode = 0;

	if (argc == 1) {
		print_usage(argv[0]);
		exit(1);
	}

	/* command */
	const char *command = argv[1];
	if (strncmp(command, "oracle", 6) == 0) {
		popt->mode = PMR_MODE_ORACLE;
	}
	else if (strncmp(command, "poracle", 7) == 0) {
		popt->mode = PMR_MODE_PERFECT_ORACLE;
	}
	else if (strncmp(command, "correct", 7) == 0) {
		popt->mode = PMR_MODE_ERROR_CORRECTION;
	}

	if (popt->mode == 0) {
		message(stderr, __FILE__, __LINE__,
				ERROR_MSG, INVALID_CMD_ARGUMENT, 
				"invalid command: \"%s\"\n", argv[1]);

		print_usage(argv[0]);
		exit(1);
	}

	/* parse user specified options */
	char *ptr_optarg;
	char c;
	int opt_index;
	while ((c = getopt_long(argc - 1, &argv[1], "k:d:o:g:E:G:R:AP:B:Q:I:t:pVc", pmr_options,
					&opt_index)) != -1) {
		switch (c) {
			case 0:	/* currently nonfunctional */
				if (strncmp(pmr_options[opt_index].name, "diagnose", 
							8) == 0) {

					if (optarg && *optarg == 'f')
						popt->diagnose_mode = 1;
					else if (optarg && *optarg == 'l')
						popt->diagnose_mode = 0;
					else if (optarg == 0)
						popt->diagnose_mode = 0;
					else {
						/* INVALID FORMAT */
					}

					INFO("[OPTION] Diagnose mode <%s> enabled by user.\n",
							popt->diagnose_mode ? "full" : "lite");
				}
				else if (strncmp(pmr_options[opt_index].name, 
							"penalized-em", 12) == 0) {
					popt->em_penalty_func = strtol(optarg, &ptr_optarg, 0);
					if (popt->em_penalty_func < 0 || popt->em_penalty_func > 2) {
						/* error */
					}
				}
				else if (strncmp(pmr_options[opt_index].name,
							"partial-decoding", 16) == 0) {
					popt->partial_decoding_lvl = strtol(optarg, &ptr_optarg, 0);
					if (popt->partial_decoding_lvl > 2) 
						popt->partial_decoding_lvl = 2;
				}
				else if (strncmp(pmr_options[opt_index].name, 
							"penalized-init", 14) == 0) {
					popt->init_penalty_func = strtol(optarg, &ptr_optarg, 0);
					if (popt->init_penalty_func < 0 || popt->init_penalty_func > 2) {
						/* error */
					}
				}
				else if (strncmp(pmr_options[opt_index].name,
							"qscore-bin-width", 16) == 0) {
					popt->qscore_bin_width = strtol(optarg, &ptr_optarg, 0);
				}
				else if (strncmp(pmr_options[opt_index].name,
							"error-rate", 10) == 0) {
					popt->error_rate = strtod(optarg, &ptr_optarg);
					//init_par_watch_map(optarg, popt->kmer_length);
				}
				else if (strncmp(pmr_options[opt_index].name,
							"error-watch", 11) == 0) {
					popt->enable_error_watch = 1;
					//init_par_watch_map(optarg, popt->kmer_length);
				}
				else if (strncmp(pmr_options[opt_index].name,
							"phred", 5) == 0) {
					// Initialize quality scores as if they accurately measure
					// sequencing error rate, according to PHRED quality score
					// formula.
					//
					// Use in Oracle mode.
					int phred_mode = strtol(optarg, &ptr_optarg, 0);
					if (phred_mode > 2) phred_mode = 2;
					if (phred_mode < 0) phred_mode = 0;
					popt->qscore_as_phred = phred_mode; 
				}
				break;
				/* kmer-length */
			case 'k':
				popt->kmer_length = strtol(optarg, &ptr_optarg, 0);
				if (*ptr_optarg != '\0') {
					if (strncmp(ptr_optarg, ",auto", 5) == 0) {
						popt->auto_kmer_length = 1;
					}
				}

				if (ptr_optarg == optarg){
					message(stderr, __FILE__, __LINE__,
							ERROR_MSG, INVALID_CMD_ARGUMENT, 
							"-k %s\n", optarg);

					exit(INVALID_CMD_ARGUMENT);
				}

				/*
				   if (popt->kmer_length <= 0 || popt->kmer_length > 31){
				   exit(INVALID_CMD_ARGUMENT);
				   }
				   */

				break;
				/* maximum substitutional errors */
			case 'd':
				popt->max_hamming_dist = strtol(optarg, &ptr_optarg, 0);

				if (ptr_optarg == optarg){
					message(stderr, __FILE__, __LINE__,
							ERROR_MSG, INVALID_CMD_ARGUMENT, 
							"-d %s\n", optarg);

					exit(INVALID_CMD_ARGUMENT);
				}

				if (popt->max_hamming_dist < 1 || popt->max_hamming_dist > 15 ||
						popt->max_hamming_dist > popt->kmer_length) {
					message(stderr, __FILE__, __LINE__, ERROR_MSG, 
							INVALID_CMD_ARGUMENT, 
							"-d %d : should be within range [1, min(k, 15)].\n", 
							popt->max_hamming_dist);

					invalid_opts |= 1;
				}

				break;
			case 't':
				popt->n_threads = strtol(optarg, &ptr_optarg, 0);

				if (ptr_optarg == optarg){
					message(stderr, __FILE__, __LINE__,
							ERROR_MSG, INVALID_CMD_ARGUMENT, 
							"-t %s\n", optarg);

					exit(INVALID_CMD_ARGUMENT);
				}

				break;
			case 'I':
				popt->max_iter = strtol(optarg, &ptr_optarg, 0);
				if (ptr_optarg == optarg){
					message(stderr, __FILE__, __LINE__,
							ERROR_MSG, INVALID_CMD_ARGUMENT, 
							"-I %s\n", optarg);

					exit(INVALID_CMD_ARGUMENT);
				}

				break;

			case 'c':
				popt->mstep_debug_console = 1;
				break;

			case 'E':
				popt->penalty_eta = strtod(optarg, &ptr_optarg);

				if (!isnan(popt->penalty_eta) && popt->penalty_eta < 0) {
					message(stderr, __FILE__, __LINE__,
							ERROR_MSG, INVALID_CMD_ARGUMENT, 
							"-L %s : negative lambda value.\n", optarg);
					invalid_opts |= 1;
				}

				break;
			case 'G':
				popt->penalty_gamma = strtod(optarg, &ptr_optarg);

				if (popt->penalty_gamma < 0) {
					message(stderr, __FILE__, __LINE__,
							ERROR_MSG, INVALID_CMD_ARGUMENT, 
							"-G %s : negative gamma value.\n", optarg);
					invalid_opts |= 1;
				}

				break;

			case 'Q':
				arg_len = strlen(optarg);
				popt->qscore_pmf_file = calloc(arg_len + 1, sizeof(char));
				strncpy(popt->qscore_pmf_file, optarg, arg_len);
				break;

			case 'g':
				arg_len = strlen(optarg);
				popt->ground_truth_file = calloc(arg_len + 1, sizeof(char));
				strncpy(popt->ground_truth_file, optarg, arg_len);
				break;

			case 'R':
				popt->tol = strtod(optarg, &ptr_optarg);
				break;

			case 'A':
				popt->adaptive_hamming_dist = 0;
				//popt->max_approx_nbhd_size = strtol(optarg, &ptr_optarg, 0);
				break;

				/* output all trained parameters to local file PREFIX.param */
			case 'p':
				popt->output_params = 1;
				break;
			case 'V':
				popt->viterbi_revcomp = 1;
				break;
				/* load parameters from file */
			case 'P':
				arg_len = strlen(optarg);
				popt->params_data = calloc(arg_len + 1, sizeof(char));
				strncpy(popt->params_data, optarg, arg_len);
				break;


				/* prefix for log files */
			case 'o':
				arg_len = strlen(optarg);
				popt->output_prefix = calloc(arg_len + 1, sizeof(char));
				strncpy(popt->output_prefix, optarg, arg_len);

				break;
			case 'C':

				break;
		}
	}

	if (optind < argc - 1) {
		popt->fastq_file = argv[optind + 1];
	}
	else {
		invalid_opts |= 1;
	}

	if (isnan(popt->penalty_eta)) {
		popt->auto_penalty_eta = 1;
	}

	if ((popt->mode & 1) && popt->ground_truth_file == NULL) {
		message(stderr, __FILE__, __LINE__, ERROR_MSG, 
				INVALID_CMD_ARGUMENT, 
				"you must specify --ground-truth-file under 'oracle' mode.\n");
		invalid_opts |= 1;
	}

	/*
	if (popt->mode == PMR_MODE_ORACLE && popt->qscore_pmf_file == NULL) {
		message(stderr, __FILE__, __LINE__, ERROR_MSG, 
				INVALID_CMD_ARGUMENT, 
				"you must specify --qscore-pmf-file|-Q under 'oracle' mode.\n");
		invalid_opts |= 1;
	}
	*/

	/* validate specified options/arguments */
	if (popt->kmer_length < 4 || popt->kmer_length > 31) {
		message(stderr, __FILE__, __LINE__, ERROR_MSG, 
				INVALID_CMD_ARGUMENT, 
				"kmer size should be within range [4, 31].\n");
		invalid_opts |= 1;
	}

	if (popt->tol < 0 || popt->tol > 1) {
		message(stderr, __FILE__, __LINE__,
				ERROR_MSG, INVALID_CMD_ARGUMENT, 
				"-R %s : should be within range (0, 1).\n", optarg);
		invalid_opts |= 1;
	}

	popt->max_backtrack_dist = popt->kmer_length; 
	// popt->max_backtrack_dist = popt->kmer_length + (popt->kmer_length >> 1);

	if (invalid_opts) {
		print_usage(argv[0]);
		exit(1);
	}

	popt->max_qual_score /= popt->qscore_bin_width;

	if (popt->max_hamming_dist == 0)
		popt->max_hamming_dist = popt->kmer_length >> 1;

	popt->hd_window_size = popt->kmer_length >> 2;
	popt->emit_window_size = popt->kmer_length >> 2;

	return NO_ERROR;
}

void print_usage(const char *basename)
{
	fprintf(stderr, pmr_usage, basename);
}

void timer(int report)
{
	static struct timeval tv;
	if (report == 0)
		gettimeofday(&tv, NULL);
	else {
		struct timeval _ntv;
		gettimeofday(&_ntv, NULL);
		fprintf(stderr, 
				"[TIMER] Elapsed: %.3f\n",
				(0.0 + _ntv.tv_sec + (double) _ntv.tv_usec/1e6) - 
				(0.0 + tv.tv_sec + (double) tv.tv_usec/1e6));

		tv = _ntv;
	}
}

/* select n reads from a population of N, without replacement */
static void sample_probing_reads(data *d, int *rsel, int N, int n)
{
	int t = 0, m = 0;
	int tstar = 0, p = 0;

	while (m < n) {
		double u = (double) rand() / RAND_MAX;
		if ((N - t) * u >= (n - m)) {
			++t;
		}
		else {
			for (; tstar <= t; ++p) {
				if (rsel[p]) {
					if (tstar < t) rsel[p] = 0;

					++tstar;
				}
			}

			// find the t-th read with selector "1"
			++t;
			++m;
		}
	}

	for (++p; p < d->n_reads; ++p) {
		rsel[p] = 0;
	}
}
