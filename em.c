#include <omp.h>
#include <stdbool.h>
#include <signal.h>

#include "em.h"
#include "read.h"
#include "mstep.h"
#include "kmer.h"
#include "iodata.h"
#include "trans.h"
#include "numeric.h"
#include "quickselect.h"
#include "sort.h"
#include "kvec.h"
//#include "pardump.h"
#include "atomic.h"

//static double avg_nbhd_size = 0.0;
static int _slowexp = 0;
static int _zero_likelihood = 0;
static int _iter = 0;
//static int _max_nbhd_size = 0;

static void hmm_dump_model(data *d);
static double hmm_normalize_emission(double *emit_p, kbits_t enf, 
		kbits_t xt);
static void hmm_update_expected_vals(data *d, state_nbhd_t *curr_nbhd, 
		int t, double *t_exp_counts, double *t_exp_trans, 
		double *t_base_emit_exp, double *t_qual_emit_exp,
		int kmer_len, int ctx_shift, int qmax, 
		const char *rseq, const char *conv_q,
		double inv_sum_exp_count);
static void hmm_update_qual_scores(data *d, state_nbhd_t *curr_nbhd,
		int t, read_t *read, read_t *argmax_read, double *t_exp_counts,
		double inv_sum_exp_count);

static void hmm_viterbi_wrapper(data *d, read_t *read, void *fdata);
static void hmm_write_to_fastq(data *d, void *fdata);

//static int kmer_suffix_comp(const void *a, const void *b);
/*
static int *compute_qual_contexts(const char *conv_q, int rlen, 
		int klen, int ctxlen, int qcnt, mempool_t *mp);
*/
static void hmm_m_step(data *d, int mstep_flag);

void hmm_e_step(data *d, read_t *read, void *fdata, read_t *argmax_read);
void hmm_viterbi(data *d, read_t *read, void *fdata);
static void dump_nbhd(state_nbhd_t *nbhd, FILE *fp);

static void dump_nbhd(state_nbhd_t *nbhd, FILE *fp)
{
	for (int i = 0; i < nbhd->size; ++i) {
		fprintf(fp, "%lu/%u,",
				nbhd->states_sorted_sfx[i],
				nbhd->kmer_ptrs[i]->trans_flag);
	}
}

void EM(data *d, struct timeval *ts)
{
	double pll = 0.0, delta = 0.0;
	int iter = 1;
	_iter = iter;

	/* initialize HMM parameters */
	hmm_init_model_params(d);

	if (d->opt->max_iter != 0) {
#ifdef PMR_APPROX_HMM
		INFO("Running iterative approximate Baum-Welch algorithms... "
				" (Maximal neigborhood size: %zu)\n", 
				d->opt->max_approx_nbhd_size);
#else
		INFO_TEE(d->log_fp, "Running iterative Baum-Welch algorithms...\n");
#endif
	}
	else {
		INFO("Skipped HMM fitting (using Baum-Welch algorithm).\n");
	}


	do {
		if (d->opt->max_iter >= 0 && iter >= d->opt->max_iter) break;

		d->loglik = 0.0;

		_slowexp = 0;
		_zero_likelihood = 0;

		INFO_TEE(d->log_fp, "E-step / iteration: %2d\n", iter);

		hmm_init_e_step(d);
		io_iterate_reads_fastq(d, d->opt->n_threads, d->opt->is_tty, hmm_e_step_wrapper, NULL, NULL);
		//INFO_TEE(d->log_fp, "Slow __slowexp invocations: %d\n", _slowexp);
		INFO_TEE(d->log_fp, "Reads with zero likelihoods: %d\n", _zero_likelihood);

/*
		if (d->opt->mstep_debug_console) {
			debug_console(d);
		}
*/

		INFO_TEE(d->log_fp, "M-step / iteration: %2d\n", iter);
		hmm_m_step(d, 0);

		delta = d->loglik - pll;

		if (_iter == 1) {
#if 0
			char *_dump_fname = io_cat_output_filename(d->opt, "", "timecov");	
			FILE *ftime_dump = fopen(_dump_fname, "wb+");
			fwrite(d->read_em_time_cov, sizeof(double), d->n_reads << 1, ftime_dump);

			free(_dump_fname);
			fclose(ftime_dump);
#endif
		}

		INFO_TEE(d->log_fp, "[EM/ITER %2d] LL: " DOUBLE_E_FMT " PLL: " DOUBLE_E_FMT "("
				DOUBLE_E_FMT ")\n", 
				iter, d->loglik, pll, delta);
		++iter;
		_iter = iter;

		pll = d->loglik;
	} while (FABS(delta / d->loglik) > d->opt->tol);

	if (d->opt->output_params)
		hmm_dump_model(d);

	/*
	INFO_TEE(d->log_fp, "Fitting a new model with labels swapped [E-step].\n");

	hmm_init_e_step(d);
	io_iterate_reads_fastq(d, d->opt->n_threads, d->opt->is_tty, hmm_e_step_wrapper, NULL);
	*/

	/*
	INFO_TEE(d->log_fp, "M-step with swapped labels.\n");
	// hmm_m_step(d, MSTEP_FLIP_LABEL);
	hmm_m_step(d, 0);
	*/

	/*
	INFO_TEE(d->log_fp, "Removing strand-specific erroneous transitions.\n");

	hmm_init_e_step(d);
	io_iterate_reads_fastq(d, d->opt->n_threads, d->opt->is_tty, hmm_e_step_wrapper, NULL);

	mstep_strand_specific_errors(d);
	*/

	//hmm_m_step(d, MSTEP_FLIP_LABEL | MSTEP_NO_EMISSION);

	/*
	if (d->opt->mstep_debug_console) {
		hmm_m_step_debug_console(d);
	}
	*/

	//debug_console(d);

#if !PMR_HMMSBG
	INFO_TEE(d->log_fp, "Starting error-correction using Viterbi algorithm...\n");
	char **fastq_output = NULL;
	if (d->opt->output_ordered) fastq_output = (char **) 
		calloc(PMR_BATCH_GROUP_SIZE, sizeof(*fastq_output));
	io_iterate_reads_fastq(d, d->opt->n_threads, d->opt->is_tty, hmm_viterbi_wrapper, 
			d->opt->output_ordered ? hmm_write_to_fastq : NULL, fastq_output);
	free(fastq_output);
#endif

#if PMR_HMMSBG
	if (d->opt->max_iter != 0) {
		INFO_TEE(d->log_fp, "Readjust kmer labeling...\n");
		readjust_kmer_labeling(d);
	}
#endif

}

/***
 * Maximization step in EM algorithm to update model parameters.
 */
void hmm_m_step(data *d, int mstep_flag)
{
	/* Estimate new set of parameters by maximizing the (penalized) 
	 * likelihood */
	if (d->opt->mode == PMR_MODE_ERROR_CORRECTION) {
		cmstep_penalized_transition(d, d->opt->em_penalty_func, 
				mstep_flag & MSTEP_FLIP_LABEL);
	}

	if (!(mstep_flag & MSTEP_NO_EMISSION)) {
#if PMR_ALTERNATIVE_EMISSION	
		// skip updating emission parameters if --phred is specified
		if (d->opt->qscore_as_phred == 0 || d->opt->mode == PMR_MODE_ERROR_CORRECTION) {
			mstep_alt_emission(d);
		}
#else
		mstep_emission(d);
#endif
	}

	//mstep_extend_all_kmers(d);
	//mstep_initial(d);
}

void hmm_probe_nbhd_size(data *d, read_t *read, void *fdata)
{
	int actg_trans_count[4];
	int cum_actg_trans_count[4];
	int total_hmm_trans = 0;

	const int kmer_len = d->opt->kmer_length;
	const int num_kmers = read->length - kmer_len + 1;
	const int shift = (kmer_len - 1) << 1;
	const int bwidth = d->opt->qscore_bin_width;
	const double _uniq_tp = 0.0;

	int tmin = 0;
	int tmax = read->length - kmer_len;
	if (tmin >= tmax) return; 

	int *rsel = (int *) fdata;
	if (rsel != NULL && !rsel[read->id]) return;

	mempool_t *mp = mempool_create(MEMPOOL_DEFAULT_SIZE);

	dict *kdict = d->kmer_dict;

	/* sequence */
	char *rseq = read->sequence;
	char *conv_q = convert_quality_scores(read->qscore, read->length, 
			d->opt->qual_score_offset, bwidth, mp);

	state_nbhd_t *T_nbhds = mempool_nalloc(mp,
			sizeof(*T_nbhds) * num_kmers, 16);

	/* flag for N (not determined) base, first kmer only */
	kbits_t obs_n_flag = kmer_n_base_flag(rseq, kmer_len);
	kbits_t next_obs_n_flag = kmer_n_base_flag(rseq + 1, kmer_len);

	/* convert string literal to numeric representation */
	kbits_t observed_kmers[BITS_TO_U64(read->length << 1)];
	read_seq_to_numeric(rseq, observed_kmers, read->length, kmer_len);

	/* set up first kmer neighborhood */
	kbits_t obs_kid = bitarray_get(observed_kmers, tmin << 1, kmer_len << 1);
	//ksubstr_t obs_base = kmer_effective_base(rseq[tmin+kmer_len-1]);
	ksubstr_t obs_next_base = kmer_effective_base(rseq[tmin+kmer_len]);

	// set the first neighborhood to be the observed kmer itself
	T_nbhds[tmin].size = hmm_load_preconstructed_nbhd(
			d, read->id, obs_kid, 
			obs_n_flag, conv_q, &T_nbhds[tmin], kmer_len, 0, kdict, mp, 
			&_uniq_tp);

	ksubstr_t dmax = d->read_dmax[read->id];
	int orig_dmax = (int) dmax;

	/* initial state distribution */
	for (int t = tmin; t < tmax; t++) {

		obs_next_base = kmer_effective_base(rseq[t + kmer_len]);
		obs_kid = bitarray_get(observed_kmers, t << 1, kmer_len << 1);

		kbits_t _new_nf = ((obs_next_base >> 2) << 1) | 
			(obs_next_base >> 2);
		next_obs_n_flag = (next_obs_n_flag >> 2) | (_new_nf << shift);

		state_nbhd_t *curr_nbhd = &T_nbhds[t],
					 *next_nbhd = &T_nbhds[t+1];

		//const ksubstr_t dmax = d->opt->max_hamming_dist;

		//int t_nbhd_size = curr_nbhd->size;

		//total_nbhd_size += t_nbhd_size;

		kbits_t *t_states_sorted_sfx = T_nbhds[t].states_sorted_sfx;

		/* ----- II. compute one iteration of forward algorithm ----- */

		/* --- II a. set up kmer neighborhood at next position --- */
		memset(actg_trans_count, 0, sizeof(int) << 2);

		int n_distinct_suffixes = 0;
		int tnext_nbhd_size = 0;
		int try_incr_dmax = 1;
		while (true) {
			/* one linear scan to determine # of distinct suffix, 
			 * compute downstream actg nucleotide counts etc. */
			n_distinct_suffixes = hmm_count_distinct_suffixes(curr_nbhd, 
					actg_trans_count, dmax, obs_kid, obs_n_flag,
					obs_next_base, d->transition_p);

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

			if (tnext_nbhd_size == 0 && try_incr_dmax &&
					curr_nbhd->max_hd >= dmax && dmax <= 12) {
				++dmax;
				try_incr_dmax = 0;
				
				/*
				fprintf(stderr, "[%d/%2d] <NS: %d HD: %d> DMAX -> %d\n", read->id, t, 
						curr_nbhd->size, max_hd, dmax);
				getchar();
				*/
			}
			else {
				break;
			}
		} 

		if (tnext_nbhd_size == 0) {
			goto destroy;
		}

		/* allocate memory given tnext_nbhd_size */
		char *pmem = mempool_nalloc(mp,
				nbhd_alloc_size(tnext_nbhd_size), 16);

		next_nbhd->size = tnext_nbhd_size;
		hmm_setup_nbhd_ptrs(next_nbhd, tnext_nbhd_size, pmem);

		int kmer_idx = 0, index_pfx = 0;
		kbits_t max_hd = 0;
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
				next_nbhd->kmer_ptrs[index_sfx] = tnext_kmer;

				total_hmm_trans += ups_trans_packed & 7;

				/* set up transition flag */
				bitarray_set_pwr2(next_nbhd->ba_kmer_trans_flag,
						index_sfx, tnext_kmer->trans_flag, 2);
				/* set up Hamming distance */

				kbits_t _hd = sfx_hd + ((mismatch_flag >> dns_base) & 1);
				bitarray_set_pwr2(next_nbhd->ba_hamming_dist, 
						index_sfx, _hd, 3);

				if (_hd > max_hd) 
					max_hd = _hd;

				next_nbhd->suffixes_order[index_pfx] = index_sfx;
				next_nbhd->states_sorted_sfx[index_sfx] = tnext_kid;
				index_pfx++;
			}

			kmer_idx += n_common_sfx_kmers;
		}
		next_nbhd->max_hd = max_hd;

		//obs_base = obs_next_base;
		obs_n_flag = next_obs_n_flag;
	}

	atomic_add_u32(&d->model_complexity, total_hmm_trans);
	atomic_add_u32(&d->model_obs_complexity, tmax + 1);

	if (dmax > orig_dmax) {
		/*
		fprintf(stderr, "[%d] %s DMAX %d -> %d\n",
				read->id, read->identifier, orig_dmax, (int) dmax);
		*/
		d->read_dmax[read->id] = dmax << 1;
	}
	else {

	}

destroy:
	mempool_destroy(mp);
}

void hmm_e_step_wrapper(data *d, read_t *read, void *fdata)
{
	int rerun;
	do {
		rerun = false;
		// rerun E-step on the same read if the rerun flag is set.
		hmm_e_step(d, read, (void *) &rerun, NULL);
	} while(rerun);
}


void hmm_kcov_quantile(data *d, read_t *read, void *fdata)
{
	const int kmer_len = d->opt->kmer_length;
	int tmax = read->length - kmer_len + 1;
	dict *kdict = d->kmer_dict;

	/* convert string literal to numeric representation */
	kbits_t observed_kmers[BITS_TO_U64(read->length << 1)];
	read_seq_to_numeric(read->sequence, observed_kmers, read->length, kmer_len);

	double avg_k_cexp = 0.0;
	double min_k_cexp = 1e6;

	for (int t = 0; t < tmax; t++) {
		kbits_t obs_kid = bitarray_get(observed_kmers, t << 1, kmer_len << 1);
		kmer_t *obs_pk = (kmer_t *) &(dict_find(kdict,
					kbits_cast_to_ptr(obs_kid))->value);

		double obs_k_cexp = EXP(d->avg_initial_p[kmer_initial_index(obs_pk,
			d->uniq_kmer_init_off)]) * d->total_kmer_exp_counts;
		avg_k_cexp += obs_k_cexp;

		if (obs_k_cexp < min_k_cexp) min_k_cexp = obs_k_cexp;
	}

	avg_k_cexp /= tmax;
	d->read_em_time_cov[read->id] = avg_k_cexp;
}


/***
 * Some general algorithm in forward-backward algorithm.
 *
 */
void hmm_e_step(data *d, read_t *read, void *fdata, read_t *argmax_read)
{
	int dump_trellis = false;

	int actg_trans_count[4];
	int cum_actg_trans_count[4];
	size_t total_nbhd_size = 0;

	double ts = 0.0, te = 0.0;
	double sum_kmer_cexp_cov = 0.0;
	int sum_nbhd_size = 0;
	const int mode = d->opt->mode;

	int update_qual_score = (argmax_read != NULL);

	if (_iter == 1) ts = omp_get_wtime();

	const int kmer_len = d->opt->kmer_length;
	const int num_kmers = read->length - kmer_len + 1;
	const int shift = (kmer_len - 1) << 1;
	const int ctx_shift = (kmer_len - d->opt->base_emit_context) << 1;
	const int qmax = d->opt->max_qual_score + 1;
	const int qmax2 = qmax << 1;
	const int bwidth = d->opt->qscore_bin_width;
	const double _uniq_tp = 0.0;
#ifdef PMR_APPROX_HMM
	const size_t max_approx_nbhd_size = d->opt->max_approx_nbhd_size;
#endif

	//int iter = *(int *) fdata;
	int *rerun = fdata;

	int tmin = d->read_start_positions[read->id];
	int tmax = read->length - kmer_len;
	if (tmin >= tmax) {
		++_zero_likelihood;
		return; 
	}

	mempool_t *mp = mempool_create(MEMPOOL_DEFAULT_SIZE);

	dict *kdict = d->kmer_dict;

	/* sequence */
	char *rseq = read->sequence;
	char *conv_q = convert_quality_scores(read->qscore, read->length, 
			d->opt->qual_score_offset, bwidth, mp);

	state_nbhd_t *T_nbhds = mempool_nalloc(mp,
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
	//ksubstr_t obs_base = kmer_effective_base(rseq[tmin+kmer_len-1]);
	ksubstr_t obs_next_base = kmer_effective_base(rseq[tmin+kmer_len]);

	T_nbhds[tmin].size = hmm_load_preconstructed_nbhd(
			d, read->id, obs_kid, 
			obs_n_flag, conv_q, &T_nbhds[tmin], kmer_len, 1, kdict, mp, 
			&_uniq_tp);

	if (dump_trellis) {
		fprintf(d->log_fp, "NBHD	%zu	0	", read->id);
		dump_nbhd(&T_nbhds[tmin], d->log_fp);
		fprintf(d->log_fp, "\n");
	}

	sum_nbhd_size += T_nbhds[tmin].size;

	int wsize = d->opt->hd_window_size;
	uint8_t *win_dmax = &d->read_dmax[read->id * d->n_hd_windows];

	/* initial state distribution */
	size_t sz_base_emit_p_off = kmer_len * 20;
	size_t sz_qual_emit_p_off = kmer_len * qmax2;
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

		if (mode == PMR_MODE_ERROR_CORRECTION) {
			double obs_k_cexp = EXP(d->avg_initial_p[kmer_initial_index(obs_pk,
				d->uniq_kmer_init_off)]) * d->total_kmer_exp_counts;
			sum_kmer_cexp_cov += obs_k_cexp;
		}

		kbits_t _new_nf = ((obs_next_base >> 2) << 1) | 
			(obs_next_base >> 2);
		next_obs_n_flag = (next_obs_n_flag >> 2) | (_new_nf << shift);

		state_nbhd_t *curr_nbhd = &T_nbhds[t],
					 *next_nbhd = &T_nbhds[t+1];

		//const ksubstr_t dmax = d->opt->max_hamming_dist;
		ksubstr_t dmax = win_dmax[t / wsize];
		//dmax >>= 1;
		//if (dmax > 4) dmax = 4;

		int t_nbhd_size = curr_nbhd->size;

		total_nbhd_size += t_nbhd_size;

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

		if (_iter == 1) {
			sum_nbhd_size += tnext_nbhd_size;

			/*
			if (tnext_nbhd_size > _max_nbhd_size) {
				_max_nbhd_size = tnext_nbhd_size;
				fprintf(stderr, "Max nbhd size: %d\n", _max_nbhd_size);
			}
			*/
		}

		if (tnext_nbhd_size == 0 && _iter == 1) {
#if PMR_SHIFT_TMIN
			// try to shift the starting position to the right, by finding the
			// first "valid" kmer
			if (tmin < tmax - 1) {
				double thresh = d->opt->penalty_eta / log1p(1/d->opt->penalty_gamma);
				int valid_kmer_located = false;
				int s = tmin + 1;
				for (; s < tmax; ++s) {
					kbits_t _kid = kmer_seq_to_id(rseq + s, kmer_len);
					kmer_t *_pks = (kmer_t *) &(dict_find(kdict, 
								kbits_cast_to_ptr(_kid))->value);

					double avg_k_cexp = EXP(d->avg_initial_p[kmer_initial_index(_pks,
								d->uniq_kmer_init_off)]) * d->total_kmer_exp_counts;

					if (avg_k_cexp > thresh) {
						valid_kmer_located = true;
						break;
					}
				}

				if (valid_kmer_located) {
					*rerun = true;
					d->read_start_positions[read->id] = s;
				}
			}
			else {
				d->read_start_positions[read->id] = tmax;
			}
#endif

			if (!*rerun) {
				++_zero_likelihood;
				d->read_em_time_cov[(read->id << 1) + 1] = 1; 
			}

			goto destroy; 

			/* if somehow all pathways terminated up to this point, 
			 * treat the endpoint as tmax, and perform backward algorithm. */
			// break;
		}

		/* allocate memory given tnext_nbhd_size */
		char *pmem = mempool_nalloc(mp,
				nbhd_alloc_size(tnext_nbhd_size), 16);

		next_nbhd->size = tnext_nbhd_size;
		hmm_setup_nbhd_ptrs(next_nbhd, tnext_nbhd_size, pmem);

		/* FIXME: is there any means to wrap following code block into a 
		 * macro template? */

		int kmer_idx = 0, index_pfx = 0;

		//double t_max_alpha = tnext_max_alpha;

		tnext_max_alpha = NAN;
		tnext_argmax_alpha = 0;
		//int tnext_nonzero_states = 0;

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

				/* fixed on the above value, the emission
				 * density is also proportional to the expected relative 
				 * frequency of the observed kmer, \mu(\omega), 
				 * where \omega is the denoted as the observed kmer. */
				/*
				emit_dens = PRODUCT(emit_dens, 
						d->avg_initial_p[kmer_initial_index(obs_pk,
							d->uniq_kmer_init_off)]);
				*/
#endif

				double alpha = hmm_compute_alpha(
						ups_trans_packed, 
						kmer_idx,
						curr_nbhd->kmer_ptrs,
						t_states_alphas,
						t_kmer_trans_p,
						emit_dens,
						dns_base);

				// approximation.
				/*
				if (t > 0 && t < tmax - 1) {
					if (t_max_alpha - alpha > PMR_NBHD_THIN_THRESHOLD) {
						alpha = NAN;
					}
				}
				*/

#if PMR_EXTEND_READS
				//if (t == tmax - 1) {
					// for the last kmer, incorporate the probability of read
					// extension.
					double pext = d->prob_kmer_extension[
							kmer_initial_index(tnext_kmer, d->uniq_kmer_init_off)];
					alpha = PRODUCT(alpha, pext);

				if (t == tmax - 1) {
					// for the last kmer, incorporate the probability of read
					// initialize the backward algorithm recursion properly.
					next_nbhd->betas[index_sfx] = pext;
				}
#endif

				// if (!isnan(alpha)) tnext_nonzero_states++;

				update_max_quantity(alpha, index_sfx,
						tnext_max_alpha, tnext_argmax_alpha);

				next_nbhd->kmer_ptrs[index_sfx] = tnext_kmer;
				next_nbhd->kmer_trans_prob[index_sfx] = 
					tnext_kmer->unique_trans ? &_uniq_tp : 
					d->transition_p + d->nonuniq_kmers[tnext_kmer->index].poff;

				/* set up transition flag */
				bitarray_set_pwr2(next_nbhd->ba_kmer_trans_flag,
						index_sfx, 
						tnext_kmer->trans_flag, // IS_ZERO(alpha) ? 0 : tnext_kmer->trans_flag,
						2);

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

		/*
		if (t < tmax - 1) {
			for (int i = 0; i < tnext_nbhd_size; ++i) {
				if (tnext_max_alpha - next_nbhd->alphas[i] > 2 * PMR_NBHD_THIN_THRESHOLD) {
					next_nbhd->alphas[i] = NAN;
					bitarray_set_pwr2(next_nbhd->ba_kmer_trans_flag, i, 0, 2);
				}
			}
		}
		*/

#ifdef PMR_APPROX_HMM
		if (tnext_nonzero_states > max_approx_nbhd_size) {
			state_nbhd_t capped_nbhd;
			double sum_alphas = log_sum(tnext_nbhd_size, 
					next_nbhd->alphas, tnext_argmax_alpha, tnext_max_alpha);
			double kth_largest_alpha = quick_select_dbl(next_nbhd->alphas,
					tnext_nbhd_size, max_approx_nbhd_size);
			
			/* extra loop to count suffixes */
			int capped_actg_tc[4] = {0};
			int capped_cum_actg_tc[4] = {0};

			double inv_sum_approx_alphas = NAN;
			for (int i = 0; i < tnext_nbhd_size; i++) {
				double _alpha = next_nbhd->alphas[i];
				if (!isnan(_alpha) && _alpha >= kth_largest_alpha) {
					inv_sum_approx_alphas = log_add(inv_sum_approx_alphas, _alpha);

					register kbits_t st = next_nbhd->states_sorted_sfx[i] >> shift;
					++capped_actg_tc[st];
				}
			}

			capped_cum_actg_tc[0] = capped_actg_tc[0];
			for (int i = 1; i < 4; i++) {
				capped_cum_actg_tc[i] = capped_cum_actg_tc[i-1] +
					capped_actg_tc[i];
			}

			/* the resulting reduced neighborhood may be slightly larger than
			 * max_approx_nbhd_size if some alpha values tie with each other */
			int tnext_approx_nbhd_size = capped_cum_actg_tc[3];

			char *pcap_mem = mempool_nalloc(mp,
					nbhd_alloc_size(tnext_approx_nbhd_size), 16);

			capped_nbhd.size = tnext_approx_nbhd_size; 
			capped_nbhd.n_uniq_sfx = 0;
			hmm_setup_nbhd_ptrs(&capped_nbhd, tnext_approx_nbhd_size, pcap_mem);

			inv_sum_approx_alphas = -inv_sum_approx_alphas;

			int t_idx = 0, tn_aprx_pfx_idx = 0, tn_aprx_sfx_idx = 0;
			int tnext_pfx_idx = 0;
			for (int i = 0; i < n_distinct_suffixes; i++) {
				kbits_t _tf = bitarray_get_pwr2(curr_nbhd->ba_pfx_sfx_flags,
						i, 3);
				kbits_t ups_trans_packed = trans_flag2packed(_tf & 15);
				kbits_t dns_trans_packed = trans_flag2packed(_tf >> 4);

				register int n = dns_trans_packed & 7;
				register int j;
				register uint32_t _dns_tf = 0;

				for (j = 0, dns_trans_packed >>= 3; j < n; 
						++j, dns_trans_packed >>= 2) {
					kbits_t dns_base = dns_trans_packed & 3;

					int tnext_idx = next_nbhd->suffixes_order[tnext_pfx_idx];

					register double _alpha = next_nbhd->alphas[tnext_idx];
					if (!isnan(_alpha) && _alpha >= kth_largest_alpha) {
						int tn_aprx_sfx_idx = capped_cum_actg_tc[dns_base] - 
							capped_actg_tc[dns_base];
						--capped_actg_tc[dns_base];

						_dns_tf |= (1U << dns_base);	
						/* normalize alpha */
						capped_nbhd.alphas[tn_aprx_sfx_idx] = PRODUCT(_alpha,
								PRODUCT(sum_alphas, inv_sum_approx_alphas));
						/* copy verbatimly what have been set up previously */
						capped_nbhd.emit_dens[tn_aprx_sfx_idx] = next_nbhd->emit_dens[
							tnext_idx];
						capped_nbhd.base_emit_ctx[tn_aprx_sfx_idx] = 
							next_nbhd->base_emit_ctx[tnext_idx];
						capped_nbhd.qual_emit_ctx[tn_aprx_sfx_idx] = 
							next_nbhd->qual_emit_ctx[tnext_idx];
						capped_nbhd.kmer_ptrs[tn_aprx_sfx_idx] = 
							next_nbhd->kmer_ptrs[tnext_idx];
						capped_nbhd.kmer_trans_prob[tn_aprx_sfx_idx] = 
							next_nbhd->kmer_trans_prob[tnext_idx];
						capped_nbhd.suffixes_order[tn_aprx_pfx_idx] = tn_aprx_sfx_idx;

						/* set up transition flag */
						bitarray_set_pwr2(capped_nbhd.ba_kmer_trans_flag,
								tn_aprx_sfx_idx, 
								next_nbhd->kmer_ptrs[tnext_idx]->trans_flag,	
								2);
						/* set up Hamming distance */
						bitarray_set_pwr2(capped_nbhd.ba_hamming_dist, 
								tn_aprx_sfx_idx,
								bitarray_get_pwr2(next_nbhd->ba_hamming_dist, 
									tnext_idx, 3),
								3);

						++tn_aprx_pfx_idx;
					}
					++tnext_pfx_idx;
				}

				/* reset transition flag */
				bitarray_set_pwr2(curr_nbhd->ba_pfx_sfx_flags, 
						i, (_dns_tf << 4) | (_tf & 15), 3);
			}

			/* overwrite next_nbhd */
			T_nbhds[t+1] = capped_nbhd;
		}
#endif

		if (dump_trellis) {
			fprintf(d->log_fp, "NBHD	%zu	%d	", read->id, t + 1);
			dump_nbhd(next_nbhd, d->log_fp);
			fprintf(d->log_fp, "\n");
		}

		//obs_base = obs_next_base;
		obs_n_flag = next_obs_n_flag;
	}

	tmax = t;

	/* compute log-likelihood at the last position */
	double hmm_loglik = log_sum(T_nbhds[tmax].size, T_nbhds[tmax].alphas,
			tnext_argmax_alpha, tnext_max_alpha);

#if 0
#ifdef PMR_DUMP_NBHD_SIZE
	if (iter == 0) {
		float _ns[(tmax+1) * 3];

		for (int t = tmin; t <= tmax; t++) {
			kbits_t obs_kt = bitarray_get(observed_kmers, t << 1, kmer_len << 1);

			kmer_t *pk = (kmer_t *) &(dict_find(d->kmer_dict, 
						kbits_cast_to_ptr(obs_kt))->value);

			/* # of times a kmer and its reverse complement are observed */
			double k_cexp = EXP(d->avg_initial_p[kmer_initial_index(pk,
						d->uniq_kmer_init_off)]) * d->total_kmer_exp_counts;

			_ns[t*3] = (float) t;
			_ns[t*3+1] = (float) k_cexp;
			_ns[t*3+2] = (float) T_nbhds[t].size;
		}

		fwrite(_ns, sizeof(float), (tmax+1) * 3, d->nbhd_size_fp);
	}
#endif
#endif

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
		kbits_t xt_argmax = 4L;
		int yt = (int) conv_q[tmax_kp];

		// probability that the last base is correct.
		double prob_corr = 0.0;

		if (update_qual_score) {
			xt_argmax = kmer_effective_base(argmax_read->sequence[tmax_kp]);
		}

		for (int i = 0; i < tmax_nbhd_size; i++) {
			double norm_exp_count = EXP(PRODUCT(tmax_nbhd->alphas[i],
					-hmm_loglik));

			if (!update_qual_score) {
				int bctx = tmax_nbhd->base_emit_ctx[i]; 
				int qctx = tmax_nbhd->qual_emit_ctx[i]; 

				atomic_add_dbl(&tmax_base_emit_exp[bctx * 5 + xt],
						norm_exp_count);
				atomic_add_dbl(&tmax_qual_emit_exp[qctx * qmax + yt],
						norm_exp_count);
			}
			else {
				kbits_t st = tmax_nbhd->states_sorted_sfx[i] >> shift;
				if (st == xt_argmax) {
					prob_corr += norm_exp_count;
				}
			}
		}

		if (update_qual_score) {
			double pe = 1.0 - prob_corr;
			uint8_t q = 40;
			if (pe > 1e-4) {
				q = (uint8_t) floor(-10.0 * log(pe) / log(10.0));
			}
			argmax_read->qscore[tmax_kp] = q + d->opt->qual_score_offset;
		}

	//}
#endif
	}
	else {
		//d->read_start_positions[read->id] = tmax;
		++_zero_likelihood;
		d->read_em_time_cov[(read->id << 1) + 1] = 1; 
		goto destroy;
	}

	//state_nbhd_t *next_nbhd = &T_nbhds[tmax];
	/* set BETA's at the tmax to be 1.0 */
	//memset(next_nbhd->betas, 0, 
	//			next_nbhd->size * sizeof(*next_nbhd->betas));

	/* ---- backward algorithm ---- */
	/*
	for (int t = tmax - 1; t >= 0; t--) {
		printf("%d ", T_nbhds[t].size);
	}
	printf("\n");
	*/

	for (int t = tmax - 1; t >= tmin; t--) {
		state_nbhd_t *curr_nbhd = &T_nbhds[t],
					 *next_nbhd = &T_nbhds[t+1];

		int t_nbhd_size = curr_nbhd->size;

		/* expected counts and expected transitions, the gamma's and
		 * xi's in Rabiner's notation. */
		double *t_exp_counts = mempool_alloc(mp, sizeof(*t_exp_counts) *
				t_nbhd_size);
		double *t_exp_trans = mempool_alloc(mp, sizeof(*t_exp_trans) *
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

				double _alpha = curr_nbhd->alphas[t_kmer_idx];
#if PMR_EXTEND_READS
				kmer_t *pk = curr_nbhd->kmer_ptrs[t_kmer_idx];
				double pext = d->prob_kmer_extension[
						kmer_initial_index(pk, d->uniq_kmer_init_off)];
				_alpha = PRODUCT(_alpha, pext);
#endif

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
#endif

		/* --- update expected values (super redudant comment) --- */
//#pragma omp critical
//{
		if (!update_qual_score) {
			hmm_update_expected_vals(d, curr_nbhd, t,
					t_exp_counts, t_exp_trans,
					t_base_emit_exp, t_qual_emit_exp,
					kmer_len, ctx_shift, qmax, rseq, 
					conv_q, inv_sum_exp_count);
		}
		else {
			hmm_update_qual_scores(d, curr_nbhd, t, 
					read, argmax_read, t_exp_counts, inv_sum_exp_count);
		}
//}

	}

destroy:
	if (_iter == 1) {
		te = omp_get_wtime();
		d->read_em_time_cov[read->id << 1] = (te - ts);
		// d->read_em_time_cov[(read->id << 1) + 1] = (double) sum_nbhd_size / t;
		// d->read_em_time_cov[(read->id << 1) + 1] = sum_kmer_cexp_cov;
	}
	mempool_destroy(mp);
}

void hmm_write_to_fastq(data *d, void *fdata)
{
	char **fastq_output = (char **) fdata;
	FILE *stream = d->errors_fp ? d->errors_fp : stdout;

	for (int i = 0; i < PMR_BATCH_GROUP_SIZE; ++i) {
		if (fastq_output[i]) {
			fprintf(stream, "%s", fastq_output[i]);

			free(fastq_output[i]);
			fastq_output[i] = NULL;
		}
	}
}

void hmm_viterbi_wrapper(data *d, read_t *read, void *fdata)
{
	static char BASES[5] = {'A', 'C', 'T', 'G', 'N'};

	read_t *rorig = read_create();
	if (d->opt->update_qual_score) {
		rorig->id = read->id;
		rorig->length = read->length;
		rorig->sequence = calloc(rorig->length + 1, sizeof(char));
		rorig->qscore = calloc(rorig->length + 1, sizeof(char));

		memcpy(rorig->sequence, read->sequence, rorig->length);
		memcpy(rorig->qscore, read->qscore, rorig->length);
	}

	hmm_viterbi(d, read, (void *) 0);
	char **fastq_output = (char **) fdata;

	if (d->opt->viterbi_revcomp && d->read_start_positions[read->id] > 0) {
		// if a read's starting position has been shifted, try correct the
		// errors in the 5' region with a Viterbi on the reverse complementary
		// strand.

		d->preconstructed_nbhds[read->id].n = 0;

		read_t *r1 = read_create();

		r1->id = read->id;
		r1->length = read->length;
		r1->sequence = calloc(r1->length + 1, sizeof(char));
		r1->qscore = calloc(r1->length + 1, sizeof(char));

		for (int p = 0; p < r1->length; ++p) {
			r1->sequence[p] = BASES[
				base_to_rc_bits(read->sequence[r1->length - p - 1])];
			r1->qscore[p] = read->qscore[r1->length - p - 1];
		}

		hmm_viterbi(d, r1, (void *) 1);

		/* reverse complement again */
		for (int p = 0; p < r1->length; ++p) {
			read->sequence[p] = BASES[
				base_to_rc_bits(r1->sequence[r1->length - p - 1])];
			read->qscore[p] = r1->qscore[r1->length - p - 1];
		}

		free(r1->sequence);
		free(r1->qscore);
		read_destroy(r1);
	}

	if (d->opt->update_qual_score) {
		// FIXME: need to use the original read,
		// and pass the argmax sequence as a parameter.
		hmm_e_step(d, rorig, NULL, read);
	}

	read_destroy(rorig);

	if (d->opt->output_ordered == 0) {
		FILE *stream = d->errors_fp ? d->errors_fp : stdout;
		fprintf(stream, "@%s ZL:%d\n%.*s\n+\n%.*s\n",
				read->identifier, 
				(int) d->read_em_time_cov[(read->id << 1) + 1],
				read->length, read->sequence,
				read->length, read->qscore
			);
	}
	else {
		size_t buf_len = strlen(read->identifier) + 16 // @,+,ZL:x,\n x 4
			+ read->length * 2;
		fastq_output[read->batch_id] = calloc(buf_len, 1);
		snprintf(fastq_output[read->batch_id], buf_len, 
				"@%s ZL:%d\n%.*s\n+\n%.*s\n",
				read->identifier, 
				(int) d->read_em_time_cov[(read->id << 1) + 1],
				read->length, read->sequence,
				read->length, read->qscore);
	}
}

void hmm_viterbi_debug(data *d, read_t *read, char *target_seq)
{
	static char BASES[5] = {'A', 'C', 'T', 'G', 'N'};

	/* FIXME: code reusage? right now this file has too much redundancy... */
	int actg_trans_count[4];
	int cum_actg_trans_count[4];

	const int kmer_len = d->opt->kmer_length;
	const int num_kmers = read->length - kmer_len + 1;
	const int shift = (kmer_len - 1) << 1;
	const int qmax = d->opt->max_qual_score + 1;
	const int qmax2 = qmax << 1;
	const int bwidth = d->opt->qscore_bin_width;
	const double _uniq_tp = 0.0;

	// in reverse complementary decoding mode, must use the first kmer (which
	// is the last kmer on the original read).
	int tmin = d->read_start_positions[read->id];

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

	state_nbhd_t *T_nbhds = mempool_nalloc(mp,
			sizeof(*T_nbhds) * num_kmers, 16);

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
	//ksubstr_t obs_base = kmer_effective_base(rseq[tmin+kmer_len-1]);
	ksubstr_t obs_next_base = kmer_effective_base(rseq[tmin+kmer_len]);

	T_nbhds[tmin].size = hmm_load_preconstructed_nbhd(
			d, read->id, obs_kid, 
			obs_n_flag, conv_q, &T_nbhds[tmin], kmer_len, 1, kdict, mp, 
			&_uniq_tp);

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
		int tnext_target_idx = -1;

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

		kbits_t tnext_target_kid = kmer_seq_to_id(&target_seq[t+1], kmer_len);

		//kmer_t *obs_pk = (kmer_t *) &(dict_find(kdict,
		//			kbits_cast_to_ptr(obs_kid))->value);

		register kbits_t _new_nf = ((obs_next_base >> 2) << 1) | 
			(obs_next_base >> 2);
		next_obs_n_flag = (next_obs_n_flag >> 2) | (_new_nf << shift);

		state_nbhd_t *curr_nbhd = &T_nbhds[t],
					 *next_nbhd = &T_nbhds[t+1];

		dmax = win_dmax[t / wsize]; 


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
		char *pmem = mempool_nalloc(mp,
				nbhd_alloc_size(tnext_nbhd_size), 16);

		next_nbhd->size = tnext_nbhd_size;
		hmm_setup_nbhd_ptrs(next_nbhd, tnext_nbhd_size, pmem);

		if (tnext_nbhd_size == 0) {
			/* if somehow all pathways terminated up to this point, perform
			 * backtrack immediately */
			break;
		}

		/* FIXME: is there any means to wrap following code block into a 
		 * macro template? */
		int kmer_idx = 0; //, index_pfx = 0;

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

#endif

				/* we will be use betas as the psi's following Rabiner's
				 * notation. */
				double delta;
				int argmax_delta = hmm_compute_delta(
						ups_trans_packed, 
						kmer_idx,
						curr_nbhd->kmer_ptrs,
						t_states_deltas,
						t_kmer_trans_p,
						emit_dens,
						dns_base,
						&delta);

				/* if we are at the end of the read, 
				 * check if the kmer is specially flagged as a "last kmer",
				 * i.e. a kmer that cannot be further extended. 
				 * all kmer designated as the "last kmer" will be assigned zero
				 * likelihood */
				// if (t == (tmax - 1) && tnext_kmer->last_kmer) delta = NAN;

#if PMR_EXTEND_READS
				//if (t == tmax - 1) {
					// for the last kmer, incorporate the probability of read
					// extension.
					delta = PRODUCT(delta, d->prob_kmer_extension[
							kmer_initial_index(tnext_kmer, d->uniq_kmer_init_off)]);
				//}
#endif

				if (tnext_kid == tnext_target_kid) 
					tnext_target_idx = index_sfx;

				update_max_quantity(delta, index_sfx,
						tnext_max_delta, tnext_argmax_delta);

				next_nbhd->kmer_ptrs[index_sfx] = tnext_kmer;
				next_nbhd->kmer_trans_prob[index_sfx] = 
					tnext_kmer->unique_trans ? &_uniq_tp : 
					d->transition_p + d->nonuniq_kmers[tnext_kmer->index].poff;

				/* set up transition flag */
				bitarray_set_pwr2(next_nbhd->ba_kmer_trans_flag,
						index_sfx, tnext_kmer->trans_flag, 2);
				/* set up Hamming distance */
				bitarray_set_pwr2(next_nbhd->ba_hamming_dist, 
						index_sfx, 
						sfx_hd + ((mismatch_flag >> dns_base) & 1), 3);

				next_nbhd->states_sorted_sfx[index_sfx] = tnext_kid;
				next_nbhd->alphas[index_sfx] = delta;
				/* FIXME: consider use unnamed union to avoid type casting */
				next_nbhd->betas[index_sfx] = (double) argmax_delta;
				next_nbhd->emit_dens[index_sfx] = emit_dens;
			}

			kmer_idx += n_common_sfx_kmers;
		}

		//obs_base = obs_next_base;
		obs_n_flag = next_obs_n_flag;

		kbits_t argmax_kid = next_nbhd->states_sorted_sfx[tnext_argmax_delta]; 
		printf("[T%3d P%3d / S:%d]:%c\n * Max<%d>: %lu:%c/%zu (%g)",
				t+1, tnext_kp+1, next_nbhd->size, rseq[tnext_kp], tnext_argmax_delta,
				argmax_kid,
				BASES[argmax_kid >> shift],
				bitarray_get_pwr2(next_nbhd->ba_hamming_dist,
					tnext_argmax_delta, 3),
				next_nbhd->alphas[tnext_argmax_delta]);

		kbits_t tpacked = trans_flag2packed(next_nbhd->kmer_ptrs[tnext_argmax_delta]->trans_flag);
		double *_tp = next_nbhd->kmer_trans_prob[tnext_argmax_delta];
		int n_trans = tpacked & 7;
		tpacked >>= 3;
		printf(" [ ");
		for (int b = 0; b < n_trans; ++b, tpacked >>= 2) {
			printf("%c:%.3e ", BASES[tpacked & 3], EXP(_tp[b]));	
		}

		int btidx = (int) next_nbhd->betas[tnext_argmax_delta];
		printf("] ^%.3e <- %lu/%d (%lg) [ >%c:%.3e ]\n", exp(next_nbhd->emit_dens[tnext_argmax_delta]),
				curr_nbhd->states_sorted_sfx[btidx],
				kmer_hamming_dist(obs_kid, curr_nbhd->states_sorted_sfx[btidx]),
				curr_nbhd->alphas[btidx],
				BASES[argmax_kid >> shift], 
				EXP(curr_nbhd->kmer_trans_prob[btidx][base_to_transidx(
					curr_nbhd->kmer_ptrs[btidx]->trans_flag, argmax_kid >> shift)])
				);

		printf(" * Target: %lu:%c/%zu (%g)",
				tnext_target_kid, 
				BASES[tnext_target_kid >> shift],
				tnext_target_idx < 0 ? -1 : 
					bitarray_get_pwr2(next_nbhd->ba_hamming_dist, tnext_target_idx, 3),
				tnext_target_idx < 0 ? 0.0: next_nbhd->alphas[tnext_target_idx]
			  );

		if (tnext_target_idx >= 0) {
			tpacked = trans_flag2packed(next_nbhd->kmer_ptrs[tnext_target_idx]->trans_flag);
			_tp = next_nbhd->kmer_trans_prob[tnext_target_idx];
			n_trans = tpacked & 7;
			tpacked >>= 3;
			printf(" [ ");
			for (int b = 0; b < n_trans; ++b, tpacked >>= 2) {
				printf("%c:%.3e ", BASES[tpacked & 3], EXP(_tp[b]));	
			}

			int btidx = (int) next_nbhd->betas[tnext_target_idx];
			printf("] ^%.3e <- %lu/%d (%lg) [ >%c:%.3e ]", EXP(next_nbhd->emit_dens[tnext_target_idx]),
					curr_nbhd->states_sorted_sfx[btidx],
					kmer_hamming_dist(obs_kid, curr_nbhd->states_sorted_sfx[btidx]),
					curr_nbhd->alphas[btidx],
					BASES[tnext_target_kid >> shift], 
					EXP(curr_nbhd->kmer_trans_prob[btidx][base_to_transidx(
						curr_nbhd->kmer_ptrs[btidx]->trans_flag, tnext_target_kid >> shift)])
					);
		}
		printf("\n");
	}


	if (((early_decode == 0) && (d->opt->partial_decoding_lvl == 0)) && t < tmax) {
		goto destroy;
	}

	int end_of_read = (tmax == t);

	/* identify the most likely pathway at the last position */
	/* and perform partial decoding */
	kbits_t tfin_obs_kid = (t < tmax) ? obs_kid : 
		(obs_kid >> 2) | (obs_next_base << shift);

	tmax = t;

	if (tnext_argmax_delta >= 0 && !IS_ZERO(tnext_max_delta)) {
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

		//if (!end_of_read) 
			//fprintf(stream, "# ----------\n");
	}

destroy:
	mempool_destroy(mp);
}

/** 
 * hmm_viterbi
 *
 * Implementation of the Viterbi algorithm used for error decoding.
 *
 ****/
void hmm_viterbi(data *d, read_t *read, void *fdata)
{
	/* FIXME: code reusage? right now this file has too much redundancy... */
	int actg_trans_count[4];
	int cum_actg_trans_count[4];

	const int kmer_len = d->opt->kmer_length;
	const int num_kmers = read->length - kmer_len + 1;
	const int shift = (kmer_len - 1) << 1;
	const int qmax = d->opt->max_qual_score + 1;
	const int qmax2 = qmax << 1;
	const int bwidth = d->opt->qscore_bin_width;
	const double _uniq_tp = 0.0;

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

	state_nbhd_t *T_nbhds = mempool_nalloc(mp,
			sizeof(*T_nbhds) * num_kmers, 16);

	/*
	int wpos = retrieve_watch_pos(read->identifier);
	wpos -= (kmer_len);
	*/

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
	//ksubstr_t obs_base = kmer_effective_base(rseq[tmin+kmer_len-1]);
	ksubstr_t obs_next_base = kmer_effective_base(rseq[tmin+kmer_len]);

	T_nbhds[tmin].size = hmm_load_preconstructed_nbhd(
			d, read->id, obs_kid, 
			obs_n_flag, conv_q, &T_nbhds[tmin], kmer_len, 1, kdict, mp, 
			&_uniq_tp);

	/*
	T_nbhds[tmin].size = hmm_load_preconstructed_nbhd(
			d, obs_kid, 
			obs_n_flag, conv_q, &T_nbhds[tmin], kmer_len, 1, kdict, mp, 
			&_uniq_tp);
	*/

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

		//kmer_t *obs_pk = (kmer_t *) &(dict_find(kdict,
		//			kbits_cast_to_ptr(obs_kid))->value);

		/*
		if (t == wpos) {
			kmer_t *opk = &(dict_find(kdict, 
						kbits_cast_to_ptr(obs_kid))->value);
			double *tp = opk->unique_trans ? &_uniq_tp : 
				d->transition_p + d->nonuniq_kmers[opk->index].poff;

			kbits_t tfpk = trans_flag2packed(opk->trans_flag);
			int n_trans = tfpk & 7, i = 0;
			printf("%d %d %lu ", read->id, t + kmer_len + 1, obs_kid);
			for (tfpk >>= 3; i < n_trans; tfpk >>= 2, i++) {
				kbits_t base = tfpk & 3;
				printf("%lu %lg ", base, EXP(tp[i]));
			}
			printf("\n");
		}
		*/
		/*
		kmer_t *obs_next_kmer = dict_find(kdict,
				kbits_cast_to_ptr(obs_next_kid))->value;
		*/

		register kbits_t _new_nf = ((obs_next_base >> 2) << 1) | 
			(obs_next_base >> 2);
		next_obs_n_flag = (next_obs_n_flag >> 2) | (_new_nf << shift);

		state_nbhd_t *curr_nbhd = &T_nbhds[t],
					 *next_nbhd = &T_nbhds[t+1];

		//const kbits_t dmax = d->opt->max_hamming_dist;
		 dmax = win_dmax[t / wsize]; 


		kbits_t *t_states_sorted_sfx = T_nbhds[t].states_sorted_sfx;
		double *t_states_deltas = T_nbhds[t].alphas;
		double **t_kmer_trans_p = T_nbhds[t].kmer_trans_prob;

		/* ----- II. compute one iteration of Viterbi algorithm ----- */

		/* --- II a. set up kmer neighborhood at next position --- */
		memset(actg_trans_count, 0, sizeof(int) << 2);

		/*
		if (curr_nbhd->size <= 2) {
			int max_path_idx = (curr_nbhd->size == 1) ? 0 :
				(curr_nbhd->alphas[0] > curr_nbhd->alphas[1] ? 0 : 1);
			int _hd = bitarray_get_pwr2(curr_nbhd->ba_hamming_dist, 
					max_path_idx, 3);
			if (_hd >= dmax) ++dmax;

		}
		*/

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
		char *pmem = mempool_nalloc(mp,
				nbhd_alloc_size(tnext_nbhd_size), 16);

		next_nbhd->size = tnext_nbhd_size;
		hmm_setup_nbhd_ptrs(next_nbhd, tnext_nbhd_size, pmem);

		if (tnext_nbhd_size == 0) {
			/* if somehow all pathways terminated up to this point, perform
			 * backtrack immediately */
			break;
		}

		/* FIXME: is there any means to wrap following code block into a 
		 * macro template? */
		int kmer_idx = 0; //, index_pfx = 0;

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
				int argmax_delta = hmm_compute_delta(
						ups_trans_packed, 
						kmer_idx,
						curr_nbhd->kmer_ptrs,
						t_states_deltas,
						t_kmer_trans_p,
						emit_dens,
						dns_base,
						&delta);

				/* if we are at the end of the read, 
				 * check if the kmer is specially flagged as a "last kmer",
				 * i.e. a kmer that cannot be further extended. 
				 * all kmer designated as the "last kmer" will be assigned zero
				 * likelihood */
				// if (t == (tmax - 1) && tnext_kmer->last_kmer) delta = NAN;

				update_max_quantity(delta, index_sfx,
						tnext_max_delta, tnext_argmax_delta);

				next_nbhd->kmer_ptrs[index_sfx] = tnext_kmer;
				next_nbhd->kmer_trans_prob[index_sfx] = 
					tnext_kmer->unique_trans ? &_uniq_tp : 
					d->transition_p + d->nonuniq_kmers[tnext_kmer->index].poff;

				/* set up transition flag */
				bitarray_set_pwr2(next_nbhd->ba_kmer_trans_flag,
						index_sfx, tnext_kmer->trans_flag, 2);
				/* set up Hamming distance */
				bitarray_set_pwr2(next_nbhd->ba_hamming_dist, 
						index_sfx, 
						sfx_hd + ((mismatch_flag >> dns_base) & 1), 3);

				next_nbhd->states_sorted_sfx[index_sfx] = tnext_kid;
				next_nbhd->alphas[index_sfx] = delta;
				/* FIXME: consider use unnamed union to avoid type casting */
				next_nbhd->betas[index_sfx] = (double) argmax_delta;
				/* NOTE: do no need 
				next_nbhd->suffixes_order[index_pfx] = index_sfx;
				next_nbhd->states_sorted_pfx[index_pfx++] = tnext_kid;
				next_nbhd->emit_dens[index_sfx] = emit_dens;
				*/
			}

			kmer_idx += n_common_sfx_kmers;
		}

		//obs_base = obs_next_base;
		obs_n_flag = next_obs_n_flag;
	}

	if (((early_decode == 0) && (d->opt->partial_decoding_lvl == 0)) && t < tmax) {
		goto destroy;
	}

	int end_of_read = (tmax == t);

	/* identify the most likely pathway at the last position */
	/* and perform partial decoding */
	kbits_t tfin_obs_kid = (t < tmax) ? obs_kid : 
		(obs_kid >> 2) | (obs_next_base << shift);

	tmax = t;

	if (tnext_argmax_delta >= 0 && !IS_ZERO(tnext_max_delta)) {
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

		//if (!end_of_read) 
			//fprintf(stream, "# ----------\n");
	}

destroy:
	mempool_destroy(mp);
}

void hmm_init_model_params(data *d)
{
	int kmer_len = d->opt->kmer_length;
	//int penal_init_params = d->opt->penalize_init_params;

	/* --- emission --- */
#if !PMR_USE_PRIMITIVE_EMISSION
	INFO("Initiating emission probabilities...\n");
	const double init_qual_emit_p = LOG(1.0 / 
			(d->opt->max_qual_score + 1));
	const int qmax = (d->opt->max_qual_score + 1);
	const int qmax2 = (d->opt->max_qual_score + 1) << 1;

	int _tw = kmer_len + d->n_emit_windows;
	for (int t = 0; t < _tw; ++t) {
		int base_toff = 0, qual_toff = 0;

		if (t < kmer_len) {
			base_toff = t * 20;
			qual_toff = t * qmax2;

			/* - base emission - */
			for (int i = 0; i <= BASE_G; i++) {
				double *ctx_emit_p = &d->base_emission_p[base_toff + i * 5];
				for (int obs_base = 0; obs_base <= BASE_N; obs_base++) {
					if (i == obs_base) {
#if PMR_ALTERNATIVE_EMISSION
						ctx_emit_p[obs_base] = 0.0; 
#else
						ctx_emit_p[obs_base] = LOG((1.0 - d->opt->error_rate));
#endif
					}
					else {
						/* currently not allowing emitting N base */
#if PMR_ALTERNATIVE_EMISSION
						ctx_emit_p[obs_base] = -LOG(3.0); 
#else
						ctx_emit_p[obs_base] = LOG(d->opt->error_rate / 3);
#endif
					}
				}
			}

			/* - quality score emission */
			double *t_qual_emit_p = &d->qual_emission_p[qual_toff];
#if PMR_ALTERNATIVE_EMISSION
			for (int i = 0; i < qmax; i++) {
				double erate = POW(10.0, (double) i / -10.0);

				// special treatment for quality score 2, 
				// unless --phred=2 is specified by user.
				if (i == 2 && d->opt->qscore_as_phred < 2) 
					erate = d->opt->error_rate;

				t_qual_emit_p[i] = LOG(erate); 
				t_qual_emit_p[i + qmax] = LOG(1.0 - erate);
			}
#else
			for (int i = 0; i < qmax2; i++) {
				t_qual_emit_p[i] = init_qual_emit_p;
			}
#endif
		}
		else {
			int w = (t - kmer_len);
			base_toff = kmer_len * 20 + w * d->base_ctx_radix;
			for (int i = 0; i < d->n_base_context; i++) {
				double *ctx_emit_p = &d->base_emission_p[base_toff + i * 5];

				int ctx_base = i & 3;
				//int single_stranded = i / 4;
#if !PMR_ALTERNATIVE_EMISSION
				double error_rate = d->opt->error_rate;
#endif
				// double error_rate = single_stranded ? 0.75 : d->opt->error_rate;
				for (int obs_base = 0; obs_base <= BASE_N; obs_base++) {
					if (ctx_base == obs_base) {
#if PMR_ALTERNATIVE_EMISSION
						ctx_emit_p[obs_base] = 0.0; 
#else
						ctx_emit_p[obs_base] = LOG((1.0 - error_rate));
#endif
					}
					else {
#if PMR_ALTERNATIVE_EMISSION
						ctx_emit_p[obs_base] = -LOG(3.0);
#else
						ctx_emit_p[obs_base] = LOG(error_rate/ 3);
#endif
					}
				}
			}

			/* - quality score emission */
			double *t_qual_emit_p = &d->qual_emission_p[kmer_len * qmax2 + 
				w * d->qual_ctx_radix];
#if PMR_ALTERNATIVE_EMISSION
			for (int i = 0; i < qmax; i++) {
				double erate = POW(10.0, (double) i / -10.0);

				// special treatment for quality score 2.
				if (i == 2) 
					erate = d->opt->error_rate;

				t_qual_emit_p[i] = LOG(erate); 
				t_qual_emit_p[i + qmax] = LOG(1.0 - erate);
			}
#else
			for (int i = 0; i < d->qual_ctx_radix; i++) {
				t_qual_emit_p[i] = init_qual_emit_p;
			}
#endif
		}
	}
#endif

	/* --- initial distribution --- */
#if 0
	INFO("Initiating initial distribution parameters...\n");
	if (!penal_init_params) {
		int i = d->kmer_dict->used;
		dictentry *de = d->kmer_dict->table;
		for (; i > 0; de++) {
			if (de->value == NULL) continue;
			i--;

			kmer_t *pk = de->value;
			pk->init_p = LOG(pk->count / d->total_kmer_exp_counts);
		}
	}
	else {
		mstep_initial(d);
	}
#endif
}

void hmm_init_e_step(data *d)
{
	//int kmer_len = d->opt->kmer_length;
	int qmax2 = (d->opt->max_qual_score + 1) << 1;

#if PMR_USE_PRIMITIVE_EMISSION
	d->error_rate_emit_exp[0] = 0.0;
	d->error_rate_emit_exp[1] = 0.0;
#else
	//double *t_base_emit_exp;
	//double *t_qual_emit_exp;

	//int _tw = kmer_len + d->n_emit_windows;
	/* reset the emission densities table */

	size_t alloc_base_emit_slots = (4 * 5 * d->opt->kmer_length +
			d->n_emit_windows * d->base_ctx_radix) * 
		sizeof(*d->base_emission_p);
	size_t alloc_qual_emit_slots = (qmax2 * d->opt->kmer_length +
			d->n_emit_windows * d->qual_ctx_radix) * 
		sizeof(*d->qual_emission_p);

	memset(d->base_emit_exp, 0, 
			alloc_base_emit_slots * d->opt->n_threads);
	memset(d->qual_emit_exp, 0, 
			alloc_qual_emit_slots * d->opt->n_threads);

	/*
	for (int t = 0; t < _tw; ++t) {
		if (t < kmer_len) {
			t_base_emit_exp = &d->base_emit_exp[t * 20];
			t_qual_emit_exp = &d->qual_emit_exp[t * qmax2];

			for (int i = 0; i < 20; i++) {
				t_base_emit_exp[i] = 0.0;
			}

			for (int i = 0; i < qmax2; i++) {
				t_qual_emit_exp[i] = 0.0;
			}
		}
		else {
			int w = t - kmer_len;
			t_base_emit_exp = &d->base_emit_exp[kmer_len * 20 + 
				w * d->base_ctx_radix];
			t_qual_emit_exp = &d->qual_emit_exp[kmer_len * qmax2 +
				w * d->qual_ctx_radix];

			for (int i = 0; i < d->base_ctx_radix; i++) {
				t_base_emit_exp[i] = 0.0;
			}

			for (int i = 0; i < d->qual_ctx_radix; i++) {
				t_qual_emit_exp[i] = 0.0;
			}
		}
	}
	*/
	
#endif

	/* reset expected transitions to zero */ 
	memset(d->exp_trans_p, 0, d->n_total_trans * sizeof(*d->exp_trans_p));
}

//inline static double hmm_compute_beta(kbits_t dns_trans_packed,
double hmm_compute_beta(kbits_t dns_trans_packed,
		int dns_kmer_idx, kbits_t _tf, double *tnext_state_betas, 
		double *kmer_trans_p, double *kmer_exp_trans, 
		int64_t *tnext_sfx_order, double *emit_dens, double alpha)
{
	dbl_t tails[4] = {{NAN}};
	int n_dns_states = dns_trans_packed & 7;
	kbits_t dtp = dns_trans_packed >> 3;

	dbl_t max_tail = {.dbl = NAN};
	int argmax_tail = 0;
	
	int n = 0;
	//register kbits_t _tf = pk->trans_flag;
	for (int i = 0; i < n_dns_states; i++, dtp >>= 2) {
		kbits_t dns_base = dtp & 3;
		if ((_tf & (1 << dns_base)) == 0) continue;

		int beta_idx = tnext_sfx_order[dns_kmer_idx + i];
		register int dns_bidx = _mm_popcnt_u64(_tf & ((1 << dns_base) - 1));

		tails[n].dbl = PRODUCT(
				kmer_trans_p[n],
				PRODUCT(emit_dens[beta_idx], 
					tnext_state_betas[beta_idx]));

		/* if max_tail < tails[dns_base], ge_mask is set to 0, 
		 * otherwise UINT64_MAX */
		uint64_t ge_mask = (uint64_t) (isnan(max_tail.dbl) || 
				(max_tail.dbl < tails[n].dbl)) - 1UL;
		uint64_t inv_ge_mask = ~ge_mask;

		max_tail.u64 = (ge_mask & max_tail.u64) | 
			(inv_ge_mask & tails[n].u64);
		argmax_tail = (argmax_tail & ge_mask) | 
			(n & inv_ge_mask);

		kmer_exp_trans[dns_bidx] = PRODUCT(alpha, tails[n].dbl);
		++n;
	}

	double _expsum = 0.0; 
	for (int i = 0; i < n; i++) {
		register uint64_t argmax_mask = ((i == argmax_tail) - 1UL);
		register dbl_t summand = {.dbl = EXP(tails[i].dbl - max_tail.dbl)};
		summand.u64 &= argmax_mask;

		_expsum += summand.dbl; 
	}

	return PRODUCT(max_tail.dbl, log1p(_expsum));
}

/* hmm_compute_emit_dens
 * st: s_{[t]}, the k-th (last) base of the hidden state s_t.
 * xt: x_{[t]}, the k-th (last) base of the observed state x_t.
 */
//inline static double hmm_compute_emit_dens(
double hmm_compute_emit_dens(
		kbits_t st, kbits_t xt, int yt, kbits_t enf,
		double *base_emit_p, double *qual_emit_p)
{
	return PRODUCT(qual_emit_p[yt],  /* q_t() */
			base_emit_p[xt]);
			//hmm_normalize_emission(base_emit_p, enf, xt)); /* g_t() */
}

inline static double hmm_normalize_emission(double *emit_p, kbits_t enf, 
		kbits_t xt)
{
	if (enf == 15) return LOG(emit_p[xt]);
	else {
		double sum_emit_p = emit_p[BASE_N];
		register kbits_t en_packed = trans_flag2packed(enf);
		register int i = (int) en_packed & 7;
		for (en_packed >>= 3; i > 0; i--, en_packed >>= 2) {
			kbits_t base = en_packed & 3;
			sum_emit_p += emit_p[base];
		}

		return LOG(emit_p[xt] / sum_emit_p);
	}
}

double hmm_compute_alpha(kbits_t ups_trans_packed, int ups_kmer_idx,
		kmer_t **t_kmer_ptrs, double *t_state_alphas, double **t_kmer_trans_p, 
		double emit_dens, int dns_base)
{
	double summands[4] = {NAN};
	int n_ups_states = ups_trans_packed & 7;
	kbits_t utp = ups_trans_packed >> 3;
	dbl_t max_summand = {.dbl = NAN};
	int argmax_idx = 0;

	int n = 0;
	for (int i = 0; i < n_ups_states; i++, utp >>= 2) {
		register int _idx = ups_kmer_idx + i;
		
		register kbits_t _tf = t_kmer_ptrs[_idx]->trans_flag;
		if ((_tf & (1 << dns_base)) == 0)
			continue;

		register int dns_bidx = _mm_popcnt_u64(_tf & ((1 << dns_base) - 1));
		register dbl_t _summand = {.dbl = PRODUCT(t_state_alphas[_idx], 
					t_kmer_trans_p[_idx][dns_bidx])};

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

	return PRODUCT(emit_dens,
			PRODUCT(max_summand.dbl, log1p(_expsum)));
}

/*
 * @param ups_trans_packed	packed form of upstream transitions (no. transitions in lower 3 bits, bases in upper bits)
 * @param ups_kmer_idx
 * @param t_kmer_ptrs		array of kmers (kmer_t's)
 * @param t_state_deltas
 * @param t_kmer_trans_p
 * @param emit_dens
 * @param dns_base		next base
 * @param delta			value to compute
 */
double hmm_compute_delta(kbits_t ups_trans_packed, 
		int ups_kmer_idx, kmer_t **t_kmer_ptrs, double *t_state_deltas, 
		double **t_kmer_trans_p, double emit_dens, int dns_base, double *delta)
{
	dbl_t max_tail = {.dbl = NAN};
	int argmax_tail = 0;

	int n_ups_states = ups_trans_packed & 7;
	kbits_t utp = ups_trans_packed >> 3;	/* upstream bases */

	int n = 0;
	for (int i = 0; i < n_ups_states; i++, utp >>= 2) {
		register int _idx = ups_kmer_idx + i;

		/* not allowed transition from this upstream kmer */
		register kbits_t _tf = t_kmer_ptrs[_idx]->trans_flag;
		if ((_tf & (1 << dns_base)) == 0)
			continue;

		register int dns_bidx = _mm_popcnt_u64(_tf & ((1 << dns_base) - 1));
		register dbl_t _tail = {.dbl = PRODUCT(t_state_deltas[_idx], 
					t_kmer_trans_p[_idx][dns_bidx])};

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

	*delta = PRODUCT(max_tail.dbl, emit_dens);
	return argmax_tail;	
}

inline void hmm_correct_errors(data *d, read_t *read, int tmax,
		kbits_t st, kbits_t xt, kbits_t xt_n_flag, int kmer_len)
{
	static char BASES[5] = {'A', 'C', 'T', 'G', 'N'};

	//FILE *stream = d->errors_fp ? d->errors_fp : stdout;

	for (int i = 1; i <= kmer_len; i++) {
		if (xt_n_flag & 3) {
			/*
			fprintf(stream, "%s %d %c %c\n", read->identifier,
					tmax + i, BASES[st & 3], 'N');
			*/
			read->sequence[tmax + i - 1] = 'N';
		}
		else if ((xt & 3) != (st & 3)) {
			/*
			fprintf(stream, "%s %d %c %c\n", read->identifier,
					tmax + i, BASES[st & 3], BASES[xt & 3]);
			*/
			read->sequence[tmax + i - 1] = BASES[st & 3];
		}

		st >>= 2;
		xt >>= 2;
		xt_n_flag >>= 2;
	}
}

inline void hmm_correct_one_error(data *d, read_t *read, int t,
		kbits_t st, kbits_t xt)
{
	static char BASES[5] = {'A', 'C', 'T', 'G', 'N'};
	// FILE *stream = d->errors_fp ? d->errors_fp : stdout;

	if (st != xt) {
		/*
		fprintf(stream, "%s %d %c %c\n", read->identifier, 
				t + 1, BASES[st], BASES[xt]);
		*/
		read->sequence[t] = BASES[st];
	}
}

int hmm_count_distinct_suffixes(state_nbhd_t *nbhd, 
		int *actg_trans_count, kbits_t t_dmax, kbits_t obs_kid,
		kbits_t obs_n_flag, kbits_t obs_next_base, double *unique_trans_p)
{
	register int t_nbhd_size = nbhd->size;
	kbits_t *t_sts_sfx = nbhd->states_sorted_sfx;

	kbits_t prev_k_sfx = kbits_suffix(t_sts_sfx[0]);

	kbits_t pfx_trans_flag = 0, sfx_trans_flag = 0;
	kbits_t _tfhd = 0;
	int n_distinct_suffixes = 0;
	int n_trans = 0;
	int sfx_idx_lb = 0;

	/* FIXME: consider use larger step size for i to optimize cache access */
	for (int i = 0; i < t_nbhd_size; i++) {
		/* kmers, sorted by suffix */
		register kbits_t kid = t_sts_sfx[i];
		kbits_t k_sfx = kbits_suffix(kid);

		/* change in kmer suffix */
		if (k_sfx != prev_k_sfx) {
			bitarray_set_pwr2(nbhd->ba_distinct_sfx, 
					n_distinct_suffixes, 
					(i-sfx_idx_lb-1), 1);

			/* overwrite the original Hamming distance */
			bitarray_set_pwr2(nbhd->ba_hamming_dist, 
					n_distinct_suffixes,
					_tfhd >> 4, 3);

			bitarray_set_pwr2(nbhd->ba_pfx_sfx_flags, 
					n_distinct_suffixes,
					(sfx_trans_flag << 4) | pfx_trans_flag,
					3);

			for (int j = 0; j < 4; j++) {
				actg_trans_count[j] += ((sfx_trans_flag >> j) & 1);
			}

			++n_distinct_suffixes;

			sfx_trans_flag = 0;
			pfx_trans_flag = 0;
			sfx_idx_lb = i;
		}


		/* unconstrained transitions */
		kbits_t k_unc_tf = bitarray_get_pwr2(nbhd->ba_kmer_trans_flag,
				i, 2);
		kbits_t k_hamming_dist = bitarray_get_pwr2(nbhd->ba_hamming_dist, 
				i, 3);

		/* _tfhd: 8 bits packed integer in following format:
		 * HDHD TFTF
		 * where HDHD is kmer (k-1)-suffix hamming distance(kid, obs_kid)
		 */

		n_trans += trans_flag2packed(k_unc_tf) & 7;

		_tfhd = trans_admissible_transitions(kid, obs_kid,
				obs_n_flag, k_unc_tf, obs_next_base, k_hamming_dist, t_dmax);

		/* allowed transitions */
		kbits_t k_tf = _tfhd & 15;

		/* NOTE: As per 08/28/13, the implementation for transition
		 * distribution requires renormalization no more.
		 * Refer to the note on 08/28 for details. */
#if 0
		/* Code block deprecated */
		register int uniq_idx = trans_flag2index(k_tf);
		uintptr_t tf_eq_flag = ((uintptr_t) (k_unc_tf > k_tf)) - 1UL;
		nbhd->kmer_trans_prob[i] = (double *) (
				((uintptr_t) nbhd->kmer_trans_prob[i] & tf_eq_flag) |
				((uintptr_t) &unique_trans_p[uniq_idx] & ~tf_eq_flag));
#endif

		pfx_trans_flag |= (1UL << kbits_first_base(kid));
		sfx_trans_flag |= k_tf; 

		prev_k_sfx = k_sfx;
	}

	bitarray_set_pwr2(nbhd->ba_distinct_sfx, 
			n_distinct_suffixes, 
			(t_nbhd_size-sfx_idx_lb-1), 1);

	/* overwrite the original Hamming distance */
	bitarray_set_pwr2(nbhd->ba_hamming_dist, 
			n_distinct_suffixes,
			_tfhd >> 4, 3);

	bitarray_set_pwr2(nbhd->ba_pfx_sfx_flags, 
			n_distinct_suffixes,
			(sfx_trans_flag << 4) | pfx_trans_flag,
			3);

	for (int j = 0; j < 4; j++) {
		actg_trans_count[j] += ((sfx_trans_flag >> j) & 1);
	}

	++n_distinct_suffixes;
	nbhd->n_trans = n_trans;

	return n_distinct_suffixes;
}

int hmm_load_preconstructed_nbhd(data *d, int rid,
		kbits_t obs_1st_kid, kbits_t obs_n_flag, const char *conv_q, 
		state_nbhd_t *first_nbhd, int kmer_len, int compute_emit_dens, 
		dict *kdict, mempool_t *mp, const double *uniq_tp)
{
	/*
	kmer_t *ref_kmer = (kmer_t *) &(dict_find(kdict, 
			kbits_cast_to_ptr(obs_1st_kid))->value);
	*/


	int mode = d->opt->mode;
	size_t t_nbhd_size = 1;
	if (d->preconstructed_nbhds == NULL ||
			d->preconstructed_nbhds[rid].nbhd == NULL) {
		// will use the observed kmer
		char *pmem = mempool_nalloc(mp, nbhd_alloc_size(t_nbhd_size), 16);
		hmm_setup_nbhd_ptrs(first_nbhd, t_nbhd_size, pmem);

		kmer_t *nb_k = (kmer_t *) &(dict_find(d->kmer_dict, 
					kbits_cast_to_ptr(obs_1st_kid))->value);
		first_nbhd->states_sorted_sfx[0] = obs_1st_kid;
		first_nbhd->kmer_ptrs[0] = nb_k; 

		bitarray_set_pwr2(first_nbhd->ba_hamming_dist, 0, 0, 3);
		bitarray_set_pwr2(first_nbhd->ba_kmer_trans_flag, 0, 
				nb_k->trans_flag, 2);

		first_nbhd->kmer_trans_prob[0] = (nb_k->unique_trans ? 
				uniq_tp : d->transition_p + 
				(d->nonuniq_kmers[nb_k->index].poff));

		return 1;
	}

	int tmin_shifted = d->read_start_positions[rid] > 0;
	kmer_nbhd_t knbhd = d->preconstructed_nbhds[rid];

	if (tmin_shifted || knbhd.n == 0) {
		// use observed kmer
		knbhd.n = 1;
		knbhd.nbhd[0] = (kmer_t *) &(dict_find(d->kmer_dict, 
					kbits_cast_to_ptr(obs_1st_kid))->value);
	}

	t_nbhd_size = knbhd.n;
	char *pmem = mempool_nalloc(mp, nbhd_alloc_size(t_nbhd_size), 16);
	hmm_setup_nbhd_ptrs(first_nbhd, t_nbhd_size, pmem);

	first_nbhd->size = t_nbhd_size;
	for (int i = 0; i < t_nbhd_size; i++) {
		register kmer_t *nb_k = knbhd.nbhd[i];
		if (compute_emit_dens) {
			double emit_prob = hmm_compute_kmer_emission(d, obs_1st_kid, 
					obs_n_flag, nb_k, kmer_len, conv_q);
			if (mode == PMR_MODE_ERROR_CORRECTION) {
				double init_p = d->avg_initial_p[kmer_initial_index(nb_k,
						d->uniq_kmer_init_off)];
				first_nbhd->alphas[i] = PRODUCT(emit_prob, init_p); //PRODUCT(nb_k->init_p, emit_prob);
			}
			first_nbhd->alphas[i] = emit_prob; //PRODUCT(nb_k->init_p, emit_prob);
		}
		else {
			first_nbhd->alphas[i] = NAN;
		}

		kbits_t kid = kbits_id_of(nb_k);

		first_nbhd->states_sorted_sfx[i] = kid;
		first_nbhd->kmer_ptrs[i] = nb_k; 
		first_nbhd->kmer_trans_prob[i] = (nb_k->unique_trans ? 
				uniq_tp : d->transition_p + 
				(d->nonuniq_kmers[nb_k->index].poff));

		bitarray_set_pwr2(first_nbhd->ba_hamming_dist, i, 
				kmer_n_hamming_dist(obs_1st_kid, kid, obs_n_flag), 3);
		bitarray_set_pwr2(first_nbhd->ba_kmer_trans_flag, i, 
				nb_k->trans_flag, 2);
	}

	return t_nbhd_size;
}

#if 0
static inline int hmm_load_preconstructed_nbhd(data *d, kbits_t obs_1st_kid, 
		kbits_t obs_n_flag, const char *conv_q, state_nbhd_t *first_nbhd, 
		int kmer_len, int compute_emit_dens, dict *kdict, mempool_t *mp)
{
	kmer_t *ref_kmer = dict_find(kdict, 
			kbits_cast_to_ptr(obs_1st_kid))->value;

	int t_nbhd_size = ref_kmer->n_neighbors;
	kmer_t **kmer_nbhd = ref_kmer->data;

	char *pmem = mempool_nalloc(mp, nbhd_alloc_size(t_nbhd_size), 16);
	hmm_setup_nbhd_ptrs(first_nbhd, t_nbhd_size, pmem);

	first_nbhd->size = t_nbhd_size;
	for (int i = 0; i < t_nbhd_size; i++) {
		register kmer_t *nb_k = kmer_nbhd[i];
		if (compute_emit_dens) {
			double emit_prob = hmm_compute_kmer_emission(d, obs_1st_kid, 
					obs_n_flag, nb_k, kmer_len, conv_q);
			first_nbhd->alphas[i] = PRODUCT(nb_k->init_p, emit_prob);
		}
		else {
			first_nbhd->alphas[i] = NAN;
		}

		first_nbhd->states_sorted_sfx[i] = nb_k->id;
		first_nbhd->kmer_ptrs[i] = nb_k; 
		first_nbhd->kmer_trans_prob[i] = nb_k->transition_p;

		bitarray_set_pwr2(first_nbhd->ba_hamming_dist, i, 
				kmer_n_hamming_dist(obs_1st_kid, nb_k->id, obs_n_flag), 2);
		bitarray_set_pwr2(first_nbhd->ba_kmer_trans_flag, i, 
				nb_k->trans_flag, 2);
	}

	return t_nbhd_size;
}
#endif

static void hmm_update_qual_scores(data *d, state_nbhd_t *curr_nbhd,
		int t, read_t *read, read_t *argmax_read, double *t_exp_counts,
		double inv_sum_exp_count)
{
	const int kmer_len = d->opt->kmer_length;
	const int shift = (kmer_len - 1) << 1;
	int t_kp = t + kmer_len - 1;
	int n_distinct_suffixes = curr_nbhd->n_uniq_sfx;
	int t_kmer_idx = 0;

	double prob_corr[kmer_len];
	memset(prob_corr, 0, kmer_len * sizeof(double));

	for (int i = 0; i < n_distinct_suffixes; i++) {
		kbits_t _tf = bitarray_get_pwr2(curr_nbhd->ba_pfx_sfx_flags, i, 3);
		kbits_t ups_trans_packed = trans_flag2packed(_tf & 15);

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
				norm_exp_count = EXP(_norm_cnt);
			}

			if (t > 0) {
				kbits_t st = ref_kid >> shift;
				// the error corrected base
				kbits_t xt_argmax = kmer_effective_base(argmax_read->sequence[t_kp]);

				if (xt_argmax == st)
					prob_corr[0] += norm_exp_count;
			}
			else {
				for (int l = 0; l < kmer_len; l++, ref_kid >>= 2) {
					kbits_t si = ref_kid & 3;
					kbits_t xi_argmax = kmer_effective_base(argmax_read->sequence[l]);

					if (xi_argmax == si)
						prob_corr[l] += norm_exp_count;
				}
			}

			t_kmer_idx++;
		}
	}

	int lmax = t > 0 ? 1 : kmer_len;
	int toff = (t > 0) ? t_kp : 0;
	for (int l = 0; l < lmax; ++l) {
		double pe = 1.0 - prob_corr[l];

		uint8_t q = 40;
		if (pe > 1e-4) {
			q = (uint8_t) floor(-10.0 * log(pe) / log(10.0));
		}

		argmax_read->qscore[toff + l] = q + d->opt->qual_score_offset;
	}
}

static void hmm_update_expected_vals(data *d, state_nbhd_t *curr_nbhd, 
		int t, double *t_exp_counts, double *t_exp_trans, 
		double *t_base_emit_exp, double *t_qual_emit_exp,
		int kmer_len, int ctx_shift, int qmax, 
		const char *rseq, const char *conv_q,
		double inv_sum_exp_count)
{
	int t_kp = t + kmer_len - 1;
#if PMR_USE_PRIMITIVE_EMISSION
	const int shift = (kmer_len - 1) << 1;
#endif

	int t_kmer_idx = 0;
	int n_distinct_suffixes = curr_nbhd->n_uniq_sfx;
	double *kmer_exp_trans = t_exp_trans;

	// thread-local sum
	// double *bctx_sum_ = calloc(20, sizeof(double));
	// double *qctx_sum_ = calloc(qmax << 1, sizeof(double)); 

	for (int i = 0; i < n_distinct_suffixes; i++) {
		kbits_t _tf = bitarray_get_pwr2(curr_nbhd->ba_pfx_sfx_flags,
				i, 3);
		kbits_t ups_trans_packed = trans_flag2packed(_tf & 15);

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
					bitarray_get_pwr2(curr_nbhd->ba_kmer_trans_flag, 
						t_kmer_idx, 2));

			//kbits_t k_trans_packed = trans_flag2packed(ref_k->trans_flag);

			double *ref_exp_trans_p = d->exp_trans_p + 
				(ref_k->unique_trans ? ref_k->index :
				 d->nonuniq_kmers[ref_k->index].poff);

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
					atomic_add_dbl(&ref_exp_trans_p[b], EXP(_norm_cnt));
					//ref_exp_trans_p[b] += EXP(_norm_cnt);
				//}
			}

			/* -- update expected emissions -- */
			if (t > 0) {
				kbits_t xt = kmer_effective_base(rseq[t_kp]);
				int yt = (int) conv_q[t_kp];

#if PMR_USE_PRIMITIVE_EMISSION
				kbits_t st = ref_kid >> shift;
				atomic_add_dbl(&d->error_rate_emit_exp[(xt != st)],
						norm_exp_count);
#else
				int bctx = curr_nbhd->base_emit_ctx[t_kmer_idx];
				int qctx = curr_nbhd->qual_emit_ctx[t_kmer_idx];

				// FIXME: slow!
				atomic_add_dbl(&t_base_emit_exp[bctx * 5 + xt], norm_exp_count);
				atomic_add_dbl(&t_qual_emit_exp[qctx * qmax + yt], norm_exp_count);
				// bctx_sum_[bctx * 5 + xt] += norm_exp_count;
				// qctx_sum_[qctx * qmax + yt] += norm_exp_count;
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
					atomic_add_dbl(&t_base_emit_exp[l * 20 + si * 5 + xi], 
							norm_exp_count);
					
					atomic_add_dbl(&t_qual_emit_exp[((l << 1) + (si == xi)) * qmax + yi],
							norm_exp_count);

					/*
					d->base_emit_exp[l * 20 + si * 5 + xi] += norm_exp_count;
					d->qual_emit_exp[((l << 1) + (si == xi)) * qmax + yi] +=
						norm_exp_count;
					*/
					
#endif
				}
			}

			t_kmer_idx++;
			//kmer_exp_trans += ref_k->n_trans;
			kmer_exp_trans += k_trans_packed & 7; 
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

	free(bctx_sum_);
	free(qctx_sum_);
	*/
}


double hmm_compute_kmer_emission(data *d, kbits_t obs_kid, 
		kbits_t obs_n_flag, kmer_t *pk, int kmer_len, const char *conv_q)
{
	double emit_dens = 0.0;
	kbits_t kid = kbits_id_of(pk);

#if PMR_USE_PRIMITIVE_EMISSION
	double _er = LOG(d->base_error_rate/3), _ner = LOG(1.0-d->base_error_rate);
	for (int t = 0; t < kmer_len; t++, kid >>= 2, obs_kid >>= 2) {
		register kbits_t xt = obs_kid & 3, st = kid & 3;
		emit_dens += ((xt == st) ? _ner : _er);
	}
#else
	int qmax = d->opt->max_qual_score + 1;
	for (int t = 0; t < kmer_len; t++, kid >>= 2, obs_kid >>= 2) {
		register kbits_t xt = obs_kid & 3, st = kid & 3;
		double *t_base_emit_p = &d->base_emission_p[t * 20 + st * 5];

		emit_dens += t_base_emit_p[xt]; //hmm_normalize_emission(t_base_emit_p, 15, xt);
		emit_dens += d->qual_emission_p[((t << 1) + (xt == st)) * qmax + 
			conv_q[t]];
	}
#endif

	return emit_dens;
}

// static inline void hmm_setup_nbhd_ptrs(state_nbhd_t *nbhd, int n, 
inline void hmm_setup_nbhd_ptrs(state_nbhd_t *nbhd, int n, 
		char *pmem)
{
	memset(pmem, 0, nbhd_alloc_size(n));

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
}

static void hmm_dump_model(data *d)
{
	uint64_t n_kmers = d->kmer_dict->used;
	dictentry *de = d->kmer_dict->table;

	for (int i = 0; i < n_kmers; de++) {
		if (de->key == NULL) continue;
		i++;

		kmer_t *pk = (kmer_t *) &de->value;
		kbits_t kid = kbits_id_of(pk);

		double *cur_tp = d->transition_p + (pk->unique_trans ?
				0 : d->nonuniq_kmers[pk->index].poff);
		double mu = d->avg_initial_p[kmer_initial_index(pk,
				d->uniq_kmer_init_off)];

		fprintf(d->params_fp, "%lu %d %lg ", kid, pk->label, mu);

		int idx = 0;
		for (int j = 0; j < 4; ++j) {
			if ((pk->trans_flag & (1U << j)) == 0) fprintf(d->params_fp, "nan ");
			else if (pk->unique_trans) fprintf(d->params_fp, "0.0 ");
			else {
				fprintf(d->params_fp, "%lg ", cur_tp[idx++]);
			}
		}
		fprintf(d->params_fp, "\n");
	}
}

#if 0
static void hmm_dump_model(data *d)
{
	const size_t n_kmers = d->kmer_dict->used;
	const size_t n_nonzero_kmers = n_kmers - d->n_zeroed_out_kmer;
	const int kmer_len = d->opt->kmer_length;
	const int qmax = d->opt->max_qual_score;

	dumpfile_header_t hdr;

	hdr.kmer_length = d->opt->kmer_length;
	hdr.max_qscore = d->opt->max_qual_score;
	hdr.max_read_length = d->max_read_length;
	hdr.max_hamming_dist = d->opt->max_hamming_dist;
	hdr.loglik = d->loglik;
	hdr.penalty_eta = d->opt->penalty_eta;
	hdr.penalty_gamma = d->opt->penalty_gamma;
	hdr.base_ctx_radix = d->base_ctx_radix;
	hdr.qual_ctx_radix = d->qual_ctx_radix;
	hdr.n_uniq_trans_kmer = d->n_uniq_trans_kmer; 
	hdr.n_multi_trans_kmer = n_kmers - hdr.n_uniq_trans_kmer - 
		d->n_zeroed_out_kmer;

	memcpy(hdr.identifier, "PMRPARAM", 8);
	strncpy(hdr.base_emit_cond_expr, d->opt->base_emit_cond_expr, 127);
	strncpy(hdr.qual_emit_cond_expr, d->opt->qual_emit_cond_expr, 127);

	size_t sz_trans_flags = sizeof(uint64_t) * BITS_TO_U64(
			n_nonzero_kmers << 2);
	size_t sz_trans_p = sizeof(double) * d->n_composite_trans;
	size_t sz_b_emit_p = sizeof(double) * (kmer_len * 20 +
			(hdr.max_read_length - kmer_len) * hdr.base_ctx_radix);
	size_t sz_q_emit_p = sizeof(double) * (kmer_len * (qmax << 1) +
			(hdr.max_read_length - kmer_len) * hdr.qual_ctx_radix);
	size_t sz_emit_p = sz_b_emit_p + sz_q_emit_p;

	hdr.off_trans_p = sizeof(hdr) + sz_trans_flags;
	hdr.off_emit_p = hdr.off_trans_p + sz_trans_p;
	hdr.off_kmers = hdr.off_emit_p + sz_emit_p;

	uint64_t *k_trans_flag = malloc(sz_trans_flags);
	double *k_trans_p = malloc(sz_trans_p);
	dumpfile_kmer_t *k_list = malloc(sizeof(*k_list) * 
			(n_kmers - d->n_zeroed_out_kmer));

	dictentry *de = d->kmer_dict->table;
	int tpidx = 0;
	int kidx = 0;
	for (int i = 0; i < n_kmers; de++) {
		if (de->value == NULL) continue;

		++i;
		kmer_t *pk = de->value;
		if (pk->trans_flag > 0) {
			uint64_t t_packed = trans_flag2packed(pk->trans_flag);
			int n_trans = t_packed & 7;
			/* only dump values for composite transitions probabilities 
			 * (i.e. more than 1 downstream transitions) */
			if (n_trans > 1) {
				int j = 0;
				for (t_packed >>= 3; j < n_trans; j++, t_packed >>= 2) {
					k_trans_p[tpidx++] = pk->transition_p[t_packed & 3];
				}
			}

			bitarray_set_pwr2(k_trans_flag, kidx, pk->trans_flag, 2);

			k_list[kidx].id = pk->id;
			k_list[kidx].exp_count = pk->exp_count;

			++kidx;
		}
	}

	/* header */
	fwrite(&hdr, sizeof(hdr), 1, d->params_fp);
	/* transition flags */
	fwrite(k_trans_flag, sz_trans_flags, 1, d->params_fp);
	/* transition probabilities */
	fwrite(k_trans_p, sz_trans_p, 1, d->params_fp);
	/* emission probabilities */
	fwrite(d->base_emission_p, sz_b_emit_p, 1, d->params_fp);
	fwrite(d->qual_emission_p, sz_q_emit_p, 1, d->params_fp);
	/* kmer list */
	fwrite(k_list, sizeof(dumpfile_kmer_t), (n_kmers - d->n_zeroed_out_kmer),
			d->params_fp);

	free(k_trans_flag);
	free(k_trans_p);
	free(k_list);
}
#endif

inline char *convert_quality_scores(const char *qscore, int len, 
		int offset, int bwidth, mempool_t *mp)
{
#ifdef __SSE2__
	/* FIXME: cache this result maybe */
	int i, w_m128 = (len >> 4) + ((len & 0xF) ? 1 : 0);
	char *conv_qscore = mempool_nalloc(mp, sizeof(*conv_qscore) * len, 16);
	__m128i *src_ptr = (__m128i *) qscore;
	__m128i *dst_ptr = (__m128i *) conv_qscore;
	__m128i qual_offset = _mm_set1_epi8(offset);
	//__m128i bin_width = _mm_set1_epi8(bwidth);
	__m128i xmm;

	for (i = 0; i < w_m128; i++){
		xmm = _mm_loadu_si128(src_ptr);
		xmm = _mm_sub_epi8(xmm, qual_offset);
		//xmm = _mm_div_epi8(xmm, bin_width);
		_mm_store_si128(dst_ptr, xmm);

		++dst_ptr;
		++src_ptr;
	}
#endif

	return conv_qscore;
}

void hmm_read_starting_pos(data *d, read_t *read, void *fdata)
{
	int kmer_len = d->opt->kmer_length;
	int tmax = read->length - kmer_len;

	double thresh = d->opt->penalty_eta / log1p(1/d->opt->penalty_gamma);

	/* convert string literal to numeric representation */
	kbits_t observed_kmers[BITS_TO_U64(read->length << 1)];
	read_seq_to_numeric(read->sequence, observed_kmers, read->length, kmer_len);

	int t = 0;
	for (; t <= tmax; ++t) {
		kbits_t obs_kid = bitarray_get(observed_kmers, t << 1, kmer_len << 1);
		kmer_t *pk = (kmer_t *) &(dict_find(d->kmer_dict, 
					kbits_cast_to_ptr(obs_kid))->value);

		double avg_k_cexp = EXP(d->avg_initial_p[kmer_initial_index(pk,
				d->uniq_kmer_init_off)]) * d->total_kmer_exp_counts;

		if (avg_k_cexp > thresh) break;
	}

	d->read_start_positions[read->id] = t;
}

void hmm_dump_strandness(data *d, read_t *read, void *fdata)
{
	FILE *fpout = fdata;

	int kmer_len = d->opt->kmer_length;
	int tmax = read->length - kmer_len;

	/* convert string literal to numeric representation */
	kbits_t observed_kmers[BITS_TO_U64(read->length << 1)];
	read_seq_to_numeric(read->sequence, observed_kmers, read->length, kmer_len);

	size_t n_u64 = BITS_TO_U64(tmax);
	kbits_t strandness[n_u64];
	memset(strandness, 0, n_u64 * sizeof(kbits_t));

	int t = 0;
	for (; t <= tmax; ++t) {
		kbits_t obs_kid = bitarray_get(observed_kmers, t << 1, kmer_len << 1);
		kmer_t *pk = (kmer_t *) &(dict_find(d->kmer_dict, 
					kbits_cast_to_ptr(obs_kid))->value);
		
		int label = pk->label;
		bitarray_set_pwr2(strandness, t, label, 1);
	}

	fwrite((void *) &strandness[0], sizeof(uint64_t), n_u64, fpout);
}

