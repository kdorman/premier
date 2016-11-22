#include "read.h"
#include "kmer.h"

void read_qual_select(data *d, read_t *read, void *fdata)
{
	int *rsel = (int *) fdata;
	int qoff = d->opt->qual_score_offset;
	char *qscore = read->qscore;
	double sum_qscore = 0.0;

	int is_above_thresh = 1;
	int read_len = read->length;
	for (int t = 0; t < read_len; ++t) {
		int q = qscore[t] - qoff;

		sum_qscore += q;
		is_above_thresh &= (q >= 30);
	}

	rsel[read->id] = (sum_qscore / read->length) >= 30;
	// RULE 2: must exceed Q=30 at all loci
	// rsel[read->id] = is_above_thresh;
}

/** Decompose a read into kmers, record kmer counts, construct
 * kmer neighborhood, and add kmer into the k-spectrum unless
 * the kmer has undetermined base(s).
 */
void read_parse(data *d, read_t *tmp_read, void *fdata)
{
	int kmer_len = d->opt->kmer_length;
	int read_len = tmp_read->length;
	int tmax = read_len - kmer_len;
	int qual_offset = d->opt->qual_score_offset;
	int max_qual = d->opt->max_qual_score;
	int shift = (kmer_len - 1) << 1;

	char *rseq = tmp_read->sequence;
	char *qscore = tmp_read->qscore;

	for (int p = 0; p < read_len; ++p) {
		if (qscore[p] < qual_offset)
			qual_offset = qscore[p];
		if (qscore[p] > max_qual)
			max_qual = qscore[p];
	}
	
	if (tmp_read->length > d->max_read_length)
		d->max_read_length = tmp_read->length;

	//d->total_loci_count += read_len;

	kbits_t kmer_id = kmer_seq_to_id(rseq, kmer_len);
	kbits_t kmer_rc_id = kmer_reverse_complement(kmer_id, kmer_len);
	kbits_t kmer_n_flag = kmer_n_base_flag(rseq, kmer_len);

	kbits_t obs_next_base = kmer_effective_base(rseq[kmer_len]);

	/*
	DOUBLE *ptr_qual_freq = d->qual_freq;
	for (i = 0; i < read_len; i++) {
		ptr_qual_freq[(int) qscore[i] - qual_offset] += 1;
	}
	*/

	/* increment transition count of forward transition, and
	 * set rc flag on reverse complement (no transition) for
	 * first kmer in read */
	if (kmer_n_flag == 0) {
//#pragma omp critical
//{
		kmer_incr_count_by_id(d, kmer_id, obs_next_base, FLAG_1ST_KMER);
		/* first position is the the last kmer on reverse strand */
		kmer_incr_count_by_id(d, kmer_rc_id, 0x4, FLAG_SET_RC);
				//FLAG_LAST_KMER | FLAG_SET_RC);
//}
	}

	/* for each subsequent transition */
	for (int t = 1; t < tmax; t++) {
		kbits_t kth_n_flag = (rseq[t - 1 + kmer_len] >> 3) & 1;
		kbits_t kth_base = base_to_bits(rseq[t - 1 + kmer_len]);
		obs_next_base = kmer_effective_base(rseq[t + kmer_len]);
		kbits_t rc_next_base = base_to_rc_bits(rseq[t - 1]);

		register kbits_t _nf = ((kbits_t) kth_n_flag << 1) | kth_n_flag;
		kmer_n_flag = kbits_suffix(kmer_n_flag) | (_nf << shift);

		kmer_id = kbits_suffix(kmer_id) | ((kbits_t) kth_base << shift);
		kmer_rc_id = kmer_reverse_complement(kmer_id, kmer_len);

		/* skip any kmer with at least one 'N' base */
		if (kmer_n_flag) continue;

		/* increment transition count on forward transition, and
		 * increment transition count on reverse transition, setting 
		 * reverse complement flag if first time this rc(kmer) seen. */
//#pragma omp critical
//{
		kmer_incr_count_by_id(d, kmer_id, obs_next_base, FLAG_NORMAL_KMER);
		kmer_incr_count_by_id(d, kmer_rc_id, rc_next_base, 
				FLAG_NORMAL_KMER | FLAG_SET_RC);
//}
		/*
		kmer_incr_count_by_id(d, kmer_rc_id, 0x4, //rc_next_base, 
				FLAG_NORMAL_KMER | FLAG_SET_RC);
		*/

		//printf("%d %d\n", tmp_read->id, t);
	}

	/* last kmer */
	kbits_t kth_n_flag = (rseq[tmax + kmer_len - 1] >> 3) & 1;
	register kbits_t _nf = ((kbits_t) kth_n_flag << 1) | kth_n_flag;

	kmer_n_flag = kbits_suffix(kmer_n_flag) | (_nf << shift);

	if (kmer_n_flag == 0) {
		kbits_t kth_base = (rseq[tmax + kmer_len - 1] >> 1) & 3;
		kmer_id = kbits_suffix(kmer_id) | ((kbits_t) kth_base << shift);
		kmer_rc_id = kmer_reverse_complement(kmer_id, kmer_len);

//#pragma omp critical
//{
		kmer_incr_count_by_id(d, kmer_id, 0x4, FLAG_LAST_KMER);

		kmer_incr_count_by_id(d, kmer_rc_id, base_to_rc_bits(rseq[tmax - 1]),
				FLAG_1ST_KMER | FLAG_SET_RC);
//}

		/* 0x4 is an "invalid" base */
		/*
		kmer_incr_count_by_id(d, kmer_rc_id, 0x4, //base_to_rc_bits(rseq[tmax - 1]),
				FLAG_1ST_KMER | FLAG_SET_RC);
		*/
	}

	d->opt->qual_score_offset = qual_offset;
	d->opt->max_qual_score = max_qual;
}	/* read_parse */
