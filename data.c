#include "data.h"
#include "trans.h" 
#include "mstep.h"
#include "atomic.h"
#include "iodata.h"

#include <signal.h>

static const char* file_name = "data.c";

#if 0
static int cmp_kmer_ptr_sfx(const void *ka, const void *kb)
{
	if ((*(kmer_t **) ka)->id < (*(kmer_t **) kb)->id) {
		return -1;
	}
	else {
		return 1;
	}
}
#endif

static void weigh_observed_transitions(data *d, read_t *read, void *fdata)
{
	const int kmer_len = d->opt->kmer_length;
	const int shift = (kmer_len - 1) << 1;
	const int qoff = d->opt->qual_score_offset;
	int tmax = read->length - kmer_len;

	double *_weighted_exp_trans_p = (double *) fdata;

	double weight = 1.0;
	kbits_t obs_kid = kmer_seq_to_id(read->sequence, kmer_len);
	kmer_t *obs_pk = kptr_locate(d, obs_kid);

	int poff = obs_pk->unique_trans ? obs_pk->index : 
		d->nonuniq_kmers[obs_pk->index].poff;
	double *obs_tcnt = d->exp_trans_p + poff;
	double *weighted_tcnt = _weighted_exp_trans_p + poff;

	for (int i = 0; i <= kmer_len; ++i) {
		weight *= qscore_to_prob(read->qscore[i], qoff);
	}
	// weight = qscore_to_prob(read->qscore[kmer_len], qoff);

	kbits_t obs_next_base = base_to_bits(read->sequence[kmer_len]);

	int tidx = base_to_transidx(obs_pk->trans_flag, obs_next_base);
	atomic_add_dbl(&weighted_tcnt[tidx], weight);
	atomic_add_dbl(&obs_tcnt[tidx], 1.0);

	for (int t = 1; t < tmax; ++t) {
		int tnext_kp = t + kmer_len;

		obs_kid = (obs_kid >> 2) | (obs_next_base << shift);
		obs_pk = kptr_locate(d, obs_kid);
		poff = obs_pk->unique_trans ? 
			obs_pk->index : d->nonuniq_kmers[obs_pk->index].poff;
		obs_tcnt = d->exp_trans_p + poff;
		weighted_tcnt = _weighted_exp_trans_p + poff;

		obs_next_base = base_to_bits(read->sequence[tnext_kp]);
		tidx = base_to_transidx(obs_pk->trans_flag, obs_next_base);

		/*
		weight = qscore_to_prob(read->qscore[tnext_kp], qoff);
		*/
		weight *= (qscore_to_prob(read->qscore[tnext_kp], qoff) /
				qscore_to_prob(read->qscore[t-1], qoff));

		atomic_add_dbl(&weighted_tcnt[tidx], weight);
		atomic_add_dbl(&obs_tcnt[tidx], 1.0);
	}
}

static inline double determine_threshold(int *kcnt_histo, double gam,
		int halve_threshold, double *thresh)
{
	int idx_valley = 1;
	for (int mul = idx_valley + 1; mul <= 1000; ++mul) {
		if (kcnt_histo[mul] > kcnt_histo[idx_valley]) break;
		++idx_valley;
	}
	//++idx_valley;

	// find the peak of the kmer histogram
	int idx_peak = idx_valley;
	int max_kcnt = kcnt_histo[idx_peak];
	for (int mul = idx_peak + 1; mul <= 1000; ++mul) {
		if (kcnt_histo[mul] > max_kcnt) {
			max_kcnt = kcnt_histo[mul];
			idx_peak = mul;
		}
	}

	// locate the valley
	int min_kcnt = kcnt_histo[idx_valley];
	for (int mul = idx_valley + 1; mul < idx_peak; ++mul) {
		if (kcnt_histo[mul] < 100) continue;
		
		if (kcnt_histo[mul] < min_kcnt) {
			min_kcnt = kcnt_histo[mul];
			idx_valley = mul;
		}
	}

	double _thresh = idx_valley == 1 ? 1.0 : ((double) idx_valley / (
			halve_threshold ? 2.0 : 1.0));
	if (idx_peak == 1000 || (max_kcnt - min_kcnt < 1000)) _thresh = 1.0;

	*thresh = _thresh;
	return ceil((_thresh + 0.1) * LOG(1.0+1.0/gam) + 0.5);
}

static inline uint64_t kmer_ptr2flag(data *d, kmer_t *pk,
		int *pcov, int *pnonuniq, int *pnonuniq_2s)
{
	*pcov = d->total_kmer_exp_counts *
		EXP(d->avg_initial_p[kmer_initial_index(pk, d->uniq_kmer_init_off)]);
	*pnonuniq = pk->n_trans > 1;
	*pnonuniq_2s = pk->n_trans > 1;

	return pk->trans_flag;
}

static inline uint64_t kmer_counter2flag(uint64_t kcnter, int *pcov, 
		int *pnonuniq, int *pnonuniq_2s)
{
	int cov = 0, n_nonzero = 0, n_nonzero_2s = 0;
	uint64_t tf = 0;
	uint64_t alt_tf = (kcnter >> 4) & 15;
	for (int base = 0; base < 4; ++base) {
		uint64_t _cnt = (kcnter >> ((8 + 14 * base)) & ((1 << 14) - 1));
		n_nonzero += (_cnt > 0);
		n_nonzero_2s += ((_cnt > 0) | ((alt_tf >> base) & 1));
		tf |= (_cnt > 0) << base;
		cov += _cnt; // kmer coverage
	}

	*pcov = cov;
	*pnonuniq = (n_nonzero > 1);
	*pnonuniq_2s = (n_nonzero_2s > 1);

	// combine the transition flag from the other strand.
	return tf | alt_tf;
}

static int masked_kmer_id_comp(const void *a, const void *b, const void *pmask)
{
	register kbits_t mask = ((mask_t *) pmask)->pmask;
	int shift = ((mask_t *) pmask)->mask_shift;
	register kbits_t id_a = (*(kbits_t *) a);
	register kbits_t id_b = (*(kbits_t *) b);
	register kbits_t mid_a = id_a & mask, mid_b = id_b & mask;

	kbits_t id_diff = ((~mask & id_a) >> shift) - ((~mask & id_b) >> shift);

	/* set to all 1 if id_1 != id_b, all 0 otherwise */
	register kbits_t non_equal_mask = (kbits_t) (mid_a == mid_b) - 1UL;
	register kbits_t _masked_sign = (int) (((mid_a - mid_b) >> 63) << 31) + 1;
	return (non_equal_mask & _masked_sign) | (~non_equal_mask & id_diff);
}

/* adjust_kmer_labeling
 *
 * assigns labels to the transition probabilities so that only probs. assigned
 * with label 0 are considered as free parameters. 
 * we require that all outgoing transitions from any given kmer must have the
 * same label, thus, the labeling of transitions is equivalent to the labeling
 * of kmers, which is performed by this function.
 *
 * the labeling is completed using bipartite graph partitioning. 
 * see manuscript for more details.
 */

void adjust_kmer_labeling(data *d)
{
	int n_kmers = d->kmer_dict->used;
	int kmer_len = d->opt->kmer_length;
	int shift = (kmer_len - 1) * 2;
	dictentry *de = d->kmer_dict->table;

	/* reset the "is labeled" flags for all kmers */
	for (int i = 0; i < n_kmers; de++) {
		if (de->key == NULL) continue;

		uint64_t kcnter = (uint64_t) de->value;
		/* reset the flag (lowest bit):
		 * 0 -> not processed
		 * 1 -> processed
		 * */
		kcnter &= ~1UL;

		// if a kmer has ever been marked as a "normal" kmer,
		// disqualify it as the "last" kmer
		if (kcnter & FLAG_NORMAL_KMER) kcnter &= ~4UL;

		// TODO: identify bipartite graphs upfront, 
		// so that we can perform graph partitioning later in parallel fashion?

		*((uint64_t *) &de->value) = kcnter;

		i++;
	}

	de = d->kmer_dict->table;
	// there can be up to 8 candidate kmers for labeling, 4 on each "strand".

	/* adjust flags */
	for (int i = 0; i < n_kmers; de++) {
		if (de->key == NULL) continue;

		kbits_t kid = (kbits_t) de->key - 1;
		uint64_t kcnter = (uint64_t) de->value;

		/* already labeled */
		if (kcnter & FLAG_LABELED) {
			i++;
			continue;
		}

		uint32_t graph_nodes = 0;
		uint32_t uv_trans_flags = 0;
		int uv_kcov[8] = {0};
		int uv_nonuniq[8] = {0};
		// to tell if a kmer has non-unique transitions by combining
		// information from both strands 
		int uv_nonuniq_2s[8] = {0};
		int nU = 0, nV = 0; // the two disjoint sets of kmers in the bipartite graph
		dictentry *lbl_cands[8] = {NULL};

		/* set U: */
		for (kbits_t base = 0; base < 4; ++base) {
			kbits_t _ku = (kid & ~3UL) | base;
			dictentry *_deu = dict_find(d->kmer_dict,
					kbits_cast_to_ptr(_ku));

			if (_deu != NULL) {
				++nU;
				lbl_cands[base] = _deu;
				graph_nodes |= (1 << base);
				uv_trans_flags |= (kmer_counter2flag((uint64_t) _deu->value, 
						&uv_kcov[base], &uv_nonuniq[base], &uv_nonuniq_2s[base]) << (base * 4));
			}
		}

		/* set V: */
		for (kbits_t base = 0; base < 4; ++base) {
			kbits_t _kv = kmer_reverse_complement((kid >> 2) | (base << shift), 
					kmer_len);

			kbits_t cbase = complementary_base(base);

			dictentry *_dev = dict_find(d->kmer_dict,
					kbits_cast_to_ptr(_kv));

			if (_dev != NULL) {
				++nV;
				lbl_cands[4 + cbase] = _dev;
				graph_nodes |= (1 << (4 + cbase));
				uv_trans_flags |= (kmer_counter2flag((uint64_t) _dev->value,
						&uv_kcov[4 + cbase], &uv_nonuniq[4 + cbase], &uv_nonuniq_2s[4 + base]) << ((4 + cbase) * 4));
			}
		}

		/* --- graph partitioning --- */
		while (graph_nodes) {
			// we arbitrarily call the first four bases/kmers as "forward"
			int _u = 0, _v = 0;
			// kmers that are never directly observed (in the k-spectrum only
			// through its reverse complement)
			int _phantom_u = 0, _phantom_v = 0;
			int n_fwd_nonuniq = 0, n_rev_nonuniq = 0;
			int n_fwd_nonuniq_2s = 0, n_rev_nonuniq_2s = 0;
			int fwd_cov = 0, rev_cov = 0;

			uint32_t visited = 0, pending = 0;

			int seed = 0;
			for (seed = 0; seed < 8; ++seed) if (graph_nodes & (1 << seed)) break;

			if (seed < 4) {
				++_u;
				_phantom_u += (uv_kcov[seed] == 0);
				fwd_cov += uv_kcov[seed];
				n_fwd_nonuniq += uv_nonuniq[seed];
				n_fwd_nonuniq_2s += uv_nonuniq_2s[seed];
			}
			else {
				++_v;
				_phantom_v += (uv_kcov[seed] == 0);
				rev_cov += uv_kcov[seed];
				n_rev_nonuniq += uv_nonuniq[seed];
				n_rev_nonuniq_2s += uv_nonuniq_2s[seed];
			}

			visited |= (1 << seed);
			pending |= (complementary_trans_flag(((uv_trans_flags >> (4 * seed)) & 15)) 
				<< ((seed < 4) * 4)) & graph_nodes;

			do {
				for (seed = 0; seed < 8; ++seed) if (pending & (1 << seed)) break;
				pending &= ~(1UL << seed);
				visited |= (1 << seed);

				if (seed < 4) {
					++_u;
					_phantom_u += (uv_kcov[seed] == 0);
					fwd_cov += uv_kcov[seed];
					n_fwd_nonuniq += uv_nonuniq[seed];
					n_fwd_nonuniq_2s += uv_nonuniq_2s[seed];
				}
				else {
					++_v;
					_phantom_v += (uv_kcov[seed] == 0);
					rev_cov += uv_kcov[seed];
					n_rev_nonuniq += uv_nonuniq[seed];
					n_rev_nonuniq_2s += uv_nonuniq_2s[seed];
				}

				pending |= (complementary_trans_flag(((uv_trans_flags >> (4 * seed)) & 15)) 
						<< ((seed < 4) * 4)) & graph_nodes & (~visited);
			} while(pending);

			/* visited now contains a connected component in the bipartite graph */
			int label = (n_fwd_nonuniq < n_rev_nonuniq);

			if (n_fwd_nonuniq == n_rev_nonuniq) {
				label = _u > _v;
			}

			if ((n_fwd_nonuniq == n_rev_nonuniq) && (_u == _v)) {
				if (_u == 1 && (_phantom_u != _phantom_v)) {
					label = _phantom_u < _phantom_v;
				}
				else {
					label = fwd_cov > rev_cov;
				}
			}

			/*
			if (n_fwd_nonuniq == n_rev_nonuniq) {
				label = fwd_cov > rev_cov;
			}
			*/

			// this is a pathological case because we have non-unique
			// transitions in both directions. if there are errors at both
			// sides, we can only remove errors at one side. :<
			//int ambiguous_lbl = (n_fwd_nonuniq_2s > 0 && n_rev_nonuniq_2s > 0);
			int ambiguous_lbl = (n_fwd_nonuniq > 0 && n_rev_nonuniq > 0);

			for (int base = 0; base < 8; ++base) {
				if ((visited & (1 << base)) == 0) continue;

				dictentry *deb = lbl_cands[base];
				uint64_t *_cnter = (uint64_t *) &(deb->value);

				int set_label = (base < 4)? label : 1-label;
				uint64_t _lbled_cnter = (*_cnter | FLAG_LABELED) & ~(8UL | 2UL);
				_lbled_cnter |= (set_label << 3) | (ambiguous_lbl << 1);
				*_cnter = _lbled_cnter;
			}

			graph_nodes &= ~visited;
		}

		i++;
	}
}

/*
 * Re-adjust kmer labeling, presumably after a few iterations of M-step.
 */
void readjust_kmer_labeling(data *d)
{
	int kmer_len = d->opt->kmer_length;
	int shift = (kmer_len - 1) * 2;
	dictentry *de = d->kmer_dict->table;

#if 0
	int n_kmers = d->kmer_dict->used;
	/* reset the "is labeled" flags for all kmers */
	for (int i = 0; i < n_kmers; de++) {
		if (de->key == NULL) continue;

		uint64_t kcnter = (uint64_t) de->value;
		/* reset the flag (lowest bit):
		 * 0 -> not processed
		 * 1 -> processed
		 * */
		kcnter &= ~1UL;

		// if a kmer has ever been marked as a "normal" kmer,
		// disqualify it as the "last" kmer
		if (kcnter & FLAG_NORMAL_KMER) kcnter &= ~4UL;

		// TODO: identify bipartite graphs upfront, 
		// so that we can perform graph partitioning later in parallel fashion?

		*((uint64_t *) &de->value) = kcnter;

		i++;
	}
	de = d->kmer_dict->table;
#endif

	// there can be up to 8 candidate kmers for labeling, 4 on each "strand".

	/* adjust flags */
	int n_readjusted = 0;
	int64_t n_size = d->kmer_dict->sizemask + 1;
	for (int64_t i = 0; i < n_size; ++i) {
		dictentry *de = &d->kmer_dict->table[i];
		if (de->value == NULL) continue;

		kmer_t *pk = (kmer_t *) &(de->value);
		kbits_t kid = kbits_id_of(pk);

		/* already relabeled */
		if (pk->dirty) continue;

		uint32_t graph_nodes = 0;
		uint32_t uv_trans_flags = 0;
		int uv_kcov[8] = {0};
		int uv_nonuniq[8] = {0};
		// to tell if a kmer has non-unique transitions by combining
		// information from both strands 
		int uv_nonuniq_2s[8] = {0};
		int nU = 0, nV = 0; // the two disjoint sets of kmers in the bipartite graph
		kmer_t *lbl_cands[8] = {NULL};

		/* set U: */
		for (kbits_t base = 0; base < 4; ++base) {
			kbits_t _ku = (kid & ~3UL) | base;

			dictentry *_deu = dict_find(d->kmer_dict,
					kbits_cast_to_ptr(_ku));

			if (_deu != NULL) {
				++nU;
				kmer_t *_pku = (kmer_t *) &(_deu->value); 
				lbl_cands[base] = _pku;
				graph_nodes |= (1 << base);
				uv_trans_flags |= (kmer_ptr2flag(d, _pku, 
							&uv_kcov[base], &uv_nonuniq[base], &uv_nonuniq_2s[base])) << (base * 4);
			}
		}

		/* set V: */
		for (kbits_t base = 0; base < 4; ++base) {
			kbits_t _kv = kmer_reverse_complement((kid >> 2) | (base << shift), 
					kmer_len);

			kbits_t cbase = complementary_base(base);

			dictentry *_dev = dict_find(d->kmer_dict,
					kbits_cast_to_ptr(_kv));

			if (_dev != NULL) {
				++nV;
				kmer_t *_pkv = (kmer_t *) &(_dev->value);
				lbl_cands[4 + cbase] = _pkv;

				graph_nodes |= (1 << (4 + cbase));
				uv_trans_flags |= (kmer_ptr2flag(d, _pkv,
						&uv_kcov[4 + cbase], &uv_nonuniq[4 + cbase], &uv_nonuniq_2s[4 + base]) << ((4 + cbase) * 4));
			}
		}

		/* --- graph partitioning --- */
		while (graph_nodes) {
			// we arbitrarily call the first four bases/kmers as "forward"
			int _u = 0, _v = 0;
			// kmers that are never directly observed (in the k-spectrum only
			// through its reverse complement)
			int _phantom_u = 0, _phantom_v = 0;
			int n_fwd_nonuniq = 0, n_rev_nonuniq = 0;
			int n_fwd_nonuniq_2s = 0, n_rev_nonuniq_2s = 0;
			int fwd_cov = 0, rev_cov = 0;

			uint32_t visited = 0, pending = 0;

			int seed = 0;
			for (seed = 0; seed < 8; ++seed) if (graph_nodes & (1 << seed)) break;

			if (seed < 4) {
				++_u;
				_phantom_u += (uv_kcov[seed] == 0);
				fwd_cov += uv_kcov[seed];
				n_fwd_nonuniq += uv_nonuniq[seed];
				n_fwd_nonuniq_2s += uv_nonuniq_2s[seed];
			}
			else {
				++_v;
				_phantom_v += (uv_kcov[seed] == 0);
				rev_cov += uv_kcov[seed];
				n_rev_nonuniq += uv_nonuniq[seed];
				n_rev_nonuniq_2s += uv_nonuniq_2s[seed];
			}

			visited |= (1 << seed);
			pending |= (complementary_trans_flag(((uv_trans_flags >> (4 * seed)) & 15)) 
				<< ((seed < 4) * 4)) & graph_nodes;

			do {
				for (seed = 0; seed < 8; ++seed) if (pending & (1 << seed)) break;
				pending &= ~(1UL << seed);
				visited |= (1 << seed);

				if (seed < 4) {
					++_u;
					_phantom_u += (uv_kcov[seed] == 0);
					fwd_cov += uv_kcov[seed];
					n_fwd_nonuniq += uv_nonuniq[seed];
					n_fwd_nonuniq_2s += uv_nonuniq_2s[seed];
				}
				else {
					++_v;
					_phantom_v += (uv_kcov[seed] == 0);
					rev_cov += uv_kcov[seed];
					n_rev_nonuniq += uv_nonuniq[seed];
					n_rev_nonuniq_2s += uv_nonuniq_2s[seed];
				}

				pending |= (complementary_trans_flag(((uv_trans_flags >> (4 * seed)) & 15)) 
						<< ((seed < 4) * 4)) & graph_nodes & (~visited);
			} while(pending);

			/* visited now contains a connected component in the bipartite graph */
			int label = (n_fwd_nonuniq < n_rev_nonuniq);

			if (n_fwd_nonuniq == n_rev_nonuniq) {
				label = _u > _v;
			}

			if ((n_fwd_nonuniq == n_rev_nonuniq) && (_u == _v)) {
				if (_u == 1 && (_phantom_u != _phantom_v)) {
					label = _phantom_u < _phantom_v;
				}
				else {
					label = fwd_cov > rev_cov;
				}
			}

			/*
			if (n_fwd_nonuniq == n_rev_nonuniq) {
				label = fwd_cov > rev_cov;
			}
			*/

			// this is a pathological case because we have non-unique
			// transitions in both directions. if there are errors at both
			// sides, we can only remove errors at one side. :<
			//int ambiguous_lbl = (n_fwd_nonuniq_2s > 0 && n_rev_nonuniq_2s > 0);
			int ambiguous_lbl = (n_fwd_nonuniq > 0 && n_rev_nonuniq > 0);

			for (int base = 0; base < 8; ++base) {
				if ((visited & (1 << base)) == 0) continue;

				kmer_t *pkb = lbl_cands[base];

				int set_label = (base < 4)? label : 1-label;
				if (pkb->label != set_label)
					++n_readjusted;

				pkb->label = set_label;
				pkb->dirty = 1;
				pkb->ambig_labeled = ambiguous_lbl;
			}

			graph_nodes &= ~visited;
		}
	}

	INFO_TEE(d->log_fp,
			"Labels readjusted: %d\n", n_readjusted);

	for (int64_t i = 0; i < n_size; ++i) {
		dictentry *de = &d->kmer_dict->table[i];
		if (de->value == NULL) continue;

		kmer_t *pk = (kmer_t *) &(de->value);

		pk->dirty = 0;
	}
}

/** allocate_trans_param
 * Allocate transition distribution parameters according to
 * number of transitions etc */
void allocate_trans_params(data *d, int halve_threshold)
{
	double exp_trans_p[4] = {0.0};
	//double thresh = d->opt->penalty_eta / log1p(1/d->opt->penalty_gamma);
	//double minthresh = (thresh < 1) ? 1.0 : thresh;

	size_t n_uniq_kmers = 0;
	size_t n_nonuniq_kmers = 0;
	size_t n_total_trans = 0;
	int n_kmers = d->kmer_dict->used;
	int kmer_len = d->opt->kmer_length;
	dictentry *de = d->kmer_dict->table;

	/* rough approximation */
	size_t n_tp = d->n_total_trans - n_kmers + d->n_reads;
	if (n_kmers > d->n_total_trans) {
		n_tp = d->n_total_trans;
	}

	d->transition_p = calloc((n_tp + 1), sizeof(*d->transition_p));
	double *uniq_kmer_cnts = calloc(n_kmers, sizeof(*uniq_kmer_cnts));
 
	double *kmer_trans_p = d->transition_p + 1;

	double total_trans_cnt = 0.0;
	size_t _uidx = 0;
	for (int i = 0; i < n_kmers; de++) {
		if (de->key == NULL) continue;

		uint64_t kcnter = (uint64_t) de->value;
		uint32_t alt_trans_flag = (kcnter >> 4) & 15;

		uint32_t orig_tf = 0;
		register uint32_t _tf = 0;
		register int n_nonzero_trans = 0;
		double sum_trans_cnt = 0.0;
		for (int base = 0; base < 4; base++) {
			uint32_t trans_cnt = (kcnter >> (8 + base * 14)) & ((1 << 14) - 1);
			exp_trans_p[n_nonzero_trans] = (double) trans_cnt;
			sum_trans_cnt += exp_trans_p[n_nonzero_trans];

#ifdef PMR_COMBINE_TRANS_FLAG
			/* if reverse complementary strand supports an unobserved
			 * transition, set the flag but set transition count still to zero.
			 */
			if (trans_cnt > 0 || (alt_trans_flag & (1 << base))) {
				_tf |= (1 << base);
				n_nonzero_trans++;
				if (trans_cnt > 0) orig_tf |= (1 << base);
			}
#else
			if (trans_cnt > 0) {
				_tf |= (1 << base);
				n_nonzero_trans++;
			}
#endif
		}

		/* check if the kmer only appears at the 3' end of a read */
		uint32_t is_last_kmer = ((kcnter & FLAG_LAST_KMER) >> 2) &&
			(orig_tf == 0);

		kmer_t km = {
			.trans_flag = _tf & 15, 
			.unique_trans = (n_nonzero_trans <= 1),
			.label = (kcnter >> 3) & 1,
			.n_trans = n_nonzero_trans,
			.prev_n_trans = n_nonzero_trans,
			.dirty = 1,
#if PMR_READ_EXTENSION_TWO_STRANDS
			.alt_trans_flag = _tf & 15,
#else
			.alt_trans_flag = orig_tf,
#endif
			//.alt_trans_flag = alt_trans_flag, 
			// .valid = sum_trans_cnt > minthresh,//sum_trans_cnt > n_nonzero_trans,
			.coverage = sum_trans_cnt > ((1 << 10) - 1) ? 
				((1 << 10) - 1) : sum_trans_cnt,
			.last_kmer = is_last_kmer, 
			.ambig_labeled = (kcnter >> 1) & 1,
			.single_stranded = 0,
			.reserved = 1, /* set to 1 to avoid the whole 64bit being zero */
			.index = 0
		};

		/* replace the kmer counter with kmer_t */
		*((kmer_t *) &de->value) = km;

		if (n_nonzero_trans == 1) {
			uniq_kmer_cnts[_uidx++] = sum_trans_cnt;
			n_uniq_kmers++;
			/* trans. prob. for unique transition is always 1.0 hence do 
			 * not need a dedicated parameter */
		}
		else if (n_nonzero_trans > 1) {
			n_nonuniq_kmers++;

			for (int j = 0; j < n_nonzero_trans; j++) {
				*kmer_trans_p = exp_trans_p[j];
				++kmer_trans_p;
			}

			n_total_trans += n_nonzero_trans;
		}

		total_trans_cnt += sum_trans_cnt;

		i++;
	}

	/* number of expected transitions */
	d->transition_p = realloc(d->transition_p, (n_total_trans + 1) * 
			sizeof(*d->transition_p));
	d->exp_trans_p = calloc((n_total_trans + n_uniq_kmers + 1), 
			sizeof(*d->exp_trans_p));

	memcpy(d->exp_trans_p, d->transition_p, (n_total_trans + 1) *
			sizeof(*d->transition_p));
	memcpy(d->exp_trans_p + n_total_trans + 1,
			uniq_kmer_cnts, n_uniq_kmers * sizeof(*uniq_kmer_cnts));

	d->n_total_trans = n_total_trans + n_uniq_kmers + 1;

	/* allocate memory */
	d->nonuniq_kmers = calloc(n_nonuniq_kmers, sizeof(*d->nonuniq_kmers));
	d->avg_initial_p = calloc(n_nonuniq_kmers + n_uniq_kmers + 1,
			sizeof(*d->avg_initial_p));

	d->n_nonzero_kmers = n_nonuniq_kmers + n_uniq_kmers;

	const int _uoff = n_nonuniq_kmers - n_total_trans - 1;
	d->uniq_kmer_init_off = _uoff;

	int *kcnt_histo = calloc(1001, sizeof(*kcnt_histo));
	int *kcnt_ss_histo = calloc(1001, sizeof(*kcnt_histo));

	size_t _poff = 1, _kidx = 0;
	_uidx = n_total_trans + 1;
	de = d->kmer_dict->table;
	for (int i = 0; i < n_kmers; de++) {
		if (de->value == NULL) continue;

		kmer_t *km = (kmer_t *) &de->value;

		if (km->n_trans == 1) {
			km->index = _uidx++;
			d->avg_initial_p[_uidx + _uoff] = 
				d->exp_trans_p[km->index] / total_trans_cnt;
		}
		else if (km->n_trans > 1) {
			double sum_trans_cnt = 0.0;
			for (int j = 0; j < km->n_trans; j++) {
				sum_trans_cnt += d->transition_p[_poff + j];
			}

			for (int j = 0; j < km->n_trans; j++) {
				d->transition_p[_poff + j] = LOG(
						d->transition_p[_poff + j] / sum_trans_cnt);
			}

			km->index = _kidx;
			d->nonuniq_kmers[_kidx++].poff = _poff;
			d->avg_initial_p[_kidx] = sum_trans_cnt / total_trans_cnt;

			_poff += km->n_trans;
		}
		else {
			km->index = 0;
		}

		i++;
	}

	de = d->kmer_dict->table;
	for (int i = 0; i < n_kmers; de++) {
		if (de->key == NULL) continue;

		kbits_t kid = (kbits_t) de->key - 1;
		kmer_t *km = (kmer_t *) &de->value;

		if (km->dirty == 1) {
			kbits_t rc_kid = kmer_reverse_complement(kid, kmer_len);
			kmer_t *rc_km = (kmer_t *) &(dict_find(d->kmer_dict, 
						kbits_cast_to_ptr(rc_kid))->value);

			km->dirty = 0;
			rc_km->dirty = 0;

			if (rc_km->coverage == 0 || km->coverage == 0) {
				km->single_stranded = 1;
				rc_km->single_stranded = 1;
			}

			/* using prev_n_trans if dirty or n_trans o/w:
			 * return km->index + _uoff if no. trans is 1 */
			int _idx = kmer_initial_index(km, _uoff),
				_ridx = kmer_initial_index(rc_km, _uoff);

			/*
			double mean_init_p = LOG(0.5 * (
						d->avg_initial_p[_idx] + d->avg_initial_p[_ridx]));
			*/
			/* joint coverage on both strands: record in histogram */
			int canon_kcnt = (int) round(total_trans_cnt * (d->avg_initial_p[_idx] +
						d->avg_initial_p[_ridx]));
			double mean_init_p = LOG(canon_kcnt);

			++kcnt_histo[canon_kcnt > 1000 ? 1000 : canon_kcnt];
			if (km->single_stranded) 
				++kcnt_ss_histo[canon_kcnt > 1000 ? 1000 : canon_kcnt];

			/*
			printf("%lu[L%uT%u] %lu[L%uT%u] idx: %u iidx: %d %d %lg %lg %lg\n",
					kid, km->label, km->n_trans,
					rc_kid, rc_km->label, rc_km->n_trans,
					km->index, _idx, _ridx, 
					d->avg_initial_p[_idx],
					d->avg_initial_p[_ridx],
					mean_init_p);
			*/

			/*
			if (IS_ZERO(mean_init_p)) {
				mean_init_p = 0;
			}
			*/

			if (_idx > 0) 
				d->avg_initial_p[_idx] = mean_init_p;
			if (_ridx > 0)
				d->avg_initial_p[_ridx] = mean_init_p;
		}

		i++;
	}

	d->avg_initial_p[0] = NAN;

	printf("Kmer multiplicity histogram: ");
	for (int i = 0; i < 100; ++i) {
		printf("%d ", kcnt_histo[i]);
	}
	printf("\n");

	printf("Single-stranded kmer multiplicity histogram: ");
	for (int i = 0; i < 100; ++i) {
		printf("%d ", kcnt_ss_histo[i]);
	}
	printf("\n");

	double _threshold = 0.0;
	if (d->opt->auto_penalty_eta) {
		d->opt->penalty_eta = determine_threshold(kcnt_histo,
				d->opt->penalty_gamma, halve_threshold, &_threshold);

		INFO_TEE(d->log_fp, "Parameter \\eta updated to: %lg\n",
				d->opt->penalty_eta);
	}

	free(kcnt_histo);
	free(kcnt_ss_histo);

	INFO_TEE(d->log_fp, "# Unique kmers: %zu\n", n_uniq_kmers);

	memset(d->exp_trans_p, 0, d->n_total_trans * sizeof(double));

	/* FIX: unless PMR_EXTEND_READS do not need this allocation, but
	 * must also fix weigh_observed_transitions */
	double *_exp_trans_p = calloc(d->n_total_trans, sizeof(double));
	// weight transitions by quality scores
	io_iterate_reads_fastq(d, d->opt->n_threads, d->opt->is_tty, 
			weigh_observed_transitions, NULL, _exp_trans_p);


#if PMR_EXTEND_READS
	// re-compute transition probabilities using the weighted observed counts.
	de = d->kmer_dict->table;
	for (int i = 0; i < n_kmers; de++) {
		if (de->key == NULL) continue;
		kmer_t *pk = (kmer_t *) &de->value;

		kbits_t kid = kbits_id_of(pk);

		double *k_w_cnts = _exp_trans_p + 
			(pk->unique_trans ? pk->index : d->nonuniq_kmers[pk->index].poff);
		double *k_obs_cnts = d->exp_trans_p + 
			(pk->unique_trans ? pk->index : d->nonuniq_kmers[pk->index].poff);

		double sum_cnt = 0.0;
		for (int b = 0; b < pk->n_trans; ++b) {
			sum_cnt += k_obs_cnts[b];
		}
		for (int b = 0; b < pk->n_trans; ++b) {
			k_obs_cnts[b] = LOG(k_w_cnts[b] / sum_cnt);
		}
		
		i++;
	}

	mstep_extend_all_kmers(d);
#endif
	free(d->exp_trans_p);
	d->exp_trans_p = _exp_trans_p;

	double _eta = d->opt->penalty_eta;
	d->opt->penalty_eta = _eta / 2;

	// apply the penalty to sparsify the transition parameters.
	//if (d->opt->mode == PMR_MODE_ORACLE)
		//cmstep_penalized_transition(d, 0, 0); // does not impose penalty, to preserve transition flags.
	//else
	if (d->opt->mode == PMR_MODE_ERROR_CORRECTION)
		cmstep_penalized_transition(d, d->opt->init_penalty_func, 0);

	d->opt->penalty_eta = _eta;
}

#if 0
int compute_kmer_hamming_nbhds(data *d)
{
	int dmax = d->opt->preconstructed_nbhd_hd;
	int num_perms = d->n_choose_k[d->opt->kmer_length][dmax];
	int n_kmers = d->kmer_dict->used;

	kbits_t *first_kmer_ids = malloc(n_kmers * sizeof(kbits_t));

	dictentry *de;
	register int i, n = 0;
	for (i = n_kmers, de = d->kmer_dict->table; i > 0; de++) {
		if (de->value == NULL) continue;
		i--;

		kmer_t *pk = (kmer_t *) de->value;
		pk->n_neighbors = 1;
		pk->alloc_nbhd_size = 2;
		pk->data = malloc(sizeof(kmer_t *) << 2);
		/* a kmer is always a neighbor of itself */
		*(kmer_t **) pk->data = pk;

		kmer_ids[n++] = pk->id;
	}

	mask_t pos_mask = {0, 0};
	for (i = 0; i < num_perms; i++) {
		INFO("Neighborhood construction iteration %2d/%2d...\n", i+1, num_perms);

		int shift = (i << 1);
		kbits_t pmask = d->position_masks[i];
		pos_mask.mask_shift = shift;
		pos_mask.pmask = pmask;

		/* sort masked kmer ids */
		qsort_r(kmer_ids, n_kmers, sizeof(kbits_t), 
				masked_kmer_id_comp, (void *) &pos_mask);

		kbits_t *base_loc = kmer_ids;

		/* identify neighborhood for all first kmers */
		for (int j = 0; j < n_kmers; j++) {
			kmer_t *pk = dict_find(d->kmer_dict, kbits_cast_to_ptr(kmer_ids[j]))->value;

			kbits_t kid = pk->id;
			register kbits_t masked_kid = pk->id & pmask;

			/* find the location of the kmer id in sorted list */
			kbits_t *id_loc = base_loc;
			for (; *id_loc != kid; id_loc++);
			base_loc = id_loc;
		
			/* identify the neighborhood by exploring kmers within locality */
			kbits_t *ptr_lb, *ptr_ub;
			kbits_t masked_base = (~pmask & kid) >> shift;
			for (ptr_lb = id_loc - masked_base; ptr_lb <= id_loc + (3UL - masked_base) && 
					masked_kid != (pmask & *(ptr_lb)); ptr_lb++);
			for (ptr_ub = ptr_lb + 1; masked_kid == (pmask & *(ptr_ub)); 
					ptr_ub++);

			int n_masked_nbhd_size = ptr_ub - ptr_lb - 1;
			if (n_masked_nbhd_size > 0) {
				int pnn = pk->n_neighbors;
				pk->n_neighbors += n_masked_nbhd_size;

				int size_nbhd = pk->n_neighbors;
				int realloc_size = 1 << pk->alloc_nbhd_size;
				if (realloc_size < size_nbhd) {
					for (; realloc_size < size_nbhd; realloc_size <<= 1, 
							pk->alloc_nbhd_size++);
					pk->data = realloc(pk->data, 
							sizeof(kmer_t *) << pk->alloc_nbhd_size);
				}

				kmer_t **nb_ptrs = (kmer_t **) pk->data;
				for (kbits_t *pid = ptr_lb; pid < ptr_ub; pid++) {
					register kbits_t _id = *pid;
					if (_id != kid) {
						nb_ptrs[pnn++] = dict_find(d->kmer_dict, 
								kbits_cast_to_ptr(_id))->value;
					}
				}
			}
		}
	}

	free(kmer_ids);
}
#endif

void kmer_id2seq(kbits_t kid, int size)
{
	int i;
	for (i = 0; i < size; i++){
		switch((kid >> (i * 2)) & 3){
			case 0:
				putchar('A');
				break;
			case 1:
				putchar('C');
				break;
			case 2:
				putchar('T');
				break;
			case 3:
				putchar('G');
				break;
		}
	}
}

/* --- precomputed tables --- */

/* Table that converts an (offset) Phred quality score to the corresponding 
 * error probability of base calling. 
 * For display, the quality score is usually added to an offset (33 by
 * default).
 */
int compute_phred_to_probability(data *d)
{
	int i; 
	d->phred_to_prob = malloc(sizeof(*(d->phred_to_prob)) * 
			(d->opt->max_qual_score + 1));

	if (d->phred_to_prob == NULL)
		return message(stderr, file_name, __LINE__, ERROR_MSG,
				MEMORY_ALLOCATION, NULL);

	for (i = 0; i <= d->opt->max_qual_score; i++)
			d->phred_to_prob[i] = POW(10., -((double)i/10));

	return NO_ERROR;
}	/* compute_phred_to_probablity */

/** 
 * Compute binomial coefficients, i.e. N choose K using Pascal triangle.
 * The largest N to compute is kmer_length + 1.
 */
int compute_binomial_coefficients(data *d)
{
	int n = d->opt->kmer_length + 1;
	int i, j;

	d->n_choose_k = malloc(n * sizeof *d->n_choose_k);
	if (d->n_choose_k == NULL)
		return message(stderr, file_name, __LINE__,
				ERROR_MSG, MEMORY_ALLOCATION, NULL);

	for (i = 0; i < n; i++){
		d->n_choose_k[i] = calloc(i+1, sizeof **d->n_choose_k);
		if (d->n_choose_k[i] == NULL)
			goto ERROR_RETURN;

		/* Pascal's triangle */
		for (j = 0; j <= i; j++){
			if (j == 0 || j == i)
				d->n_choose_k[i][j] = 1;
			else
				d->n_choose_k[i][j] = d->n_choose_k[i-1][j] +
					d->n_choose_k[i-1][j-1];
		}
	}

	return NO_ERROR;

ERROR_RETURN:
	for (n = 1; n < i; n++)
		free(d->n_choose_k[n]);
	free(d->n_choose_k);
	return message(stderr, file_name, __LINE__, ERROR_MSG,
			MEMORY_ALLOCATION, NULL);
}	/* compute_binomial_coefficients */

/* Permute all (k choose d) combinations for d erroneous positions in any kmer,
 * d = 1, 2, ..., d_max.
 * 
 * All computed results are stored in two-dimenional array: d->position_perms.
 * The 1st dimension of this array corresponds to d = 1, 2, ..., d_max.
 * The 2nd dimension for given d, is of length [(k choose d) * d]. 
 * Every d consecutive positions store a combination of d positions.
 */
int compute_error_positions(data *d)
{
	int ne, i, j;	/* number of errors */
	int num_perms;
	char *ptr_slots;	/* pointer to current d slots */
	int kmer_len = d->opt->kmer_length;
	int dmax = d->opt->preconstructed_nbhd_hd;
	kbits_t mask;

	d->position_perms = malloc(dmax * sizeof *d->position_perms);
	d->position_masks = malloc(d->n_choose_k[kmer_len][dmax] * 
			sizeof *d->position_masks);

	if (d->position_perms == NULL || d->position_masks == NULL)
		return message(stderr, file_name, __LINE__, ERROR_MSG,
				MEMORY_ALLOCATION, NULL);

	for (ne = 1; ne <= dmax; ne++){
		num_perms = d->n_choose_k[kmer_len][ne];
		d->position_perms[ne-1] = malloc(num_perms * ne * 
				sizeof **d->position_perms);

		if (d->position_perms[ne-1] == NULL)
			goto ERROR_RETURN;

		ptr_slots = d->position_perms[ne-1];

		/* initialization */
		mask = KBITS_MASK_MAX;
		for (i = 0; i < ne; i++){
			ptr_slots[i] = i + 1;
			mask &= ~(3UL << (i << 1));
		}

		if (ne == dmax)
			d->position_masks[0] = mask;

		memcpy(ptr_slots + ne, ptr_slots, ne * sizeof(*ptr_slots));
		ptr_slots += ne;	/* move pointer to next permutation */

		/* In total we have (n choose k) permutations, among which the initial
		 * permutation (1, 2, ..., d) counts 1. 
		 * Following code permute the remaining (n choose k) - 1 combinations.
		 */
		for (i = 1; i < num_perms; i++){
			ptr_slots[ne-1]++;	/* increase the rightmost position by 1 */
			for (j = ne - 1; j >= 0; j--){
				if (ptr_slots[j] > kmer_len - (ne - 1 - j)){
					/* carray on */
					ptr_slots[j-1]++;
					ptr_slots[j]  = ptr_slots[j-1] + 1;
					/* check if we can reset the next radix */
					if (j < ne - 1 && ptr_slots[j+1] > (kmer_len-(ne-2-j))){
						ptr_slots[j+1] = ptr_slots[j] + 1;
					}
				}
				else
					break;
			}

			if (ne == dmax){
				mask = KBITS_MASK_MAX; 
				for (j = 0; j < ne; j++){
					mask &= ~(3UL << ((ptr_slots[j] - 1) << 1));
				}
				d->position_masks[i] = mask;
			}

			if (ptr_slots[0] > (kmer_len - ne + 1))
				break;

			if (i < num_perms - 1){
				memcpy(ptr_slots + ne, ptr_slots, ne * sizeof(*ptr_slots));
				ptr_slots += ne;
			}
		}
	}

	return NO_ERROR;

ERROR_RETURN:
	for (i = 0; i < ne - 1; i++)
		free(d->position_perms[i]);
	free(d->position_perms);

	return message(stderr, file_name, __LINE__, ERROR_MSG,
			MEMORY_ALLOCATION, NULL);
}	/* compute_error_positions */

/* compute all (k choose d) x 4^d error patterns, 
 * for d = 1, 2, ..., dmax.
 *
 */
int compute_all_errors(data *d)
{
	int i, j, k;
	int offset, num_perms, num_err_combs, num_states;
	kbits_t error_mask, dmer_mask, pos_shift;
	int kmer_len = d->opt->kmer_length;
	int dmax = d->opt->preconstructed_nbhd_hd;

	num_perms = d->n_choose_k[kmer_len][dmax];
	num_err_combs = pow(4, dmax);
	num_states = num_err_combs * num_perms;

	d->error_perms = malloc(sizeof(*d->error_perms) * num_states);
	if (d->error_perms == NULL)
		return MEMORY_ALLOCATION;

	for (i = 0; i < num_perms; i++){
		offset = i * dmax;
		for (j = 0; j < num_err_combs; j++){
			dmer_mask = j;
			error_mask = 0;

			for (k = 0; k < dmax; k++){
				pos_shift = (d->position_perms[dmax-1][offset+k] - 1) << 1;
				error_mask |= ((dmer_mask & 3) << pos_shift);

				dmer_mask >>= 2;
			}

			d->error_perms[i*num_err_combs+j] = error_mask;
		}
	}

	return NO_ERROR;
}
