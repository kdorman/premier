#include "kmer.h"
#include "data.h"
#include "mempool.h"
#include "kmerdict.h"
#include "atomic.h"

void kmer_incr_count_by_id(data *d, kbits_t kid, kbits_t next_base, 
		uint64_t flag)
{
	/* dictionary value:
	 *
	 * GGGGGGGGGGGGGG TTTTTTTTTTTTTT CCCCCCCCCCCCCC AAAAAAAAAAAAAA RCTF FLAG
	 * |____________| |____________| |____________| |____________| |__| |__| 
	 *       ^                                                      ^
	 *  number of transitions into G base                           |
	 *                                                    transitions seen
	 *                                                    at rev. comp. strand
	 * FLAG indicates first kmer, normal kmer, first observed on reverse complement
	 * RCTF = Reverse Complement Transition Flag
	 */

	static uint64_t base_masks[5] = {
		((1UL << 14) - 1) << 8, 
		((1UL << 14) - 1) << (8 + 14), 
		((1UL << 14) - 1) << (8 + 28), 
		((1UL << 14) - 1) << (8 + 42),
		0UL
	};

	const uint64_t one_rep4 = (1UL << 8) | (1UL << 22) | (1UL << 36) | 
		(1UL << 50);

	dictentry *de;
	de = dict_addkey_find(d->kmer_dict, kbits_cast_to_ptr(kid));

	uint64_t *kcnter = (uint64_t *) &(de->value);

	/* if FLAG_SET_RC is set and not PMR_USE_JOINT_INITIALIZATION,
	 * do not add the counter; set the corresponding RCTF instead.
	 */

	int update_total_trans = 0;
	uint64_t _kold, _ksnap;
	do {
		_kold = *kcnter;
		_ksnap = _kold;

		uint64_t rc_label = flag & FLAG_SET_RC;

#if PMR_USE_JOINT_INITIALIZATION 
		uint64_t _mask = base_masks[next_base];
#else
		uint64_t rc_mask = (rc_label >> 3) - 1UL;
		uint64_t _mask = base_masks[(next_base & rc_mask) + (~rc_mask & 4)];
#endif

		/* assign label only when the kmer is seen for the first time */
		if (_ksnap == 0) {
			_ksnap |= rc_label; 
		}

		/* cap trans. flag at 1111 */
#if PMR_USE_JOINT_INITIALIZATION 
		_ksnap |= ((15 & (1UL << next_base)) << 4);
#else
		_ksnap |= (~rc_mask & ((15 & (1UL << next_base)) << 4));
#endif

		uint64_t tcnt = (_ksnap & _mask);
		/* count new transition to non-N base 1st time encountered */
		if (tcnt == 0 && _mask > 0) {
			update_total_trans = 1;
			//d->n_total_trans += 1;
		}

		_ksnap |= (flag & 0x7);
		/* guard overflow */
		if (tcnt < _mask) {
			_ksnap += (one_rep4 & _mask);
		}

	/* to avoid race condition */
	} while(!__sync_bool_compare_and_swap(kcnter, _kold, _ksnap));

	if (update_total_trans) 
		atomic_add_u64(&d->n_total_trans, 1);

	//*kcnter = _ksnap;
}
