#ifndef __PREMIER_TRANS_H__
#define __PREMIER_TRANS_H__

#include <stdint.h>

#include "premier.h"

extern const uint16_t trans_flag2packed_lookup_table[];
/* translate unique transition flag to corresponding index/offset in 
 * data->transition_p. */
extern const int trans_flag2index_lookup_table[];

#define trans_flag2packed(f) trans_flag2packed_lookup_table[(f)]
#define trans_flag2index(f) trans_flag2index_lookup_table[(f)]

/* count number of admissible transitions implied by the number of 
 * 1 bits in the supplied flag */
#define trans_flag_popcnt(f) (trans_flag2packed((f)) & 7)


/*
 * @param kid 		bit array of kmer in neighborhood of observed kmer
 * @param obs_kid	observed kmer bit array
 * @param obs_n_flag	indicate N base in obs_kid
 * @param k_trans_flag	transition flags of kid
 * @param obs_next_base	observed next base
 * @param hd		hamming distance of kid
 * @param dmax		maximum allowed hamming distance
 * @return uint64_t suffix_hamming_distance_above_4_bits | allowed_next_bases
 */
static inline ksubstr_t trans_admissible_transitions(kbits_t kid, 
		kbits_t obs_kid, kbits_t obs_n_flag, ksubstr_t k_trans_flag, 
		ksubstr_t obs_next_base, ksubstr_t hd, ksubstr_t dmax)
{
	/* do the kid and obs kid differ at first base? */
#ifdef PMR_USE_SCALAR_KMER_TYPE
	ksubstr_t mut_flag = ((kid ^ obs_kid) | obs_n_flag) & 3UL;
#else
	ksubstr_t mut_flag = ((kid.u64[0] ^ obs.u64[0]) | obs_n_flag.u64[0]) & 3UL;
#endif
	/* hamming distance of k-1 suffix */
	kbits_t sfx_hd = hd - ((mut_flag & 1UL) | (mut_flag >> 1));

	/* 1111 (if hd already at max) or next base, which does not add to hd, in 4 bits */
	ksubstr_t adm_trans_mask = (((uint64_t) (sfx_hd >= dmax)) - 1UL) |
		(1UL << obs_next_base);

	return ((sfx_hd << 4) | (k_trans_flag & adm_trans_mask));
}

#endif
