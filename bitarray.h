#ifndef __BITARRAY_H__
#define __BITARRAY_H__

#ifdef __POPCNT__
#include <nmmintrin.h>
#endif

#include <stdint.h>
#include <stdlib.h>

#define BITS_TO_U64(nb) (((nb) >> 6) + (((nb) & ((1 << 6) - 1)) > 0))

static inline void bitarray_reset(uint64_t *ba, int nbits) 
{
	memset(ba, 0, BITS_TO_U64(nbits) << 3);
}

static inline void bitarray_set_pwr2(uint64_t *ba, int off_unit, uint64_t val, 
		int pwr)
{
	register int off = off_unit << pwr; 
	register int u64_off = off >> 6;
	register int shift = off & 63;

	register uint64_t len_mask = (1UL << (1 << pwr)) - 1;
	register uint64_t capped_val = val & len_mask; 

	ba[u64_off] &= ~(len_mask << shift);
	ba[u64_off] |= (capped_val << shift);
}

static inline uint64_t bitarray_get_pwr2(uint64_t *ba, int off_unit, int pwr)
{
	register int len = (1 << pwr);
	register int off = off_unit << pwr; 
	register int u64_off = off >> 6;
	register int shift = off & 63;

	register uint64_t masked_bits = (ba[u64_off] >> shift) & ((1UL << len) - 1);

	return masked_bits;
}

static inline void bitarray_set(uint64_t *ba, int off, uint64_t val, int len)
{
	register int u64_loff = off >> 6;
	register int u64_uoff = (off + len - 1) >> 6;

	register uint64_t len_mask = (1UL << len) - 1;
	register int capped_val = val & len_mask; 
	register int lshift = off & 63;
	register int ushift = (off + len) & 63;

	register uint64_t u64_lmask = ~(len_mask << lshift);
	register uint64_t u64_umask = ~(len_mask >> ushift);

	ba[u64_loff] = (ba[u64_loff] & u64_lmask) | (capped_val << lshift);
	ba[u64_uoff] = (ba[u64_uoff] & u64_umask) | (capped_val >> ushift);
}

static inline uint64_t bitarray_get(uint64_t *ba, int off, int len)
{
	register int u64_loff = off >> 6;
	register int u64_uoff = (off + len - 1) >> 6;

	register int shift = off & 63;

	register uint64_t masked_bits = ((ba[u64_loff] >> shift) | 
		(ba[u64_uoff] << (64 - shift))) & ((1UL << len) - 1);

	return masked_bits;
}

static inline void bitarray_set1(uint64_t *ba, int off)
{
	register int u64_offset = off >> 6;
	ba[u64_offset] |= 1UL << (off & 63);
}

static inline uint64_t bitarray_get1(uint64_t *ba, int off)
{
	register int u64_offset = off >> 6;
	register int shift = off & 63;
	return (ba[u64_offset] >> shift) & 1UL;
}

static inline int bitarray_countn(uint64_t *ba, int off, int len)
{
#ifdef __POPCNT__
	register uint64_t u64 = bitarray_get(ba, off, len);
	register int n = _mm_popcnt_u64(u64);
	return n;
#else
	/* TODO: fallback code */
	return 0;
#endif
}

#endif
