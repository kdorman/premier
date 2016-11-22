#ifndef __PMR_ATOMIC_H__
#define __PMR_ATOMIC_H__

#include <stdint.h>

/* atomic addition using GCC's atomic compare_and_swap intrinsics */

/* atomic version of expression:
 *
 *		*ptr = *ptr + val 
 * */
static inline double atomic_add_dbl(double *ptr, double val)
{
	double oldval;
	double newval;
	do {
		oldval = *ptr;
		newval = oldval + val;
	} while(!__sync_bool_compare_and_swap((uint64_t *) ptr, 
				*(uint64_t *) &oldval, *(uint64_t *) &newval));

	return newval;
}

static inline void atomic_cas_dbl(double *ptr, double newval)
{
	double oldval;
	do {
		oldval = *ptr;
	} while(!__sync_bool_compare_and_swap((uint64_t *) ptr, 
				*(uint64_t *) &oldval, *(uint64_t *) &newval));
}

static inline uint64_t atomic_add_u64(uint64_t *ptr, uint64_t val)
{
	uint64_t oldval;
	uint64_t newval;
	do {
		oldval = *ptr;
		newval = oldval + val;
	} while(!__sync_bool_compare_and_swap(ptr, oldval, newval));

	return newval;
}

static inline uint32_t atomic_add_u32(uint32_t *ptr, uint32_t val)
{
	uint32_t oldval;
	uint32_t newval;
	do {
		oldval = *ptr;
		newval = oldval + val;
	} while(!__sync_bool_compare_and_swap(ptr, oldval, newval));

	return newval;
}

#endif
