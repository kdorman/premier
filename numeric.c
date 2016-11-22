/**
 * @file numeric.c
 * Mathematical functions, particularly for implementing log calculations.
 */
#include "numeric.h"

//static const char *file_name = "numeric.c";

typedef union {
	DOUBLE dbl;
	uint64_t u64;
} dbl_t;

DOUBLE log_product(DOUBLE l1, DOUBLE l2) {
	if (isnan(l1) || isnan(l2))
		return NAN;
	else
		return l1+l2;
} /* log_product */

/* atomic log add, which equals to,
 *		l1 = SUM(l1, l2).
 * this function utilizes atomic compare_and_swap builtin to update l1 with
 * the sum in a lock-free fashion.
 *
 * @return (l1+l2)
 */
inline DOUBLE log_sum(int n, DOUBLE *summands, int argmax, DOUBLE maxval)
{
	register DOUBLE _expsum = 0.0;
	for (int i = 0; i < n; i++) {
		uint64_t argmax_mask = (uint64_t) (i == argmax) - 1UL;
		dbl_t summand = {.dbl = EXP(summands[i] - maxval)};
		summand.u64 &= argmax_mask;

		_expsum += summand.dbl;
	}
	return PRODUCT(maxval, log1p(_expsum));
}

DOUBLE log_add(DOUBLE l1, DOUBLE l2) {
	if (isnan(l1) || isnan(l2)) {
		if (isnan(l1))
			return l2;
		else
			return l1;
	} else {
		if (l1 > l2)
			return l1 + log1p(REAL_EXP(l2-l1));
		else
			return l2 + log1p(REAL_EXP(l1-l2));
	}
} /* log_sum */

/**
 * Overflow-safe log.
 * @param d value to be logged
 * @param logged value or NAN
 */
DOUBLE elog(DOUBLE d)
{
	if (d == 0)
		return NAN;
	else if (d>0)
		return REAL_LOG(d);
	else {
		//message(stderr, file_name, "elog", __LINE__, ERROR_MSG,
	//		CUSTOM_ERROR, "invalid negative argument\n");
		return NAN;
	}
} /* elog */

DOUBLE eexp(DOUBLE d) {
	if (isnan(d))
		return 0.0;
	else 
		return REAL_EXP(d);
} /* eexp */

DOUBLE emax(DOUBLE a, DOUBLE b) {
	if (isnan(a)) {
		/* really doesn't matter if b is NAN or not */
		return b;
	}
	else {
		return (isnan(b)) ? a : ((a < b) ? b : a);
	}
}

