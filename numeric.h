#ifndef __NUMERIC_H__
#define __NUMERIC_H__

#include <stdio.h>
#include <stdint.h>
#include <math.h>

#define USE_DOUBLE

#define USE_LOG_SCALE 1

#define MAX(X,Y) ( (X) > (Y) ? (X) : (Y) )
#define MIN(X,Y) ( (X) < (Y) ? (X) : (Y) )

#define update_max_quantity(val, idx, max_val, amax_val)                       \
	do {                                                                   \
		if (isnan(max_val) || (max_val < val)) {                       \
			max_val = val;                                         \
			amax_val = idx;                                        \
		}                                                              \
	} while (0);

typedef unsigned char boolean;
/* defined in Rmath.h: enum {FALSE, TRUE}; */

#define DOUBLE double
#define DOUBLE_F_FMT "%lf"
#define DOUBLE_F_NFMT "lf"
#define DOUBLE_G_FMT "%lg"
#define DOUBLE_G_NFMT "lg"
#define DOUBLE_E_FMT "%le"
#define DOUBLE_E_NFMT "le"
#define FABS fabs
#define POW pow
#ifndef USE_LOG_SCALE
#	define LOG(a)          log((a))
#	define EXP(a)          exp((a))
#endif

#ifndef isnan
# define isnan(x) \
	(sizeof (x) == sizeof (long double) ? isnan_ld (x) \
	: sizeof (x) == sizeof (double) ? isnan_d (x) \
	: isnan_f (x))
static inline int isnan_f  (float       x) { return x != x; }
static inline int isnan_d  (double      x) { return x != x; }
static inline int isnan_ld (long double x) { return x != x; }
#endif

#if USE_LOG_SCALE
/* we don't really need the function log_product because the if condition is 
 * redundant --- if either a or b is NAN, a+b is guranteed to be NAN by
 * IEEE 754 standard (is this claim true?) */
#	define PRODUCT(a,b)    ((a)+(b))
#	define SUM(a,b)        log_add((a),(b))
#	define LOG(a)          elog((a))
#	define EXP(a)          eexp((a))
#	define REAL_LOG(a)     log((a))
#	define REAL_EXP(a)     exp((a))
#	define IS_ZERO(a)      isnan((a))
#else
#	define PRODUCT(a,b)    ((a)*(b))
#	define SUM(a,b)        ((a)+(b))
#	define REAL_LOG(a)     LOG((a))
#	define REAL_EXP(a)     EXP((a))
#	define IS_ZERO(a)      ((a) == 0)
#endif

/** Various mathematical constants expressed with
 * long double precision
 */
static const DOUBLE PI = 3.141592653589793238;
static const DOUBLE LOG_PI = 1.144729885849400135;
static const DOUBLE SQRT_2PI = 2.506628274631000453;
static const DOUBLE LOG_SQRT_2PI = 0.918938533204672722;
static const DOUBLE LOG_HALF = -0.693147180559945309;

DOUBLE log_product(DOUBLE l1, DOUBLE l2);
DOUBLE log_add(DOUBLE l1, DOUBLE l2);
DOUBLE log_sum(int n, DOUBLE *summands, int argmax, DOUBLE maxval);
DOUBLE elog(DOUBLE d);
DOUBLE eexp(DOUBLE d);

#endif
