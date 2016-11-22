#include "numeric.h"

#include <stdlib.h>

/* TODO: wrap this function into a type generic macro */

/* TODO: optimize memory usage, avoid using calloc() and free() frequently.
 * The optimization may be based on the fact that |Al| < |A|, |As| < |A|.
 */

/* QuickSelect(A, k)
 *   Finds the k-th largest element within arrary A of length n.
 */
double quick_select_dbl(double *A, int n, int k)
{
	int iter = 0;
	do {
		int r = (int) (((double) rand() / (RAND_MAX)) * n);
		double pivot = A[r];

		double *Al = calloc(n, sizeof(*Al)); /* larger elements */
		double *As = calloc(n, sizeof(*As)); /* smaller elements */
		register int nl = 0, ns = 0;

		for (int i = 0; i < n; i++) {
			register double elem = A[i];
			if (isnan(pivot)) {
				if (!isnan(elem))
					Al[nl++] = elem;
			}
			else {
				if (isnan(elem) || elem < pivot)
					As[ns++] = elem;
				else if (elem > pivot)
					Al[nl++] = elem;
			}
		}

		if (iter++ > 0) free(A);

		if (k <= nl) {
			/* kth largest element is within pile Al */
			A = Al;
			n = nl;
			free(As);
		}
		else if (k > n - ns) {
			/* kth largest element is within pile As */
			A = As;
			k = k - (n - ns);
			n = ns;
			free(Al);
		}
		else {
			free(Al);
			free(As);
			return pivot;
		}
	} while (1);
}
