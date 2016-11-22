#ifndef __QTRF_SORT_H__
#define __QTRF_SORT_H__

typedef struct tuple_char_s {
	int order;
	char value;
} tuple_char_t;

typedef struct tuple_dbl_s {
	int64_t order;
	double value;
} tuple_dbl_t;

static int cmp_tuple_char(const void *a, const void *b) {
	const tuple_char_t *ta = (const tuple_char_t *) a;
	const tuple_char_t *tb = (const tuple_char_t *) b;

	return ((int) ta->value - tb->value);
}

static int cmp_tuple_dbl(const void *a, const void *b) {
	const tuple_dbl_t *ta = (const tuple_dbl_t *) a;
	const tuple_dbl_t *tb = (const tuple_dbl_t *) b;

	return (ta->value < tb->value) ? 1 : -1;
}

static int cmp_double(const void *a, const void *b) {
	return (*(double *) a < *(double *) b) ? 1 : -1;
}

static inline int *sort_order_char(char *y, int n, char **sorted_y)
{
	tuple_char_t *y_enum = (tuple_char_t *) malloc(n * sizeof(*y_enum));
	for (int i = 0; i < n; i++) {
		y_enum[i].value = y[i];
		y_enum[i].order = i;
	}

	qsort(y_enum, n, sizeof(*y_enum), cmp_tuple_char);
	char *srty = NULL;

	if (sorted_y != NULL) {
		srty = (char *) malloc(n * sizeof(*y));
	}

	int *sord = (int *) malloc(n * sizeof(int));

	for (int i = 0; i < n; i++) {
		if (srty) srty[i] = y_enum[i].value;
		sord[i] = y_enum[i].order;
	}

	free(y_enum);

	if (sorted_y != NULL) {
		*sorted_y = srty;
	}

	return sord;
}

static inline int64_t *sort_order(double *y, int n, double **sorted_y)
{
	tuple_dbl_t *y_enum = (tuple_dbl_t *) malloc(n * sizeof(*y_enum));
	for (int i = 0; i < n; i++) {
		y_enum[i].value = y[i];
		y_enum[i].order = i;
	}

	qsort(y_enum, n, sizeof(*y_enum), cmp_tuple_dbl);
	double *srty = y;

	if (sorted_y != NULL) {
		srty = (double *) malloc(n * sizeof(*y));
	}

	int64_t *sord = (int64_t *) malloc(n * sizeof(int64_t));

	for (int i = 0; i < n; i++) {
		srty[i] = y_enum[i].value;
		sord[i] = y_enum[i].order;
	}

	free(y_enum);

	if (sorted_y != NULL) {
		*sorted_y = srty;
	}

	return sord;
}

#if 0
static inline int *sort_rank(int *y_ord, int n)
{
	tuple_int_t *o_enum = malloc(n * sizeof(*o_enum));
	for (int i = 0; i < n; i++) {
		o_enum[i].order = y_ord[i];
		o_enum[i].rank = i;
	}

	qsort(o_enum, n, sizeof(*o_enum), cmp_tuple_int);
	int *srnk = malloc(n * sizeof(int));

	for (int i = 0; i < n; i++) {
		srnk[i] = o_enum[i].rank;
	}

	free(o_enum);

	return srnk;
}
#endif

#endif
