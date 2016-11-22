#ifndef __PREMIER_IODATA_H__
#define __PREMIER_IODATA_H__

#include <stdio.h>
#include <time.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <sys/time.h>
#include <unistd.h>
#include <ctype.h>
#include <math.h>

#include "premier.h"

#define SKIP_TO_NEWLINE(p)      \
	for(; *(p) != '\n'; (p)++); \
	(p)++;                      \

#define IO_LOCATE_NEWLINE                                                      \
	for(pnr = p;  pnr < buf_max && *pnr != '\n'; pnr++);                   \
	if (pnr == buf_max || *(pnr) != '\n') {                                \
		partial_read_size = (pnr - pread + 1);                         \
		break;                                                         \
	}                                                                      \

/*
int setup_worker_entries(data *d);
*/

char *io_cat_output_filename(options *opt, const char *suffix, 
		char *ext);

void io_draw_progress_bar(int i, int n, double elapsed, int calc_eta, int finished,
		const char *unit);

/*
int io_mmap_fastq_file(data *d, const char *filename);
void io_munmap_fastq_file(data *d);
*/
int io_iterate_reads_fastq(data *d, int nthreads,
		int progress_bar,
		void (*mapper)(data *d, read_t *read, void *fdata), 
		void (*reducer)(data *d, void *fdata), 
		void *fdata);
		
#endif
