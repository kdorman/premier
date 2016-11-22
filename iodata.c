#include <omp.h>

#include "iodata.h"
#include "read.h"
#include "kvec.h"

#define PMR_PROG_BAR_LEN 25
#define PMR_PROG_REPORT_INTERVAL 500

void io_draw_progress_bar(int i, int n, double elapsed, int calc_eta, int finished,
		const char *unit)
{
	const char *prog_unit = unit == NULL ? "Reads" : unit;

	i = (i > n) ? n : i;

	if (!calc_eta || finished) {
		fprintf(stderr, "\r%s proc'ed: %10d [", prog_unit, 
				finished ? n : i);

		if (!finished) {
			int progress = i * PMR_PROG_BAR_LEN / n;

			for (int j = 0; j <= progress; ++j) {
				fputc('.', stderr);
			}
			for (int j = progress + 1; j < PMR_PROG_BAR_LEN; ++j) {
				fputc(' ', stderr);
			}
		}
		else {
			for (int j = 0; j < PMR_PROG_BAR_LEN; ++j) {
				fputc('.', stderr);
			}
		}

		double perc = finished ? 100.0 : (double) i / n * 100.0;
		fprintf(stderr, "] %5.1f%% Elapsed: %02d:%02d:%06.3f", 
				perc,
				(int) elapsed / 3600, ((int) elapsed % 3600) / 60, 
				elapsed - ((int) elapsed / 60) * 60
			);

		if (finished) fputc('\n', stderr); 
	}
	else {
		// ETA, in seconds
		int eta = elapsed * ((double) n / i) - elapsed;
		int progress = i * PMR_PROG_BAR_LEN / n;

		fprintf(stderr, "\r%s proc'ed: %10d [", prog_unit, i);

		for (int j = 0; j <= progress; ++j) {
			fputc('.', stderr);
		}
		for (int j = progress + 1; j < PMR_PROG_BAR_LEN; ++j) {
			fputc(' ', stderr);
		}

		fprintf(stderr, "] %5.1f%% ETA: %02d:%02d:%02d ", 
				(double) i / n * 100.0,
				eta < 0 ? 99 : eta / 3600, 
				eta < 0 ? 99 : (eta % 3600) / 60, 
				eta < 0 ? 99 : eta % 60);
	}
}

char *io_cat_output_filename(options *opt, const char *suffix, 
		char *ext)
{
	time_t t;
	struct tm *lts;
	char *fname;
	char time_string[32];
	char eta_string[32];

	fname = calloc(PMR_LINE_BUFFER_SIZE, sizeof(char));
	if (fname == NULL) return NULL;

	/* generate time info */
	time(&t);
	lts = localtime(&t);

	if (ext == NULL) {
		ext = "txt";
	}

	if (opt->filename_append_timeinfo) {
		strftime(time_string, 31, "_%m%d_%H%M", lts);
	}
	else {
		time_string[0] = '\0';
	}

	if (!isnan(opt->penalty_eta)) {
		snprintf(eta_string, 31, "%0.f", opt->penalty_eta);
	}
	else {
		snprintf(eta_string, 31, "%s", "auto");
	}

	if (opt->filename_append_params) {
		snprintf(fname, PMR_LINE_BUFFER_SIZE - 1, 
				"%s%s_K%dD%d_E%sG%1.0eR%1.0e%s.%s", 
				opt->output_prefix, time_string, 
				opt->kmer_length, opt->max_hamming_dist, 
				eta_string, opt->penalty_gamma, opt->tol,
				suffix, ext);
	}
	else {
		snprintf(fname, PMR_LINE_BUFFER_SIZE - 1, 
				"%s%s%s.%s", 
				opt->output_prefix, time_string, 
				suffix, ext);
	}

	return fname;
}

/*
int io_mmap_fastq_file(data *d, const char *filename)
{
	char *p;
	int fd;
	struct stat fs;
	size_t fastq_size;

	if ((fd = open(filename, O_RDONLY)) == -1){
		return message(stderr, __FILE__, __LINE__, ERROR_MSG,
				FILE_NOT_FOUND, filename);
	}
	else {
		if (fstat(fd, &fs) != -1)
			fastq_size = fs.st_size;
		else
			return message(stderr, __FILE__, __LINE__, ERROR_MSG, 
				FILE_NOT_FOUND, filename);
	}

	if (!fastq_size)
		return message(stderr, __FILE__, __LINE__, ERROR_MSG,
				FILE_NOT_FOUND, filename);

	p = mmap(0, fastq_size, PROT_READ, MAP_SHARED, fd, 0);
	if (p == MAP_FAILED)
		return message(stderr, __FILE__, __LINE__, ERROR_MSG,
				MEMORY_ALLOCATION, NULL);

	d->fastq_mmap = p;
	d->fastq_size = fastq_size;
	d->fastq_fd = fd;

	return NO_ERROR;
}
*/

/*
void io_munmap_fastq_file(data *d)
{
	if (d->fastq_fd) {
		munmap(d->fastq_mmap, d->fastq_size);
		close(d->fastq_fd);
	}
}
*/

int io_iterate_reads_fastq(data *d, int nthreads,
		int progress_bar,
		void (*mapper)(data *d, read_t *read, void *fdata), 
		void (*reducer)(data *d, void *fdata), 
		void *fdata)
{
	int fd;
	struct stat fs;
	struct timeval tv;
	int64_t fastq_size;
	size_t size_ident;
	int read_idx = 0;
	int read_processed = 0;
	char *pnr = NULL;

	gettimeofday(&tv, NULL);

	if ((fd = open(d->opt->fastq_file, O_RDONLY)) == -1){
		return message(stderr, __FILE__, __LINE__, ERROR_MSG,
				FILE_NOT_FOUND, d->opt->fastq_file);
	}
	else {
		if (fstat(fd, &fs) != -1)
			fastq_size = (int64_t) fs.st_size;
		else 
			return message(stderr, __FILE__, __LINE__, ERROR_MSG, 
				FILE_NOT_FOUND, d->opt->fastq_file);
	}

	close(fd);
	
	/* create a temporary read */
	read_t *tmp_read = read_create();
	tmp_read->identifier = NULL; // buf_ident;

	int chunk = 0;
	int chunks_to_read = (int) ceil((double) fastq_size / PMR_READ_BUFFER_SIZE);

	// pad 16 bytes to avoid invalid read with SSE intrinsics (_mm_loadu_xx).
	char *_buffer = malloc(16 + (fastq_size < PMR_READ_BUFFER_SIZE ? fastq_size :
			PMR_READ_BUFFER_SIZE));
	char *buf_loc = _buffer;

	double wt_start = 0.0;

	/* process in blocks */
	FILE *fp = fopen(d->opt->fastq_file, "r");
	while (fastq_size > 0) {
		size_t chunk_size = fastq_size + (buf_loc - _buffer) < PMR_READ_BUFFER_SIZE ? 
			fastq_size : PMR_READ_BUFFER_SIZE - (buf_loc - _buffer);

		size_t buf_size = fastq_size + (buf_loc - _buffer) < PMR_READ_BUFFER_SIZE ? 
			fastq_size + (buf_loc - _buffer) :
			PMR_READ_BUFFER_SIZE;

		kvec_t(read_t) read_batch;
		kv_init(read_batch);
		kv_resize(read_t, read_batch, buf_size / (128 << 2));

		char *buf_max = _buffer + buf_size - 1;

		char *p = _buffer;
		++chunk;

		if (!progress_bar) { 
			INFO_TEE(d->log_fp, "Reading in FASTQ data in chunk [%2d/%2d] (%zu / REM: %zu)...\n", chunk, 
					chunks_to_read, chunk_size, fastq_size);
		}

		fread(buf_loc, chunk_size, 1, fp);

		size_t partial_read_size = 0;
		char *pread = p;

		/* locate the last read in buffer */
		//size_t rem_buf_size = buf_size;
		while ((p - _buffer) < buf_size) {
			pread = p;
			buf_loc = _buffer;

			IO_LOCATE_NEWLINE;
			/* process header */
			char *p_ident = p + 1;
			for (; *p_ident != ' ' && p_ident < pnr; p_ident++);

			size_ident = (size_t) (p_ident - p - 1);
			tmp_read->identifier = p + 1;

			/*
			strncpy(buf_ident, p + 1, size_ident);
			buf_ident[size_ident] = '\0';
			*/

			/* read sequence */
			p = pnr + 1;
			tmp_read->sequence = p;

			IO_LOCATE_NEWLINE;

			int read_length = pnr - p;

			tmp_read->length = read_length;
			/* skip quality header line */
			p = pnr + 1;
			if (*p == '+') {
				IO_LOCATE_NEWLINE;
			}

			/* WAS BUG: break from above block did not escape as needed */
			if (pnr == buf_max || *pnr != '\n')
				break;

			/* quality score */
			p = pnr + 1;
			tmp_read->qscore = p;

			IO_LOCATE_NEWLINE;

			tmp_read->id = read_idx++;

			// overwrite buffer raw data.
			tmp_read->identifier[size_ident] = '\0';
			kv_push(read_t, read_batch, *tmp_read);

			p = pnr + 1;
		}

		size_t batch_size = kv_size(read_batch);
		size_t batch_groups = batch_size / PMR_BATCH_GROUP_SIZE;
		if (batch_groups * PMR_BATCH_GROUP_SIZE < batch_size) ++batch_groups;

		double wt_end = 0.0;
		// reads processed by all threads in current batch
		int batch_processed = 0;

		//INFO_TEE(d->log_fp, "Processing reads in chunk %d...\n", chunk);
		if (wt_start == 0.0) wt_start = omp_get_wtime();

		omp_set_num_threads(nthreads);
		#pragma omp parallel private(wt_end)
		{ //omp

		for (int b = 0; b < batch_groups; ++b) {
			int imin = b * PMR_BATCH_GROUP_SIZE;
			int imax = (b+1) * PMR_BATCH_GROUP_SIZE > batch_size ? batch_size :
				(b+1) * PMR_BATCH_GROUP_SIZE;

			#pragma omp for schedule(dynamic) reduction(+:batch_processed)
			for (int i = imin; i < imax; ++i) {
				// invoke mapping callback in parallel.
				read_t *read = &kv_A(read_batch, i);
				read->batch_id = i - imin;
				mapper(d, read, fdata);

				++batch_processed;
			}

			// try resize, only performed at the master thread.
			#pragma omp master 
			{
				if (progress_bar && batch_processed % PMR_PROG_REPORT_INTERVAL == 0) {
					wt_end = omp_get_wtime();
					io_draw_progress_bar(read_processed + batch_processed, 
							d->n_reads, wt_end - wt_start, 1, 0, NULL);
				}

				read_processed += batch_processed;
				batch_processed = 0;

				if (reducer != NULL) reducer(d, fdata);

				/*
				INFO_TEE(d->log_fp, "Batch processed, time elapsed: %.3f\n", 
						wt_end - wt_start);
				INFO_TEE(d->log_fp, "Dictionary size: %lu/%lu\n", 
						d->kmer_dict->used, d->kmer_dict->sizemask + 1);
				*/
			}

			#pragma omp barrier
		}

		wt_end = omp_get_wtime();
		#pragma omp master
		{
			if (!progress_bar) {
				INFO_TEE(d->log_fp, "Threads finished iterating, time elapsed: %.3f\n", 
						wt_end - wt_start);
			}
		}

		} //end of omp region

		// move the partial read to the front of the buffer
		if (partial_read_size > 0) {
			memcpy(_buffer, pread, partial_read_size);
			buf_loc = _buffer + partial_read_size;
		}

		kv_destroy(read_batch);

		fastq_size -= chunk_size;
	}


	struct timeval _ntv;
	gettimeofday(&_ntv, NULL);
	double elapsed = (0.0 + _ntv.tv_sec + (double) _ntv.tv_usec/1e6) - 
				(0.0 + tv.tv_sec + (double) tv.tv_usec/1e6);

	if (progress_bar) {
		io_draw_progress_bar(0, d->n_reads, elapsed, 0, 1, NULL);
	}

	INFO_TEE(d->log_fp, "All (%9d) reads processed in %.3f seconds.\n",
			read_idx, elapsed);

	d->n_reads = read_idx;

	read_destroy(tmp_read);

	fclose(fp);
	free(_buffer);

	return NO_ERROR;
}
