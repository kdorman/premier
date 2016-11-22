/**
 * @file message.c
 * Functions for outputting error and debugging messages.
 */

#include <time.h>
#include <sys/time.h>
#include "message.h"

/**
 * Write uniformly formatted message to requested stream.
 * Write uniformly formatted debug or error messages to requested stream.
 * Certain common errors/problems issue default messages so they need not be
 * typed repeatedly and nonstandardly.  Returns #format_id so that the function
 * can be conveniently called as return message(...) or exit(message(...)).
 *
 * @param fp file handle where message to be written
 * @param file_name name of offending source file
 * @param fxn_name name of offending function
 * @param line line number of offending function
 * @param format_type urgency of message
 * @param format_id identify one of a default set of canned messages
 * @param format formatted string to output via vprintf
 * @return format_id
 */
int message(FILE *fp, const char *file_name, int line,
	int msg_type, int msg_id, const char *format, ...) {

	static char timestr[32];
	time_t rawtime;
	struct timeval tv;

	time(&rawtime);
	gettimeofday(&tv, NULL);

	struct tm *timeinfo = localtime(&rawtime);
	snprintf(timestr, 31, "%4d-%02d-%02d %02d:%02d:%02d,%03d",
			1900+timeinfo->tm_year, timeinfo->tm_mon, timeinfo->tm_mday,
			timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec,
			(int) tv.tv_usec / 1000);
	//strftime(timestr, 31, "%c", timeinfo);

	va_list args;
	fprintf(fp, "%s [%s] (%s:%d): ",
		timestr,
		msg_type == INFO_MSG ? "INFO" : msg_type == DEBUG_MSG ? "DEBUG"
			: msg_type == WARNING_MSG ? "WARNING" : "ERROR", 
		file_name, line);

	if (msg_id == NO_ERROR) {
		va_start(args, format);
		vfprintf(fp, format, args);
		va_end(args);
	} else {
		switch(msg_id) {
			case MEMORY_ALLOCATION: 
				if (format) {
					fprintf(fp, "could not allocate ");
					va_start(args, format);
					vfprintf(fp, format, args);
					va_end(args);
				} else
					fprintf(fp, "memory allocation error\n");
				break;
			case INVALID_CMD_OPTION:
				fprintf(fp, "unrecognized command option");
				if (format) {
					fprintf(fp, ": ");
					va_start(args, format);
					vfprintf(fp, format, args);
					va_end(args);
				} else
					fprintf(fp, "\n");
				break;
			case INVALID_CMD_ARGUMENT:
				fprintf(fp, "invalid argument to command option");
				if (format) {
					fprintf(fp, ": ");
					va_start(args, format);
					vfprintf(fp, format, args);
					va_end(args);
				} else
					fprintf(fp, "\n");
				break;
			case FILE_NOT_FOUND:
				fprintf(fp, "file \"%s\" not found\n", format);
				break;
			case FILE_FORMAT_ERROR:
				fprintf(fp, "invalid file format in \"%s\"\n", format);
				break;
			case END_OF_FILE:
				fprintf(fp, "unexpected end of file in file \"%s\"\n", format);
				break;
			case INTERNAL_MISMATCH:
				fprintf(fp, "[internal mismatch] %s\n", format);
				break;
			default:
				va_start(args, format);
				vfprintf(fp, format, args);
				va_end(args);
		}
	}
	return msg_id;
} /* message */
