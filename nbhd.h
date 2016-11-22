#ifndef __PMR_NBHD_H__
#define __PMR_NBHD_H__

#ifdef __cplusplus
#define EXTERNC extern "C"
#else
#define EXTERNC
#endif

EXTERNC void hmm_build_1st_nbhd(data *d, read_t *read, void *fdata);
EXTERNC void hmm_determine_dmax(data *d, read_t *read, void *fdata);
EXTERNC void hmm_adaptive_dmax(data *d, read_t *read, void *fdata);

#undef EXTERNC

#endif
