#ifndef __PMR_ORACLE_H__
#define __PMR_ORACLE_H__

#ifdef __cplusplus
#define EXTERNC extern "C"
#else
#define EXTERNC
#endif

EXTERNC void run_oracle(data *d);
EXTERNC void oracle_init_em(data *d);

#undef EXTERNC

#endif
