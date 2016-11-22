#ifndef __PMR_KMER_DICT_H__
#define __PMR_KMER_DICT_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <stdint.h>

enum {DICT_OK, DICT_ERR, DICT_KEY_EXISTS, DICT_ADDED, DICT_REPLACED};

#define DICT_MINSIZE 8
#define PERTURB_SHIFT 5

#define DICT_SIZE(d) ((d)->sizemask + 1)
#define DICT_COMPARE_KEYS(key1, key2) ((key1) == (key2))
#define DICT_HASH_KEY(key) _dict_uint64_hash_raw((key))

typedef uint64_t dictsize_t;
typedef uint64_t hash_t;

typedef struct dicttype_ dicttype;
typedef struct dict_ dict;
typedef struct dictentry_ dictentry;

struct dictentry_ {
	void *key;
	void *value;
};

/*
struct dicttype_ {
	hash_t (*hashfunc)(const void *key);
	int (*keycomp)(const void *key1, const void *key2);
};
*/

struct dict_ {
	dictsize_t sizemask;
	dictsize_t used;
	/*
	dicttype type;
	*/
	dictentry *table;
};


/* ------ Functions ----- */
dict *dict_create(register dictsize_t init_size);
void dict_destroy(dict *d);
void dict_reset(dict *d);
void dict_try_resize(dict *d);
int dict_add(dict *d, void *key, void *value);
dictentry *dict_add_find(dict *d, void *key, void *val);
dictentry *dict_addkey_find(dict *d, void *key);
dictentry *dict_find(dict *d, void *key);

#endif
