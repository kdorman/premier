#ifndef __MEMPOOL_H
#define __MEMPOOL_H

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#define __MEMPOOL_DEBUG 1

#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>

#include "message.h"

/**
 * Header file for memory pool.
 */

#define MEMPOOL_DEFAULT_SIZE (16 * 1024)
#define MEMPOOL_ALIGNMENT 16

#define MEMPOOL_MAX_LARGE_BLOCKS_ATTACHED 3

#define ALIGN_PTR(p, a)                                                   \
    (char *) (((uintptr_t) (p) + ((uintptr_t) a - 1)) & ~((uintptr_t) a - 1))

typedef struct mempool_data_t_ {
	char *last;
	char *end;
	struct mempool_t_ *next;
	uint32_t failed;
} mempool_data_t;

typedef struct mempool_large_t_ {
	struct mempool_large_t_ *next;
	void *alloc;
} mempool_large_t;

typedef struct mempool_t_ {
	mempool_data_t data;
	size_t max;
	struct mempool_t_ *current;
	mempool_large_t *large;
	/*void (*pool_destructor)(mempool_data_t *data);*/
#ifdef __MEMPOOL_DEBUG
	size_t allocated;
#endif
} mempool_t;


mempool_t *mempool_create(size_t size);
void mempool_destroy(mempool_t *pool);
void * mempool_alloc(mempool_t *pool, size_t size);
void * mempool_nalloc(mempool_t *pool, size_t size, size_t alignment);

#endif
