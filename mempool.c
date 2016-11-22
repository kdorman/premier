#include "mempool.h"

static void * mempool_alloc_block(mempool_t *pool, size_t size);
static void * mempool_alloc_large(mempool_t *pool, size_t size);
static void * _mempool_memalign(size_t alignment, size_t size);

mempool_t *mempool_create(size_t size)
{
	mempool_t *pool;
	static int mp_pagesize = 0;

	if (!mp_pagesize)
		mp_pagesize = getpagesize();

	pool = (mempool_t *) _mempool_memalign(MEMPOOL_ALIGNMENT, size);

	if (pool == NULL)
		return NULL;

	/* Full pool size minus the space occupied by mempool_t struct itself */
	pool->data.last = (char *) pool + sizeof(mempool_t);
	pool->data.end = (char *) pool + size;
	pool->data.next = NULL;
	pool->data.failed = 0;

	size = size - sizeof(mempool_t);
    pool->max = (size < mp_pagesize - 1) ? size : mp_pagesize - 1; 

	pool->current = pool;
	pool->large = NULL;

#ifdef __MEMPOOL_DEBUG
	pool->allocated = 0;
#endif

	return pool;
}

void mempool_destroy(mempool_t *pool)
{
	mempool_large_t *large;
	mempool_t *p, *n;

	/* release all large memory blocks */
	for (large = pool->large; large; large = large->next){
		if (large->alloc)
			free(large->alloc);
	}

	for (p = pool, n = pool->data.next; ; p = n, n = n->data.next){
		free(p);

		if (n == NULL)
			break;
	}
}

void * mempool_alloc(mempool_t *pool, size_t size)
{
	char *m;
	mempool_t *p;

#ifdef __MEMPOOL_DEBUG
	pool->allocated += size;
#endif
	if (size <= pool->max){
		p = pool->current;	/* find current pool that is available */
		do {
			m = ALIGN_PTR(p->data.last, MEMPOOL_ALIGNMENT);

			/* find enough space to allocate */
			if ((size_t) (p->data.end - m) >= size){
#ifdef __MEMPOOL_DEBUG
				pool->allocated += (size_t) ((intptr_t) m - 
						(intptr_t) p->data.last);
#endif
				p->data.last = m + size;
				return m;
			}

			p = p->data.next;
		} while (p);

		return mempool_alloc_block(pool, size);
	}
	else
		return mempool_alloc_large(pool, size);

}

/* allocate memory block by given alignment size, which must be a power of 2 */
void * mempool_nalloc(mempool_t *pool, size_t size, size_t alignment)
{
	char *m;
	mempool_t *p;

	if (size <= pool->max) {
		p = pool->current;

		do {
			m = ALIGN_PTR(p->data.last, alignment);

			if ((size_t) (p->data.end - m) >= size) {
				p->data.last = m + size;
				return m;
			}

			p = p->data.next;
		} while(p);

		return mempool_alloc_block(pool, size);
	}

	return mempool_alloc_large(pool, size);
}

static void * mempool_alloc_block(mempool_t *pool, size_t size)
{
	char *m;
	size_t pool_size;
    mempool_t *p, *new, *current;

    pool_size = (size_t) (pool->data.end - (char *) pool);

	/* create a new pool of the same size */
	m = _mempool_memalign(MEMPOOL_ALIGNMENT, pool_size);
	if (m == NULL)
		return NULL;

    new = (mempool_t *) m;

    new->data.end = m + pool_size;
    new->data.next = NULL;
    new->data.failed = 0;

	/* allocate memory block on new pool */
    m += sizeof(mempool_data_t);
    m = ALIGN_PTR(m, MEMPOOL_ALIGNMENT);
    new->data.last = m + size;

    current = pool->current;

	/* increase all failed counter for previous pools by 1.
	 * forgo allocating on any pool with more than 4 times of failures.
	 */
    for (p = current; p->data.next; p = p->data.next) {
        if (p->data.failed++ > 4) {
            current = p->data.next;
        }
    }

	/* insert current pool by the end of the chain list */
    p->data.next = new;

	/* decide which pool should be the default one for allocation */
    pool->current = current ? current : new;

    return m;
}

static void * mempool_alloc_large(mempool_t *pool, size_t size)
{
	void *p;
	mempool_large_t *large;
	int n;

	p = calloc(1, size);
	if (p == NULL){
		message(stderr, __FILE__, __LINE__, ERROR_MSG,
				MEMORY_ALLOCATION, NULL);
		return NULL;
	}

	n = 0;

	for (large = pool->large; large; large = large->next) {
		/* try find an available spot to attach the memory block 
		 * (previously released?)
		 */
		if (large->alloc == NULL){
			large->alloc = p;
			return p;
		}

		if (n++ > MEMPOOL_MAX_LARGE_BLOCKS_ATTACHED)
			break;
	}

	/* if we exceed the maximum number of attached large memory blocks, try
	 * allocate a pointer to large memory block within pool */
	large = mempool_alloc(pool, sizeof(mempool_large_t));
	if (large == NULL){
		free(p);
		message(stderr, __FILE__, __LINE__, ERROR_MSG,
				MEMORY_ALLOCATION, NULL);
		return NULL;
	}

	/* and place the next large memory block at the head of the chain list */
	large->alloc = p;
	large->next = pool->large;
	pool->large = large;

	return p;
}

/* release large memory block */
int mempool_free(mempool_t *pool, void *p)
{
	mempool_large_t *large;
	for (large = pool->large; large; large = large->next) {
		if (p == large->alloc){
			free(p);
			large->alloc = NULL;	/* set to NULL for further attachment */
			return NO_ERROR;
		}
	}

	return RESOURCE_NOT_FOUND; 
}

/* --- privation functions --- */
static void * _mempool_memalign(size_t alignment, size_t size)
{
	void *p;
	int err;

	err = posix_memalign(&p, alignment, size);
	if (err)
		return NULL;

	return p;
}
