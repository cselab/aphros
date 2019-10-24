void *memory_malloc(int);
void memory_free(void *ptr);
void *memory_realloc(void*, int);

#define MALLOC(n, p)							\
    do {								\
	*(p) = memory_malloc(n*sizeof(**(p)));				\
	if (*(p) == NULL)						\
	    ERR(("alloc failed, n = %d", n));				\
    } while(0)

#define REALLOC(n, p)							\
    do {								\
	*(p) = memory_realloc(n*sizeof(**(p)));				\
	if (*(p) == NULL)						\
	    ERR(("realloc failed, n = %d", n));				\
    } while(0)

#define FREE(p) memory_free(p);
