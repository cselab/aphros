void *memory_malloc(int);
void memory_free(void *ptr);

#define MALLOC(n, p)							\
    do {								\
	*(p) = memory_malloc(n*sizeof(**(p)));				\
	if (*(p) == NULL)						\
	    ERR(("alloc failed, n = %d", n));				\
    } while(0)

#define FREE(p) memory_free(p);
