#define MALLOC(n, p)							\
    do {								\
	*(p) = malloc(n*sizeof(**(p)));					\
	if (*(p) == NULL) {						\
	    fprintf(stderr, "%s:%d: alloc failed, n = %d\n", __FILE__, __LINE__, n); \
	    exit(2);							\
	}								\
    } while(0)

#define FREAD(n, p, f) \
    do {								\
	if ((int)fread(p, sizeof(*(p)), (n), (f)) != (n)) {		\
	    fprintf(stderr, "%s:%d: failt to read, n = %d\n", __FILE__, __LINE__, n); \
	    exit(2);							\
	}								\
} while(0)
