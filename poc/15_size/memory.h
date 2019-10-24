#define MALLOC(n, p)							\
    do {								\
    *(p) = malloc(n*sizeof(**(p)));					\
    if (*(p) == NULL)							\
	ERR(("alloc failed, n = %d", n));				\
    } while(0)

#define FREAD(n, p, f)							\
    do {								\
	if ((int)fread(p, sizeof(*(p)), (n), (f)) != (n))		\
	    ERR(("fread failed, n = %d", n));				\
    } while(0)
