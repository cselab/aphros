#define T Table
typedef struct T *T;
extern T table_new(int hint,
                   int cmp(const void *x, const void *y),
                   unsigned hash(const void *key));
extern void table_free(T *);
extern int table_length(T);
extern void *table_put(T, const void *key, void *value);
extern void *table_get(T, const void *key);
extern void *table_remove(T, const void *key);
extern void table_map(T,
                      void apply(const void *key, void **value, void *cl),
                      void *cl);
extern void **table_toArray(T, void *end);
