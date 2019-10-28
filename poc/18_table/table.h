struct Table;
struct Table *table_new(int hint,
                        int cmp(const void *x, const void *y),
                        unsigned hash(const void *key));
void table_free(struct Table *);
int table_length(struct Table *);
void *table_put(struct Table *, const void *key, void *value);
void *table_get(struct Table *, const void *key);
void *table_remove(struct Table *, const void *key);
void table_map(struct Table *,
               void apply(const void *key, void **value, void *cl),
               void *cl);
#undef T
