struct Table;
enum {NONE = -1};
struct Table *table_new(int hint,
                        int cmp(const void *x, const void *y),
                        unsigned hash(const void *key));
void table_free(struct Table *);
int table_length(struct Table *);
int table_put(struct Table *, const void *key, int value);
int table_get(struct Table *, const void *key);
int table_remove(struct Table *, const void *key);
void table_map(struct Table *,
               void apply(const void *key, int *value, void *cl),
               void *cl);
#undef T
