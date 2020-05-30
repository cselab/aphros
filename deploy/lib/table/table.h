struct Table;
enum { TABLE_EMPY = -1 };
struct Table* table_ini(int hint);
int table_fin(struct Table*);
int table_length(struct Table*);
int table_put(struct Table*, int key, int value);
int table_get(struct Table*, int key, int* value);
int table_remove(struct Table*, int key);
int* table_array(struct Table*);

#undef T
