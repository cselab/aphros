#ifdef __cplusplus
extern "C" {
#endif
int aphros_main(int, const char**);
struct aphros_Parser;
struct aphros_Parser *aphros_parser_file_ini(const char *);
int aphros_parser_fin(struct aphros_Parser *);
int aphros_parser_print_vars(struct aphros_Parser*);
int aphros_parser_status(struct aphros_Parser*);
int aphros_parser_int(struct aphros_Parser *, const char *);

#ifdef __cplusplus
}
#endif
