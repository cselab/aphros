enum { CSV_MAX_NF = 999 };
struct CSV {
  char* name[CSV_MAX_NF];
  double* data[CSV_MAX_NF];
  int nf, nr;
};

struct CSV* csv_ini(int);
struct CSV* csv_read(FILE*);
int csv_fin(struct CSV*);
double* csv_field(struct CSV*, const char*);
int csv_nf(struct CSV*);
int csv_nr(struct CSV*);
int csv_write(struct CSV*, FILE*);
int csv_write_space(struct CSV*, FILE*);
int csv_add(struct CSV*, const char*);
