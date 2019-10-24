struct CSV {
  char **name;
  double **data;
  int nf, nr;
};

int csv_read(FILE *, struct CSV *);
int csv_fin(struct CSV *);
double *csv_get(struct CSV *, const char *);
int csv_delete(struct CSV *, const char *);
int csv_write(struct CSV *, FILE *);
