struct March {
/* TODO:   double origin[3]; */
    int (*vertex) (double[3], void *);
    void *cdata;
};
int march_cube(struct March *, double *cube);
