struct March {
    double (*f) (double[3], void *);
    void *fdata;
    int size[3];
    double spacing;
/* TODO:   double origin[3]; */
    int (*normal) (double[3], void *);
    int (*vertex) (double[3], void *);
    void *cdata;
};
int march_cube(struct March *);
