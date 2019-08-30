struct March {
    double (*f) (double, double, double, void *);
    void *fdata;
    int size[3];
    double spacing;
/* TODO:   double origin[3]; */
    int (*normal) (double, double, double, void *);
    int (*vertex) (double, double, double, void *);
    void *cdata;
};
int march_cube(struct March *);
int march_tetrahedron(struct March *);
