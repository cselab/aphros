struct March {
    double (*f) (double, double, double, void *);
    void *fdata;
    int n;
    double spacing;
/* TODO:   double origin[3]; */
    int (*normal)(double, double, double, void *);
    int (*vertex)(double, double, double, void *);
    void *cdata;
};
