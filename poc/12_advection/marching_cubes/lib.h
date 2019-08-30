struct Vec {
    double x;
    double y;
    double z;
};
double GetOffset(double, double);
void GetColor(double, double, double, double*, double*, double*);
void PrintHelp(void);
void GetNormal(struct Vec *, double, double, double,
	       double (*f) (double, double, double, void *), void *);
