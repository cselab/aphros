struct Vec {
    double x;
    double y;
    double z;
};
struct Vec;
double GetOffset(double, double);
void GetColor(struct Vec *, struct Vec *);
void PrintHelp(void);
void GetNormal(struct Vec *, double, double, double,
	       double (*f) (double, double, double, void *), void *);
