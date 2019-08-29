struct Sample;
struct Sample* sample_ini(void);
void sample_fin(struct Sample*);
double sample_f(struct Sample*, double, double, double);
int sample_time(struct Sample*, double);
int sample_next(struct Sample*);
