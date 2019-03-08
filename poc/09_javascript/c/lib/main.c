typedef double Scal;
static Scal n[2*25], a[25];

int segment_ini(void) {
    return 0;
}

int segment_get(const Scal alpha[25], /**/ Scal **pn, Scal **pa) {
    *pn = n; *pa = a;
    return 0;
}
