// Created by Sergey Litvinov on 02.11.2018
// Copyright 2018 ETH Zurich

#ifdef __cplusplus
extern "C" {
#endif

typedef struct YoungParam YoungParam;
struct YoungParam {
  double rhov, rhou, muv, muu, hv, hu, gamma1, gamma0, T1, T0;
  double R;
};

int young_set(YoungParam*);
int young_ini(YoungParam);
int young_fields(
    double, double, /**/ double* vy, double* vz, double* p, double* T);

#ifdef __cplusplus
}
#endif
