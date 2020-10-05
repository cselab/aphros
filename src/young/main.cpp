// Created by Sergey Litvinov on 02.11.2018
// Copyright 2018 ETH Zurich

#include <math.h>
#include <stddef.h>

#include "young.h"

static YoungParam g_param;

#define SET(key, val) \
  if ((key) != NULL) *(key) = (val)

static double sq(double x) {
  return x * x;
}
static double cube(double x) {
  return x * x * x;
}

static int yz2sph(double y, double z, double* pr, double* pth) {
  double r, th;
  r = sqrt(sq(y) + sq(z));
  th = (r == 0) ? 0 : acos(z / r);
  if (y < 0) th = -th;
  *pr = r;
  *pth = th;
  return 0;
}

static int sph2yz(double th, double vr, double vth, double* vy, double* vz) {
  SET(vy, vr * sin(th) + vth * cos(th));
  SET(vz, vr * cos(th) - vth * sin(th));
  return 0;
}
static int u(
    YoungParam param, double r, double th, double* vr, double* vp, double* p,
    double* T) {
  double rhov, rhou, muv, muu, hv, hu, gamma1, gamma0, T1, T0;
  double R, kv, ku, av, g, au, a0, Tc;

  rhov = param.rhov;
  rhou = param.rhou;
  muv = param.muv;
  muu = param.muu;
  hv = param.hv;
  hu = param.hu;
  gamma1 = param.gamma1;
  gamma0 = param.gamma0;
  T1 = param.T1;
  T0 = param.T0;
  R = param.R;

  Tc = (3 * T1 * hv) / (2 * hv + hu);

  av = (sq(R) * Tc * gamma1 * muv) / (6 * (muv + muu));
  au = (10 * Tc * gamma1 * muu) / (3 * R * (muv + muu));
  a0 = (2 * (T0 * gamma1 + gamma0)) / R;
  kv = (cube(R) * (hv - hu)) / (2 * hv + hu);
  ku = (3 * hv) / (2 * hv + hu);
  g = (Tc * gamma1 * muv) / (R * (muv + muu) * (rhov - rhou));

  (void)av;
  (void)kv;

  SET(vr, (au * (sq(r) - sq(R)) * cos(th)) / (10 * muu));
  SET(vp, -(au * (2 * sq(r) - sq(R)) * sin(th)) / (10 * muu));
  SET(p, (-g * r * rhou * cos(th)) + au * r * cos(th) + a0);
  SET(T, T1 * ku * r * cos(th) + T0);

  return 0;
}
static int v(
    YoungParam param, double r, double th, double* vr, double* vp, double* p,
    double* T) {
  double rhov, rhou, muv, muu, hv, hu, gamma1, gamma0, T1, T0;
  double R, kv, ku, g, av, au, a0, Tc;

  rhov = param.rhov;
  rhou = param.rhou;
  muv = param.muv;
  muu = param.muu;
  hv = param.hv;
  hu = param.hu;
  gamma1 = param.gamma1;
  gamma0 = param.gamma0;
  T1 = param.T1;
  T0 = param.T0;
  R = param.R;

  Tc = (3 * T1 * hv) / (2 * hv + hu);

  av = (sq(R) * Tc * gamma1 * muv) / (6 * (muv + muu));
  au = (10 * Tc * gamma1 * muu) / (3 * R * (muv + muu));
  a0 = (2 * (T0 * gamma1 + gamma0)) / R;
  kv = (cube(R) * (hv - hu)) / (2 * hv + hu);
  ku = (3 * hv) / (2 * hv + hu);
  g = (Tc * gamma1 * muv) / (R * (muv + muu) * (rhov - rhou));

  (void)a0;
  (void)au;
  (void)ku;

  SET(vr, (av * (2 / r - (2 * sq(R)) / cube(r)) * cos(th)) / muv);
  SET(vp, -(av * (sq(R) / cube(r) + 1 / r) * sin(th)) / muv);
  SET(p, (2 * av * cos(th)) / sq(r) - g * r * rhov * cos(th));
  SET(T, T1 * (kv / sq(r) + r) * cos(th) + T0);

  return 0;
}
int young_set(YoungParam* p) {
  double rhov, rhou, muv, muu, hv, hu, gamma1, gamma0, T1;
  double T0, R;

  R = 0.2;
  rhov = 1;
  rhou = 0.04;

  muv = 0.0505964425627;
  muu = muv * 0.04;
  gamma1 = -0.01;
  gamma0 = 0.064;

  hv = hu = T1 = 1;
  T0 = 0;

  p->rhov = rhov;
  p->rhou = rhou;
  p->muv = muv;
  p->muu = muu;
  p->hv = hv;
  p->hu = hu;
  p->gamma1 = gamma1;
  p->gamma0 = gamma0;
  p->T1 = T1;
  p->T0 = T0;
  p->R = R;
  return 0;
}

int young_ini(YoungParam p) {
  g_param = p;
  return 0;
}

int young_fields(
    double y, double z, /**/ double* pvy, double* pvz, double* p, double* T) {
  double R, r, th;
  double vr, vth, vy, vz;

  R = g_param.R;
  yz2sph(y, z, &r, &th);

  if (r > R)
    v(g_param, r, th, &vr, &vth, p, T);
  else
    u(g_param, r, th, &vr, &vth, p, T);

  sph2yz(th, vr, vth, &vy, &vz);

  SET(pvy, vy);
  SET(pvz, vz);

  return 0;
}
