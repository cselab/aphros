// Created by Sergey Litvinov on 06.01.2020
// Copyright 2020 ETH Zurich

struct Adamar {
  double mu0;
  double mu1;
  double a;
  double pi;

  double U;
  double distance;
};

int adamar_fields(
    struct Adamar*, double x, double y, double z, double* vx, double* vy,
    double* vz, double* p);
