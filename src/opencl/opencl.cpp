// Created by Petr Karnakov on 14.04.2021
// Copyright 2021 ETH Zurich

#include "opencl.ipp"

#define X(dim) template class OpenCL<MeshCartesian<double, dim>>;
MULTIDIMX

