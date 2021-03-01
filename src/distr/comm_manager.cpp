// Created by Petr Karnakov on 02.01.2021
// Copyright 2021 ETH Zurich

#include "util/macros.h"

#if USEFLAG(MPI)
#include "comm_manager.ipp"
#else
#include "comm_manager_seq.ipp"
#endif

template struct CommManager<1>;
template struct CommManager<2>;
template struct CommManager<3>;
template struct CommManager<4>;
