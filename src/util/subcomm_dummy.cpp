// Created by Petr Karnakov on 05.06.2020
// Copyright 2020 ETH Zurich

#include <stdexcept>

#include "subcomm.h"
#include "util/logger.h"

void SetAffinity(int) {
  fassert(false, "not implemented");
}

void PrintStats(MPI_Comm, MPI_Comm, MPI_Comm) {
  fassert(false, "not implemented");
}

void SubComm(MPI_Comm&, MPI_Comm&, MPI_Comm&) {
  fassert(false, "not implemented");
}
