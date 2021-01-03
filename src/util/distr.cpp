// Created by Petr Karnakov on 12.10.2020
// Copyright 2020 ETH Zurich

#ifdef _OPENMP
#include <omp.h>
#endif
#include "geom/idx.h"

#include "distr.ipp"

MpiWrapper::MpiWrapper(int* argc, const char*** argv, MPI_Comm comm)
    : comm_(comm) {
#ifdef _OPENMP
  omp_set_dynamic(0);
#endif
  int required = MPI_THREAD_FUNNELED;
  int provided;
  MPICALL(MPI_Init_thread(argc, (char***)argv, required, &provided));
  fassert_equal(required, provided, ", mismatch in thread support level");
}

MpiWrapper::~MpiWrapper() {
  MPI_Finalize();
}

MPI_Comm MpiWrapper::GetComm() const {
  return comm_;
}

int MpiWrapper::GetCommSize(MPI_Comm comm) {
  int size;
  MPICALL(MPI_Comm_size(comm, &size));
  return size;
}

int MpiWrapper::GetCommRank(MPI_Comm comm) {
  int rank;
  MPICALL(MPI_Comm_rank(comm, &rank));
  return rank;
}

int MpiWrapper::GetCommSize() const {
  return GetCommSize(comm_);
}

int MpiWrapper::GetCommRank() const {
  return GetCommRank(comm_);
}

bool MpiWrapper::IsRoot(MPI_Comm comm) {
  return GetCommRank(comm) == 0;
}

bool MpiWrapper::IsRoot() const {
  return IsRoot(comm_);
}

template class Subdomains<generic::MIdx<1>>;
template class Subdomains<generic::MIdx<2>>;
template class Subdomains<generic::MIdx<3>>;
template class Subdomains<generic::MIdx<4>>;
