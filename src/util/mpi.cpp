// Created by Petr Karnakov on 02.02.2021
// Copyright 2021 ETH Zurich

#include "mpi.h"
#include "logger.h"

#if USEFLAG(MPI)

MpiWrapper::MpiWrapper(int* argc, const char*** argv, MPI_Comm comm)
    : finalize_(true), comm_(comm) {
  int required = MPI_THREAD_FUNNELED;
  int provided;
  MPICALL(MPI_Init_thread(argc, (char***)argv, required, &provided));
  fassert_equal(required, provided, ", mismatch in thread support level");
}

MpiWrapper::~MpiWrapper() {
  if (finalize_) {
    MPI_Finalize();
  }
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

#else

MpiWrapper::MpiWrapper(int* argc, const char*** argv, MPI_Comm comm)
    : finalize_(true), comm_(comm) {}

MpiWrapper::~MpiWrapper() {}

int MpiWrapper::GetCommSize(MPI_Comm comm) {
  return 1;
}

int MpiWrapper::GetCommRank(MPI_Comm comm) {
  return 0;
}

#endif // USEFLAG(MPI)

MpiWrapper::MpiWrapper(MPI_Comm comm) : finalize_(false), comm_(comm) {}

MPI_Comm MpiWrapper::GetComm() const {
  return comm_;
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
