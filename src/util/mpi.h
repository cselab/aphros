// Created by Petr Karnakov on 02.02.2021
// Copyright 2021 ETH Zurich

#pragma once

#include "macros.h"

#if USEFLAG(MPI)

#include <mpi.h>
#include <stdexcept>
#include <string>
#define MPICALL(x)                                                    \
  do {                                                                \
    int errorcode;                                                    \
    errorcode = x;                                                    \
    if (errorcode != MPI_SUCCESS) {                                   \
      char string[MPI_MAX_ERROR_STRING];                              \
      int resultlen;                                                  \
      MPI_Error_string(errorcode, string, &resultlen);                \
      throw std::runtime_error(                                       \
          std::string() + __FILE__ + ":" + std::to_string(__LINE__) + \
          ": mpi failed: " + string);                                 \
    }                                                                 \
  } while (0)

#else

using MPI_Comm = int;
#define MPI_COMM_WORLD 0
#define MPICALL(x)

#endif // USEFLAG(MPI)

#include <stdexcept>

class MpiWrapper {
 public:
  MpiWrapper(int* argc, const char*** argv, MPI_Comm comm = MPI_COMM_WORLD);
  MpiWrapper(MPI_Comm comm = MPI_COMM_WORLD);
  MpiWrapper(const MpiWrapper&) = delete;
  MpiWrapper(MpiWrapper&&) = delete;
  ~MpiWrapper();
  MpiWrapper& operator=(const MpiWrapper&) = delete;
  MpiWrapper& operator=(MpiWrapper&&) = delete;
  MPI_Comm GetComm() const;
  static int GetCommSize(MPI_Comm);
  static int GetCommRank(MPI_Comm);
  static bool IsRoot(MPI_Comm);
  int GetCommSize() const;
  int GetCommRank() const;
  bool IsRoot() const;

 private:
  bool finalize_;
  MPI_Comm comm_;
};
