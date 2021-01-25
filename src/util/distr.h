// Created by Petr Karnakov on 12.10.2020
// Copyright 2020 ETH Zurich

#pragma once

#include <mpi.h>
#include <stdexcept>
#include <string>
#include <vector>

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

#define MPICALL(x)                                                    \
  do {                                                                \
    int errorcode;                                                    \
    errorcode = x;                                                    \
    if (errorcode != MPI_SUCCESS) {                                   \
      char string[MPI_MAX_ERROR_STRING];                              \
      int resultlen;                                                  \
      MPI_Error_string(errorcode, string, &resultlen);                \
      throw std::runtime_error(FILELINE + ": mpi failed: " + string); \
    }                                                                 \
  } while (0)

// Description of subdomains with blocks.
template <class MIdx>
struct Subdomains {
  // Partitions mesh to equal rectangular subdomains with blocks
  // to minimize the communication area.
  // mesh_size: global mesh size in cells
  // block_size: block size in cells
  // nproc: number of processors
  Subdomains(MIdx mesh_size, MIdx block_size, size_t nproc);

  // Generates config for DistrSolver.
  // Example:
  // with `blocks_size={16, ...}`, `procs={2, ...}` and `blocks={4, ...}`
  // Returns:
  // set int bsx 16
  // set int px 2
  // set int bx 4
  // ... (same for y and z)
  std::string GetConfig() const;

  static std::vector<MIdx> GetValidProcs(
      MIdx mesh_size, MIdx block_size, size_t nproc);

  struct Info {
    MIdx mesh_size; // mesh size
    MIdx block_size; // block size
    size_t nproc; // total number of processors
    MIdx procs; // number of processors in each direction
    MIdx blocks; // number of blocks in the subdomain for each processor
    size_t area; // communication area in cells per processor
  };
  Info info;
};
