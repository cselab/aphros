// Created by Petr Karnakov on 20.07.2018
// Copyright 2018 ETH Zurich

#include <mpi.h>
#include <unistd.h>
#include <fstream>
#include <sstream>

#include "sysinfo.h"
#include "util/logger.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#if USEFLAG(AMGX)
#include <cuda_runtime.h>
#endif

namespace sysinfo {

// Returns value from /proc/self/status found by key k
size_t GetProcStat(std::string k) {
  std::ifstream f("/proc/self/status");

  f >> std::skipws;
  while (f) {
    std::string s;
    f >> s;
    if (s == k + ":") {
      size_t r = 0;
      f >> r;
      return r;
    }
  }
  return 0;
}

size_t GetMem() {
  return GetProcStat("VmRSS") << 10;
}

bool HasHyperthreads() {
  bool has_ht = false;
  std::ifstream f("/proc/cpuinfo");
  f >> std::skipws;
  while (f) {
    std::string s;
    f >> s;
    if (s == "flags") {
      break;
    }
  }
  std::string line;
  std::getline(f, line);
  std::istringstream is(line);
  while (!is.eof()) {
    std::string s;
    is >> s;
    if (s == "ht") {
      has_ht = true;
      break;
    }
  }
  return has_ht;
}

std::string GetHostname() {
  const size_t kMaxLength = 4096;
  char buf[kMaxLength];
  int err = gethostname(buf, kMaxLength);
  fassert_equal(err, 0);
  return buf;
}

Info GetInfo(InfoSelect select) {
  Info info;
  if (select.mpi) {
    MPI_Comm_size(MPI_COMM_WORLD, &info.comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &info.comm_rank);
  }

#ifdef _OPENMP
  if (select.openmp) {
    info.max_threads = omp_get_max_threads();
  }
#endif

#if USEFLAG(AMGX)
  if (select.cuda) {
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    info.cuda_enabled = true;
    info.cuda_uuid = *(uint16_t*)(&prop.uuid);
    info.cuda_mem = prop.totalGlobalMem;
  }
#endif

  info.hostname = GetHostname();
  return info;
}

} // namespace sysinfo
