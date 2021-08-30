// Created by Petr Karnakov on 20.07.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace sysinfo {

// Returns resident memory [bytes]
size_t GetMem();
// Returns true if Intel Hyperthreads are enabled
bool HasHyperthreads();

std::string GetHostname();

struct Info {
  int comm_rank = 0;
  int comm_size = 1;

  int omp_num_threads = 1;
  int omp_max_threads = 1;

  bool cuda_enabled = false;
  uint16_t cuda_uuid;
  size_t cuda_mem;

  std::string hostname;
};
struct InfoSelect {
  bool mpi = true;
  bool openmp = true;
  bool cuda = true;
};
Info GetInfo(InfoSelect);

struct Misc {
  Misc() = default;
  Misc(int argc, const char** argv);

  int argc = 0;
  int arg_after_double_dash = 0; // Index immediately after first `--`
  const char** argv = nullptr;
};

extern Misc misc; // Miscellaneous global variables

} // namespace sysinfo
