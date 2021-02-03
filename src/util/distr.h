// Created by Petr Karnakov on 12.10.2020
// Copyright 2020 ETH Zurich

#pragma once

#include <string>
#include <vector>

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
