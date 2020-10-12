// Created by Petr Karnakov on 12.10.2020
// Copyright 2020 ETH Zurich

#include <algorithm>
#include <iostream>
#include <vector>

#include "format.h"
#include "logger.h"

#include "distr.h"

template <class MIdx>
std::vector<MIdx> Subdomains<MIdx>::GetValidProcs(
    MIdx mesh_size, MIdx block_size, size_t nproc) {
  const auto ms = mesh_size;
  const auto bs = block_size;

  if (mesh_size % block_size != MIdx(0)) {
    return {};
  }
  const auto gblocks = ms / bs; // global blocks
  auto is_valid = [&](MIdx procs) -> bool {
    if (int(nproc) != procs.prod()) {
      return false;
    }
    if (gblocks % procs != MIdx(0)) {
      return false;
    }
    return true;
  };

  std::vector<size_t> divisors; // divisors of nproc
  for (size_t px = 1; px <= nproc; ++px) {
    if (nproc % px == 0) {
      divisors.push_back(px);
    }
  }
  std::vector<MIdx> res;
  for (auto px : divisors) {
    for (auto py : divisors) {
      auto pz = nproc / (px * py);
      const MIdx procs(px, py, pz);
      if (is_valid(procs)) {
        res.push_back(procs);
      }
    }
  }
  return res;
}

template <class MIdx>
Subdomains<MIdx>::Subdomains(MIdx mesh_size, MIdx block_size, size_t nproc) {
  info.mesh_size = mesh_size;
  info.block_size = block_size;
  info.nproc = nproc;
  const auto ms = mesh_size;
  const auto bs = block_size;

  fassert(ms > MIdx(0), util::Format("Mesh size {} must be positive", ms));

  fassert(
      ms % bs == MIdx(0),
      util::Format("Mesh size {} not divisible by block size {}", ms, bs));

  auto valid_procs = GetValidProcs(mesh_size, block_size, nproc);
  const auto gblocks = ms / bs;

  auto get_area = [&](MIdx procs) {
    auto blocks = gblocks / procs;
    size_t res = 0;
    for (size_t i = 0; i < MIdx::dim; ++i) {
      res += (bs.prod() / bs[i]) * (blocks.prod() / blocks[i]);
    }
    return res;
  };

  auto pos = std::min_element(
      valid_procs.begin(), valid_procs.end(),
      [&](MIdx pa, MIdx pb) { return get_area(pa) < get_area(pb); });

  fassert(
      pos != valid_procs.end(),
      util::Format(
          "No valid partitions found for mesh size {}, block size {}, nproc {}",
          ms, bs, nproc));
  const MIdx procs = *pos;
  info.procs = procs;
  info.blocks = gblocks / procs;
  info.area = get_area(procs);
}

template <class MIdx>
std::string Subdomains<MIdx>::GetConfig() const {
  std::string res;
  res += "set int px " + std::to_string(info.procs[0]);
  res += "\nset int py " + std::to_string(info.procs[1]);
  res += "\nset int pz " + std::to_string(info.procs[2]);

  res += "\n\nset int bx " + std::to_string(info.blocks[0]);
  res += "\nset int by " + std::to_string(info.blocks[1]);
  res += "\nset int bz " + std::to_string(info.blocks[2]);

  res += "\n\nset int bsx " + std::to_string(info.block_size[0]);
  res += "\nset int bsy " + std::to_string(info.block_size[1]);
  res += "\nset int bsz " + std::to_string(info.block_size[2]);
  res += "\n";
  return res;
}
