// Created by Petr Karnakov on 12.10.2020
// Copyright 2020 ETH Zurich

#include <algorithm>
#include <iostream>
#include <vector>

#include "format.h"
#include "geom/dir.h"
#include "geom/range.h"
#include "logger.h"

#include "distr.h"

template <class MIdx>
std::vector<MIdx> Subdomains<MIdx>::GetValidProcs(
    MIdx mesh_size, MIdx block_size, size_t nproc) {
  const size_t dim = MIdx::dim;
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

  auto divisors = [](size_t n) {
    std::vector<size_t> res;
    for (size_t d = 1; d * d <= n; ++d) {
      if (n % d == 0) {
        res.push_back(d);
        const size_t r = n / d;
        if (r != d) {
          res.push_back(r);
        }
      }
    }
    return res;
  };
  std::vector<MIdx> res;
  for (auto px : divisors(nproc)) {
    for (auto py : divisors(dim > 1 ? nproc / px : 1)) {
      for (auto pz : divisors(dim > 2 ? nproc / (px * py) : 1)) {
        for (auto pw : divisors(dim > 3 ? nproc / (px * py * pz) : 1)) {
          MIdx procs;
          if (dim > 0) procs[0] = px;
          if (dim > 1) procs[1] = py;
          if (dim > 2) procs[2] = pz;
          if (dim > 3) procs[3] = pw;
          if (is_valid(procs)) {
            res.push_back(procs);
          }
        }
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
  static const size_t dim = MIdx::dim;
  std::string res;

  auto dirs = generic::Range<size_t>(dim);
  auto cmd = [&](std::string prefix, size_t d, int value) {
    return "set int " + prefix + GDir<dim>(d).letter() + ' ' +
           std::to_string(value) + '\n';
  };
  for (auto d : dirs) {
    res += cmd("p", d, info.procs[d]);
  }
  for (auto d : dirs) {
    res += cmd("b", d, info.blocks[d]);
  }
  for (auto d : dirs) {
    res += cmd("bs", d, info.block_size[d]);
  }
  return res;
}
