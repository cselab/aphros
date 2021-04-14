// Created by Petr Karnakov on 10.07.2020
// Copyright 2020 ETH Zurich

#include <fstream>

#include "geom/mesh.h"

template <class Scal>
void DumpCsv(
    const std::vector<std::pair<std::string, std::vector<Scal>>>& data,
    std::string path) {
  if (data.empty()) {
    return;
  }
  for (auto& d : data) {
    fassert_equal(d.second.size(), data[0].second.size());
  }
  std::ofstream o(path);
  o.precision(16);
  {
    bool first = true;
    for (auto& d : data) {
      first || o << ',', first = false;
      o << d.first;
    }
  }
  o << std::endl;
  for (size_t i = 0; i < data[0].second.size(); ++i) {
    bool first = true;
    for (auto& d : data) {
      first || o << ',', first = false;
      o << d.second[i];
    }
    o << "\n";
  }
}

template <class M>
void DumpCsv(
    const std::vector<std::pair<std::string, std::vector<typename M::Scal>>>&
        indata,
    std::string path, M& m) {
  using Scal = typename M::Scal;
  auto sem = m.GetSem("dumpcsv");
  struct {
    std::vector<std::pair<std::string, std::vector<Scal>>> data;
  } * ctx(sem);
  auto& data = ctx->data;
  if (sem("gather")) {
    for (auto& d : indata) {
      data.push_back({d.first, d.second});
    }
    for (auto& d : data) {
      m.Reduce(&d.second, Reduction::concat);
    }
  }
  if (sem("write") && m.IsRoot()) {
    DumpCsv(ctx->data, path);
  }
}

using Scal = double;
using M = MeshCartesian<Scal, 3>;

template void DumpCsv<Scal>(
    const std::vector<std::pair<std::string, std::vector<Scal>>>& indata,
    std::string path);

template void DumpCsv<M>(
    const std::vector<std::pair<std::string, std::vector<typename M::Scal>>>&
        indata,
    std::string path, M& m);
