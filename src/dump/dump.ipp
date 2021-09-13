// Created by Petr Karnakov on 10.08.2021
// Copyright 2021 ETH Zurich

#include <fstream>
#include <sstream>

#include "dump/dump.h"
#include "geom/mesh.h"
#include "util/format.h"

namespace dump {

template <class Scal>
std::vector<std::pair<std::string, std::vector<Scal>>> ReadCsv(
    std::istream& in, char delim) {
  std::string line;
  std::getline(in, line);
  const std::vector<std::string> header = SplitByDelimiter(line, delim);
  std::vector<std::pair<std::string, std::vector<Scal>>> res;
  for (auto h : header) {
    res.push_back({h, {}});
  }
  int lineno = 1;
  while (std::getline(in, line)) {
    ++lineno;
    const std::vector<std::string> values = SplitByDelimiter(line, delim);
    fassert_equal(
        header.size(), values.size(),
        util::Format(" while processing line {:}\n{:}", lineno, line));
    for (size_t i = 0; i < header.size(); ++i) {
      res[i].second.push_back(std::stod(values[i]));
    }
  }
  return res;
}

template <class Scal>
std::vector<std::pair<std::string, std::vector<Scal>>> ReadCsv(
    const std::string& path, char delim) {
  std::ifstream fin(path);
  fassert(fin.good(), "Can't open file '" + path + "' for reading");
  return ReadCsv<Scal>(fin, delim);
}

template <class Scal>
void DumpCsv(
    const std::vector<std::pair<std::string, std::vector<Scal>>>& data,
    std::string path, char delim) {
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
      first || o << delim, first = false;
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
    std::string path, M& m, char delim) {
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
    DumpCsv(ctx->data, path, delim);
  }
}

} // namespace dump
