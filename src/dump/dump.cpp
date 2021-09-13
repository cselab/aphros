// Created by Petr Karnakov on 10.07.2020
// Copyright 2020 ETH Zurich

#include "dump/dump.ipp"

namespace dump {

template std::vector<std::pair<std::string, std::vector<double>>> ReadCsv(
    std::istream&, char delim);
template std::vector<std::pair<std::string, std::vector<double>>> ReadCsv(
    const std::string&, char delim);

template void DumpCsv<double>(
    const std::vector<std::pair<std::string, std::vector<double>>>& indata,
    std::string path, char delim);

#define XX(M)                                                             \
  template void DumpCsv<M>(                                               \
      const std::vector<                                                  \
          std::pair<std::string, std::vector<typename M::Scal>>>& indata, \
      std::string path, M& m, char delim);

#define COMMA ,
#define X(dim) XX(MeshCartesian<double COMMA dim>)
MULTIDIMX
#undef X
#undef XX

} // namespace dump
