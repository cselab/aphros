// Created by Petr Karnakov on 01.02.2021
// Copyright 2021 ETH Zurich

#include "util/macros.h"

#include "raw.ipp"

namespace dump {

#define XX(M)                                                                 \
  template class Raw<M>;                                                      \
  template void Raw<M>::Read(                                                 \
      FieldCell<typename M::Scal>&, const Meta& meta, std::string, M&);       \
  template void Raw<M>::Write(                                                \
      const std::string& path, const std::vector<typename M::MIdx>& starts,   \
      const std::vector<typename M::MIdx>& sizes,                             \
      const std::vector<std::vector<typename M::Scal>>& data,                 \
      typename M::MIdx global_size, Type type, const MpiWrapper& mpi);        \
  template void Raw<M>::Write(                                                \
      const FieldCell<typename M::Scal>&, const Meta& meta, std::string, M&); \
  template void Raw<M>::WriteWithXmf(                                         \
      const FieldCell<typename M::Scal>&, std::string, std::string, M&);

#define COMMA ,
#define X(dim) XX(MeshCartesian<double COMMA dim>)
MULTIDIMX

} // namespace dump
