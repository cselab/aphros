// Created by Petr Karnakov on 01.02.2021
// Copyright 2021 ETH Zurich

#include "util/macros.h"

#include "raw.ipp"

namespace dump {

#define XX(M)                                                                  \
  template class Raw<M>;                                                       \
  template void Raw<M>::Read(                                                  \
      FieldCell<typename M::Scal>&, const Meta& meta, const std::string&, M&); \
  template void Raw<M>::Write(                                                 \
      const std::string& path, const std::vector<typename M::MIdx>& starts,    \
      const std::vector<typename M::MIdx>& sizes,                              \
      const std::vector<std::vector<typename M::Scal>>& data,                  \
      typename M::MIdx global_size, Type type, const MpiWrapper& mpi,          \
      bool nompi);                                                             \
  template void Raw<M>::Write(                                                 \
      const FieldCell<typename M::Scal>&, const Meta& meta,                    \
      const std::string&, M&);                                                 \
  template void Raw<M>::WriteWithXmf(                                          \
      const FieldCell<typename M::Scal>&, const std::string&,                  \
      const std::string&, M&);                                                 \
  template void Raw<M>::WritePlainArrayWithXmf(                                \
      const std::string& rawpath, const std::string& fieldname,                \
      const typename M::Scal* data, MIdx size);

#define COMMA ,
#define X(dim) XX(MeshCartesian<double COMMA dim>)
MULTIDIMX

} // namespace dump
