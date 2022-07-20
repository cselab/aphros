// Created by Petr Karnakov on 04.07.2020
// Copyright 2020 ETH Zurich

#include "util/macros.h"

#if USEFLAG(MPI) && USEFLAG(HDF)
#include "hdf.ipp"
#else
#include "hdf_nompi.ipp"
#endif

namespace dump {

#define XX(M)                                                            \
  template class Hdf<M>;                                                 \
  template void Hdf<M>::Read(                                            \
      FieldCell<typename M::Scal>&, std::string, M&, std::string);       \
  template void Hdf<M>::Read(                                            \
      FieldCell<typename M::Vect>&, std::string, M&, std::string);       \
  template void Hdf<M>::Read(                                            \
      FieldCell<typename M::Expr>&, std::string, M&, std::string);       \
  template void Hdf<M>::Write(                                           \
      const FieldCell<typename M::Scal>&, std::string, M&, std::string); \
  template void Hdf<M>::Write(                                           \
      const FieldCell<typename M::Vect>&, std::string, M&, std::string); \
  template void Hdf<M>::Write(                                           \
      const FieldCell<typename M::Expr>&, std::string, M&, std::string);

#define COMMA ,
#define X(dim) XX(MeshCartesian<double COMMA dim>)
MULTIDIMX

} // namespace dump
