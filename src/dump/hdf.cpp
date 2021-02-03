// Created by Petr Karnakov on 04.07.2020
// Copyright 2020 ETH Zurich

#include "util/macros.h"

#if USEFLAG(MPI) && USEFLAG(HDF)
#include "hdf.ipp"
#else
#include "hdf_nompi.ipp"
#endif

using M = MeshStructured<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using Expr = typename M::Expr;
template class Hdf<M>;

template void Hdf<M>::Read(FieldCell<Scal>&, std::string, M&, std::string);
template void Hdf<M>::Read(FieldCell<Vect>&, std::string, M&, std::string);
template void Hdf<M>::Read(FieldCell<Expr>&, std::string, M&, std::string);

template void Hdf<M>::Write(
    const FieldCell<Scal>&, std::string, M&, std::string);
template void Hdf<M>::Write(
    const FieldCell<Vect>&, std::string, M&, std::string);
template void Hdf<M>::Write(
    const FieldCell<Expr>&, std::string, M&, std::string);
