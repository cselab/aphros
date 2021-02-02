// Created by Petr Karnakov on 01.02.2021
// Copyright 2021 ETH Zurich

#include "raw.ipp"

namespace dump {

using M = MeshStructured<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using Expr = typename M::Expr;
template class Raw<M>;

template void Raw<M>::Read(FieldCell<Scal>&, const Meta& meta, std::string, M&);

template void Raw<M>::Write(
    const FieldCell<Scal>&, const Meta& meta, std::string, M&);

} // namespace dump
