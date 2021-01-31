// Created by Petr Karnakov on 29.12.2019
// Copyright 2019 ETH Zurich

#include "approx.ipp"
#include "linear/linear.h"

using M = MeshStructured<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using T = Scal;
using TV = Vect;
constexpr size_t dim = 3;

template std::array<Scal, 3> GetCoeff(ConvSc);

template void SmoothenNode(FieldCell<T>& fc, M& m, size_t rep);
template void SmoothenNode(FieldNode<T>& fn, M& m, size_t iters);

template std::vector<Scal> GetGradCoeffs(Scal x, const std::vector<Scal>& z);

template std::vector<Scal> GetGradCoeffs(
    Scal x, const std::vector<Scal>& z, size_t b);

template void BcApply(
    FieldCell<T>& uc, const MapEmbed<BCond<T>>& me, const M& m);
template void BcApply(
    FieldCell<TV>& uc, const MapEmbed<BCond<TV>>& me, const M& m);

template void BcReflectAll(
    FieldCell<T>& uc, const MapEmbed<BCond<T>>& me, const M& m);
