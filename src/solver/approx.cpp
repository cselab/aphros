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

template void InterpolateI(
    const FieldCell<T>& fc, const FieldCell<typename M::Vect>& fcgp,
    const FieldFace<T>& ffw, const M& m, ConvSc sc, typename M::Scal th,
    FieldFace<T>& ff);

template void Interpolate(
    const FieldCell<T>& fc, const FieldCell<typename M::Vect>& fcgp,
    const MapCondFace& mfc, const FieldFace<T>& ffw, const M& m, ConvSc sc,
    typename M::Scal th, FieldFace<T>& ff);

template void GradientI(const FieldCell<T>& fc, const M& m, FieldFace<T>& ff);

template void GradientI(const FieldCell<TV>& fc, const M& m, FieldFace<TV>& ff);

template void GradientB(
    const FieldCell<T>& fc, const MapCondFace& mfc, const M& m,
    FieldFace<T>& ff);

template void GradientB(
    const FieldCell<TV>& fc, const MapCondFace& mfc, const M& m,
    FieldFace<TV>& ff);

template void Gradient(
    const FieldCell<T>& fc, const MapCondFace& mfc, const M& m,
    FieldFace<T>& ff);

template FieldFace<T> Interpolate(const FieldNode<T>& fn, const M& m);

template void InterpolateI(
    const FieldCell<T>& fc, FieldFace<T>& ff, const M& m);

template void InterpolateS(
    const FieldCell<T>& fc, FieldFace<T>& ff, const M& m);

template void InterpolateB(
    const FieldCell<T>& fc, const MapCondFace& mfc, FieldFace<T>& ff,
    const M& m);

template FieldFace<T> Interpolate(
    const FieldCell<T>& fc, const MapCondFace& mfc, const M& m);
template FieldFace<TV> Interpolate(
    const FieldCell<TV>& fc, const MapCondFace& mfc, const M& m);

template FieldFace<typename M::Scal> InterpolateSuperbee(
    const FieldCell<typename M::Scal>& fc,
    const FieldCell<typename M::Vect>& fcg, const MapCondFace& mfc,
    const FieldFace<typename M::Scal>& ffw, const M& m, typename M::Scal th);

template FieldCell<T> Average(const FieldFace<T>& ff, const M& m);
template FieldCell<TV> Average(const FieldFace<TV>& ff, const M& m);

template void Smoothen(
    FieldCell<T>& fc, const MapCondFace& mfc, M& m, size_t rep);
template void Smoothen(
    FieldCell<TV>& fc, const MapCondFace& mfc, M& m, size_t rep);
template void SmoothenNode(FieldCell<T>& fc, M& m, size_t rep);
template void SmoothenNode(FieldNode<T>& fn, M& m, size_t iters);

template FieldCell<typename M::Vect> Gradient(
    const FieldFace<typename M::Scal>& ff, const M& m);

template <class Field, class M>
Scal CalcDiff(const Field& fa, const Field& fb, const M& m);

template <class Idx, class M>
Scal CalcDiff(
    const GField<typename M::Vect, Idx>& fa,
    const GField<typename M::Vect, Idx>& fb, const M& m);

template std::vector<Scal> GetGradCoeffs(Scal x, const std::vector<Scal>& z);

template std::vector<Scal> GetGradCoeffs(
    Scal x, const std::vector<Scal>& z, size_t b);

template void BcApply(FieldCell<T>& uc, const MapCondFace& mfc, const M& m);
template void BcApply(FieldCell<TV>& uc, const MapCondFace& mfc, const M& m);
template void BcApply(
    FieldCell<T>& uc, const MapEmbed<BCond<T>>& me, const M& m);
template void BcApply(
    FieldCell<TV>& uc, const MapEmbed<BCond<TV>>& me, const M& m);

template void BcReflectAll(
    FieldCell<T>& uc, const MapCondFace& mfc, const M& m);
template void BcReflectAll(
    FieldCell<T>& uc, const MapEmbed<BCond<T>>& me, const M& m);
