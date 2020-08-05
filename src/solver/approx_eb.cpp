// Created by Petr Karnakov on 11.02.2020
// Copyright 2020 ETH Zurich

#include "approx_eb.ipp"
#include "linear/linear.h"

using M = MeshStructured<double, 3>;
using EB = Embed<M>;
using UEB = UEmbed<M>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using T = Scal;
using TV = Vect;
using Expr = typename M::Expr;
constexpr size_t dim = M::dim;

template struct ULinear<Scal>;
template struct UEmbed<M>;

template FieldCell<T> UEB::Interpolate(const FieldEmbed<T>& feu, const EB& eb);
template FieldCell<T> UEB::Interpolate(const FieldFace<T>& feu, const M&);

template FieldCell<T> UEB::AverageCutCells(
    const FieldCell<T>& fcu, const EB& eb);
template FieldCell<TV> UEB::AverageCutCells(
    const FieldCell<TV>& fcu, const EB& eb);

template FieldCell<T> UEB::RedistributeCutCells(
    const FieldCell<T>& fcu, const EB& eb);
template FieldCell<T> UEB::RedistributeCutCells(
    const FieldCell<T>& fcu, const M& m);

template FieldFace<T> UEB::InterpolateBilinearFaces(
    const FieldFace<T>& ffu, const EB& eb);
template FieldFace<T> UEB::InterpolateBilinearFaces(
    const FieldFace<T>& ffu, const M& m);

template FieldEmbed<T> UEB::Interpolate(
    const FieldCell<T>& fcu, const MapEmbed<BCond<T>>& mebc, const EB& eb);
template FieldEmbed<TV> UEB::Interpolate(
    const FieldCell<TV>& fcu, const MapEmbed<BCond<TV>>& mebc, const EB& eb);
template FieldEmbed<T> UEB::Gradient(
    const FieldCell<T>& fcu, const MapEmbed<BCond<T>>& mebc, const EB& eb);
template FieldEmbed<TV> UEB::Gradient(
    const FieldCell<TV>& fcu, const MapEmbed<BCond<TV>>& mebc, const EB& eb);

template FieldFace<T> UEB::Interpolate(
    const FieldCell<T>& fcu, const MapEmbed<BCond<T>>& mebc, const M& m);
template FieldFace<TV> UEB::Interpolate(
    const FieldCell<TV>& fcu, const MapEmbed<BCond<TV>>& mebc, const M& m);
template FieldFace<T> UEB::Gradient(
    const FieldCell<T>& fcu, const MapEmbed<BCond<T>>& mebc, const M& m);
template FieldFace<TV> UEB::Gradient(
    const FieldCell<TV>& fcu, const MapEmbed<BCond<TV>>& mebc, const M& m);
