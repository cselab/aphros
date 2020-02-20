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
constexpr size_t dim = M::dim;

template struct ULinear<Scal>;
template struct UEmbed<M>;

template FieldCell<T> UEB::Interpolate(const FieldEmbed<T>& feu, const EB& eb);

template FieldEmbed<T> UEB::Interpolate(
    const FieldCell<T>& fcu, const MapCondFace& mfc, size_t bc, T bcv,
    const EB& eb);
template FieldEmbed<TV> UEB::Interpolate(
    const FieldCell<TV>& fcu, const MapCondFace& mfc, size_t bc, TV bcv,
    const EB& eb);

template FieldEmbed<T> UEB::InterpolateUpwind(
    const FieldCell<T>& fcu, const FieldEmbed<Scal>& fev,
    const MapCondFace& mfc, size_t bc, T bcv, const EB& eb);

template FieldCell<T> UEB::AverageCutCells(
    const FieldCell<T>& fcu, const EB& eb);
template FieldCell<TV> UEB::AverageCutCells(
    const FieldCell<TV>& fcu, const EB& eb);

template FieldCell<T> UEB::RedistributeCutCells(
    const FieldCell<T>& fcu, const EB& eb);

template FieldEmbed<T> UEB::Gradient(
    const FieldCell<T>& fcu, const MapCondFace& mfc, size_t bc, T bcv,
    const EB& eb);

template FieldFace<T> UEB::InterpolateBilinearFaces(
    const FieldFace<T>& ffu, const EB& eb);

template FieldFace<T> UEB::InterpolateBilinear(
    const FieldCell<T>& fcu, const EB& eb);

template void UEB::InterpolateEmbedFaces(
    const FieldCell<T>& fcu, size_t bc, const MapCell<T>& mcu,
    FieldEmbed<T>& feu, const EB& eb);

template void UEB::InterpolateUpwindEmbedFaces(
    const FieldCell<T>& fcu, size_t bc, const MapCell<Scal>& mcu,
    const FieldEmbed<Scal>& fev, FieldEmbed<T>& feu, const EB& eb);

template FieldEmbed<T> UEB::InterpolateBilinear(
    const FieldCell<T>& fcu, size_t bc, const MapCell<T>& mcu, const EB& eb);
template FieldEmbed<TV> UEB::InterpolateBilinear(
    const FieldCell<TV>& fcu, size_t bc, const MapCell<TV>& mcu, const EB& eb);

template FieldEmbed<T> UEB::InterpolateBilinear(
    const FieldCell<T>& fcu, size_t bc, T bcv, const EB& eb);
template FieldEmbed<TV> UEB::InterpolateBilinear(
    const FieldCell<TV>& fcu, size_t bc, TV bcv, const EB& eb);

template FieldEmbed<T> UEB::InterpolateBilinear(
    const FieldCell<T>& fcu, const MapCondFace& mfc, size_t bc, T bcv,
    const EB& eb);
template FieldEmbed<TV> UEB::InterpolateBilinear(
    const FieldCell<TV>& fcu, const MapCondFace& mfc, size_t bc, TV bcv,
    const EB& eb);

template FieldFace<T> UEB::GradientBilinear(
    const FieldCell<T>& fcu, const EB& eb);

template FieldEmbed<T> UEB::GradientBilinear(
    const FieldCell<T>& fcu, size_t bc, const MapCell<T>& mcu, const EB& eb);

template FieldEmbed<T> UEB::GradientBilinear(
    const FieldCell<T>& fcu, const MapCondFace& mfc, size_t bc, T bcv,
    const EB& eb);
template FieldEmbed<TV> UEB::GradientBilinear(
    const FieldCell<TV>& fcu, const MapCondFace& mfc, size_t bc, TV bcv,
    const EB& eb);
