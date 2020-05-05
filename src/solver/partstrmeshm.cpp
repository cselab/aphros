// Created by Petr Karnakov on 28.08.2019
// Copyright 2019 ETH Zurich

#include "partstrmeshm.ipp"
#include "embed.h"

using M = MeshStructured<double, 3>;
template class PartStrMeshM<M>;

template void PartStrMeshM<M>::Part(
    const Multi<const FieldCell<Scal>*>& vfca,
    const Multi<const FieldCell<Vect>*>& vfcn,
    const Multi<const FieldCell<bool>*>& vfci,
    const Multi<const FieldCell<Scal>*>& vfccl,
    const FieldCell<Scal>* fc_contang, const Embed<M>& eb);

template void PartStrMeshM<M>::Part(
    const Multi<const FieldCell<Scal>*>& vfca,
    const Multi<const FieldCell<Vect>*>& vfcn,
    const Multi<const FieldCell<bool>*>& vfci,
    const Multi<const FieldCell<Scal>*>& vfccl,
    const FieldCell<Scal>* fc_contang, const M& eb);
