// Created by Petr Karnakov on 27.10.2020
// Copyright 2020 ETH Zurich

#include "mesh.ipp"

using M = MeshCartesian<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;

template class MeshCartesian<double, 3>;

template M InitUniformMesh(
    Rect<typename M::Vect> domain, typename M::MIdx begin, typename M::MIdx s,
    int halos, bool isroot, bool islead, typename M::MIdx gs, int id);

template void M::ApplyNanFaces(FieldCell<Scal>& fc);
template void M::ApplyNanFaces(FieldCell<Vect>& fc);
