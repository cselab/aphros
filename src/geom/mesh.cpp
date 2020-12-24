// Created by Petr Karnakov on 27.10.2020
// Copyright 2020 ETH Zurich

#include "mesh.ipp"

template class MeshCartesian<double, 3>;
using M = MeshCartesian<double, 3>;
template M InitUniformMesh(
    Rect<typename M::Vect> domain, typename M::MIdx begin, typename M::MIdx s,
    int halos, bool isroot, bool islead, typename M::MIdx gs, int id);
template void M::ApplyNanFaces(FieldCell<typename M::Scal>& fc);
template void M::ApplyNanFaces(FieldCell<typename M::Vect>& fc);


template class MeshCartesian<double, 4>;
using M4 = MeshCartesian<double, 4>;
template M4 InitUniformMesh(
    Rect<typename M4::Vect> domain, typename M4::MIdx begin, typename M4::MIdx s,
    int halos, bool isroot, bool islead, typename M4::MIdx gs, int id);
template void M4::ApplyNanFaces(FieldCell<typename M4::Scal>& fc);
template void M4::ApplyNanFaces(FieldCell<typename M4::Vect>& fc);
