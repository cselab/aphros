// Created by Petr Karnakov on 27.10.2020
// Copyright 2020 ETH Zurich

#include "mesh.ipp"

template class MeshStructured<double, 3>;

template M InitUniformMesh(
    Rect<typename M::Vect> domain, typename M::MIdx begin, typename M::MIdx s,
    int halos, bool isroot, bool islead, typename M::MIdx gs, int id);
