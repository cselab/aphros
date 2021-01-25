// Created by Petr Karnakov on 25.11.2019
// Copyright 2019 ETH Zurich

#pragma once

#include <mpi.h>
#include <memory>

#include "distr.h"
#include "parse/vars.h"

template <class M>
std::unique_ptr<DistrMesh<M>> CreateCubismnc(
    MPI_Comm, const KernelMeshFactory<M>&, Vars&);
