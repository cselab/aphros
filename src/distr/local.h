// Created by Petr Karnakov on 25.04.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <mpi.h>
#include <memory>

#include "distr.h"
#include "parse/vars.h"

template <class M>
std::unique_ptr<DistrMesh<M>> CreateLocal(
    MPI_Comm, const KernelMeshFactory<M>&, Vars&);
