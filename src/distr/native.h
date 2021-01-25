// Created by Petr Karnakov on 02.01.2021
// Copyright 2021 ETH Zurich

#pragma once

#include <mpi.h>
#include <memory>

#include "distr.h"
#include "parse/vars.h"

template <class M>
std::unique_ptr<DistrMesh<M>> CreateNative(
    MPI_Comm, const KernelMeshFactory<M>&, Vars&);
