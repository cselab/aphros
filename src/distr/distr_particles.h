// Created by Petr Karnakov on 28.08.2021
// Copyright 2021 ETH Zurich

#pragma once

#include <functional>
#include <memory>
#include <vector>

#include "kernel/kernelmesh.h"
#include "util/mpi.h"

template <class M>
void TransferParticles(
    const std::vector<std::unique_ptr<KernelMesh<M>>>&, MPI_Comm,
    const std::function<int(int)>&);
