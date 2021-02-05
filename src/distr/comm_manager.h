// Created by Petr Karnakov on 02.01.2021
// Copyright 2021 ETH Zurich

#pragma once

#include "distr/distrsolver.h"
#include "util/distr.h"

template <size_t dim_>
struct CommManager {
 public:
  static constexpr size_t dim = dim_;
  using MIdx = generic::MIdx<dim>;

  struct Block {
    const GBlockCells<dim>* incells; // inner cells
    // cell indexer, must be valid for cells in 5x5x5 stencil
    const GIndex<IdxCell, dim>* indexc;
  };
  struct LocalCell {
    size_t block;
    IdxCell cell;
  };
  // Communication task.
  struct Task {
    // send[rank] is the list of cells to send
    std::map<int, std::vector<LocalCell>> send;
    // recv[rank] is the list of cells to receive
    std::map<int, std::vector<LocalCell>> recv;
  };
  struct Tasks {
    // Communication tasks for various sets of neighbors
    Task full_one; // cells in 3x3x3 stencil
    Task full_two; // cells in 5x5x5 stencil
    Task direct_one; // cells in 3-point stencil in each direction
    Task direct_two; // cells in 5-point stencil in each direction
  };

  // blocks: blocks owned by current rank
  static Tasks GetTasks(
      const std::vector<Block>& blocks, std::function<int(MIdx)> cell_to_rank,
      MIdx globalsize, generic::Vect<bool, dim> is_periodic,
      const MpiWrapper& mpi);

  struct Imp;
};
