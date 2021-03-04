// Created by Petr Karnakov on 11.11.2020
// Copyright 2020 ETH Zurich

#include <omp.h>
#include <array>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "distr/comm_manager.h"
#include "parse/argparse.h"
#include "util/distr.h"

template <size_t dim>
int Run(const Vars& args, const MpiWrapper& mpi) {
  using MIdx = generic::MIdx<dim>;
  using Dir = GDir<dim>;
  MIdx meshsize;
  generic::Vect<bool, dim> is_periodic;
  for (size_t i = 0; i < dim; ++i) {
    meshsize[i] = args.Int[std::string() + 'n' + Dir(i).letter()];
    is_periodic[i] = args.Int[std::string("per") + Dir(i).letter()];
  }
  const int block = args.Int["block"];

  Subdomains<MIdx> sub(meshsize, MIdx(block), mpi.GetCommSize());

  using Manager = CommManager<dim>;
  using Block = typename Manager::Block;
  std::vector<Block> blocks;
  std::vector<GBlockCells<dim>> v_incells;
  std::vector<GIndexCells<dim>> v_indexc;

  GIndex<int, dim> index_proc(sub.info.procs);
  GIndex<int, dim> index_block(sub.info.blocks);
  const MIdx wproc = index_proc.GetMIdx(mpi.GetCommRank());
  for (auto wblock : GBlock<size_t, dim>(sub.info.blocks)) {
    const MIdx origin =
        (wproc * sub.info.blocks + wblock) * sub.info.block_size;
    v_incells.emplace_back(origin, sub.info.block_size);
    v_indexc.emplace_back(origin + MIdx(-2), sub.info.block_size + MIdx(4));
  }

  for (size_t i = 0; i < v_incells.size(); ++i) {
    blocks.push_back({&v_incells[i], &v_indexc[i]});
  }

  auto cell_to_rank = [&](MIdx w) { //
    return index_proc.GetIdx(w / sub.info.block_size / sub.info.blocks);
  };

  auto tasks =
      Manager::GetTasks(blocks, cell_to_rank, meshsize, is_periodic, mpi);
  auto myrank = mpi.GetCommRank();
  auto print = [&](auto& task, std::string name) {
    std::ofstream out(util::Format("out_{}", mpi.GetCommRank()));
    out << name << '\n';
    for (auto w : GBlockCells<dim>(meshsize)) {
      out << w << ' ' << cell_to_rank(w) << '\n';
    }
    for (auto& p : task.recv) {
      out << myrank << ": recv from rank " << p.first << "\n";
      for (auto bc : p.second) {
        out << blocks[bc.block].indexc->GetMIdx(bc.cell) << ' ';
      }
      out << '\n';
    }
    for (auto& p : task.send) {
      out << myrank << ": send to rank " << p.first << "\n";
      for (auto bc : p.second) {
        out << blocks[bc.block].indexc->GetMIdx(bc.cell) << ' ';
      }
      out << '\n';
    }
  };

  print(tasks.full_one, "full_one");
  print(tasks.direct_one, "direct_one");
  print(tasks.full_two, "full_two");
  print(tasks.direct_two, "direct_two");

  return 0;
}

int main(int argc, const char** argv) {
  MpiWrapper mpi(&argc, &argv);

  ArgumentParser parser("Test for communication maps.", mpi.IsRoot());
  parser.AddVariable<int>("--dim", 2).Help("Dimension").Options({1, 2, 3, 4});
  parser.AddVariable<int>("--nx", 2).Help("Mesh size in x-direction");
  parser.AddVariable<int>("--ny", 2).Help("Mesh size in y-direction");
  parser.AddVariable<int>("--nz", 2).Help("Mesh size in z-direction");
  parser.AddVariable<int>("--nw", 2).Help("Mesh size in w-direction");
  parser.AddVariable<int>("--perx", 1).Help("Periodic in x-direction");
  parser.AddVariable<int>("--pery", 0).Help("Periodic in y-direction");
  parser.AddVariable<int>("--perz", 1).Help("Periodic in z-direction");
  parser.AddVariable<int>("--perw", 1).Help("Periodic in w-direction");
  parser.AddVariable<int>("--block", 1).Help("Block size").Options({1, 2, 8});
  auto args = parser.ParseArgs(argc, argv);
  if (const int* p = args.Int.Find("EXIT")) {
    return *p;
  }

  const int dim = args.Int["dim"];
  switch (dim) {
    case 1:
      return Run<1>(args, mpi);
    case 2:
      return Run<2>(args, mpi);
    case 3:
      return Run<3>(args, mpi);
    case 4:
      return Run<4>(args, mpi);
    default:
      fassert(false, util::Format("Unknown dim={}", dim));
  }
}
