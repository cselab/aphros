// Created by Petr Karnakov on 11.11.2020
// Copyright 2020 ETH Zurich

#include <omp.h>
#include <algorithm>
#include <iostream>
#include <map>
#include <set>

#include "distr/commmap.h"
#include "distr/distrbasic.h"
#include "parse/argparse.h"
#include "util/distr.h"

// Communication mapper.
// Assigns blocks to processors and enumerates halo cells.
// No dependency on MPI.
// Valid only on root.
// Supports irregular domains.
// For example, to exclude blocks with all excluded cells.
//
// The computational domain is a subset of cells.
// Its bounding box is divided into equally sizes blocks.
// Blocks that contain cells from the computational domain are
// distributed among multiple processors.
// Each block is assigned to one processor,
// while one processor may have multiple blocks.
//
// Each block is characterized by:
// * extent in cell space
// * processor
// * cost
//
// The processors are assigned to blocks such that:
// * the total cost of blocks owned by a processor is uniform
// * the communication area is minimized
// * regular enumeration is a special case if all blocks have the same cost
//   and the number of blocks is divisible by the number of processors
//   in each direction

// TODO move to idx.h
template <>
struct std::less<IdxCell> {
  bool operator()(IdxCell ca, IdxCell cb) const {
    return ca.raw() < cb.raw();
  }
};

template <size_t dim_>
struct CommMapper {
 public:
  static constexpr size_t dim = dim_;
  using MIdx = generic::MIdx<dim>;

  struct Block {
    const GBlockCells<dim>* incells; // inner cells
    // cell indexer, must be valid for cells in 5x5x5 stencil
    const GIndex<IdxCell, dim>* indexc;
  };

  struct BlockCell {
    size_t block;
    IdxCell cell;
  };

  // Communication task.
  struct Task {
    // send[rank] is the list of cells to send
    std::map<int, std::vector<BlockCell>> send;
    // recv[rank] is the list of cells to receive
    std::map<int, std::vector<BlockCell>> recv;
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
      MIdx globalsize, std::array<bool, dim> is_periodic,
      const MpiWrapper& mpi) {
    Tasks res;
    res.full_two.recv =
        GetRecvCells(blocks, cell_to_rank, 2, true, globalsize, is_periodic);
    res.direct_two.recv =
        GetRecvCells(blocks, cell_to_rank, 2, false, globalsize, is_periodic);
    res.full_one.recv =
        GetRecvCells(blocks, cell_to_rank, 1, true, globalsize, is_periodic);
    res.direct_one.recv =
        GetRecvCells(blocks, cell_to_rank, 1, false, globalsize, is_periodic);

    for (auto t : {
             &res.full_two,
             &res.direct_two,
             &res.full_one,
             &res.direct_one,
         }) {
      t->send = GetSendCells(blocks, mpi, globalsize, t->recv);
    }

    return res;
  }

 private:
  static bool IsValid(
      MIdx w, MIdx globalsize, std::array<bool, dim> is_periodic) {
    for (size_t i = 0; i < dim; ++i) {
      if (!is_periodic[i] && (w[i] < 0 || w[i] >= globalsize[i])) {
        return false;
      }
    }
    return true;
  }
  static MIdx GetPeriodic(MIdx w, MIdx globalsize) {
    return (w + globalsize * 16) % globalsize;
  }
  static std::map<int, std::vector<BlockCell>> GetRecvCells(
      const std::vector<Block>& blocks, std::function<int(MIdx)> cell_to_rank,
      int halos, bool full, MIdx globalsize,
      std::array<bool, dim> is_periodic) {
    std::map<int, std::vector<BlockCell>> res;

    auto nnz = [](MIdx dw) {
      int cnt = 0;
      for (size_t i = 0; i < dw.dim; ++i) {
        cnt += (dw[i] != 0);
      }
      return cnt;
    };

    for (size_t ib = 0; ib < blocks.size(); ++ib) {
      auto& b = blocks[ib];
      for (auto wc : *b.incells) {
        GBlock<size_t, dim> offsets(MIdx(-halos), MIdx(2 * halos + 1));
        for (auto dw : offsets) {
          if (full || nnz(dw) == 1) {
            const auto w = wc + dw;
            const auto wp = GetPeriodic(w, globalsize);
            if (IsValid(w, globalsize, is_periodic) &&
                !b.incells->IsInside(w)) {
              res[cell_to_rank(wp)].push_back({ib, b.indexc->GetIdx(w)});
            }
          }
        }
      }
    }
    return res;
  }
  static std::map<int, std::vector<BlockCell>> GetSendCells(
      const std::vector<Block>& blocks, const MpiWrapper& mpi, MIdx globalsize,
      const std::map<int, std::vector<BlockCell>>& recv) {
    const int myrank = mpi.GetCommRank();
    std::vector<int> msg_count(mpi.GetCommSize(), 0);
    for (auto& p : recv) {
      auto& rank = p.first;
      ++msg_count[rank];
    }
    MPI_Allreduce(
        MPI_IN_PLACE, msg_count.data(), msg_count.size(), MPI_INT, MPI_SUM,
        mpi.GetComm());

    // msg_count[myrank] is the number of messages to receive.

    // multi-indices of inner cells to receive
    std::map<int, std::vector<MIdx>> rank_to_midx;
    for (auto& p : recv) {
      auto& rank = p.first;
      auto& midx = rank_to_midx[rank];
      for (auto bc : p.second) {
        const auto w = blocks[bc.block].indexc->GetMIdx(bc.cell);
        midx.push_back(GetPeriodic(w, globalsize));
      }
    }

    // Send cell multi-indices.
    std::vector<MPI_Request> reqs;
    reqs.reserve(rank_to_midx.size());
    const int tag = 0;
    for (auto& p : rank_to_midx) {
      auto& rank = p.first;
      auto& midx = p.second;
      reqs.emplace_back();
      MPI_Isend(
          midx.data(), midx.size() * sizeof(MIdx), MPI_CHAR, rank, tag,
          mpi.GetComm(), &reqs.back());
    }

    auto cell_to_block = [&](MIdx w) -> size_t {
      size_t cnt = 0;
      size_t ib_last = 0;
      for (size_t ib = 0; ib < blocks.size(); ++ib) {
        auto& b = blocks[ib];
        if (b.incells->IsInside(w)) {
          ++cnt;
          ib_last = ib;
        }
      }
      fassert(
          cnt != 0, util::Format("Cell {} not owned by rank {}", w, myrank));
      fassert(
          cnt == 1,
          util::Format(
              "Cell {} belongs to {} blocks on rank {}", w, cnt, myrank));
      return ib_last;
    };

    std::map<int, std::vector<BlockCell>> res;

    // Receive cell multi-indices.
    for (int i = 0; i < msg_count[myrank]; ++i) {
      MPI_Status status;
      MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, mpi.GetComm(), &status);
      int count;
      MPI_Get_count(&status, MPI_CHAR, &count);
      fassert_equal(count % sizeof(MIdx), 0);
      const int rank = status.MPI_SOURCE;
      fassert(
          res[rank].empty(),
          util::Format("Got more than one message from rank {}", rank));
      // indices of cells required by remote rank
      std::vector<MIdx> midx(count / sizeof(MIdx));
      MPI_Recv(
          midx.data(), midx.size() * sizeof(MIdx), MPI_CHAR, rank, tag,
          mpi.GetComm(), MPI_STATUS_IGNORE);
      for (auto w : midx) {
        auto ib = cell_to_block(w);
        auto c = blocks[ib].indexc->GetIdx(w);
        res[rank].push_back({ib, c});
      }
    }
    MPI_Waitall(reqs.size(), reqs.data(), MPI_STATUSES_IGNORE);
    return res;
  }

  GIndex<IdxCell, dim> index_;
  GBlock<IdxCell, dim> inner_;
};

template <size_t dim>
int Run(const Vars& args, const MpiWrapper& mpi) {
  using MIdx = generic::MIdx<dim>;
  using Dir = GDir<dim>;
  MIdx meshsize;
  std::array<bool, dim> is_periodic;
  for (size_t i = 0; i < dim; ++i) {
    meshsize[i] = args.Int[std::string() + 'n' + Dir(i).letter()];
    is_periodic[i] = args.Int[std::string("per") + Dir(i).letter()];
  }
  const int block = args.Int["block"];

  Subdomains<MIdx> sub(meshsize, MIdx(block), mpi.GetCommSize());

  using Mapper = CommMapper<dim>;
  using Block = typename Mapper::Block;
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
      Mapper::GetTasks(blocks, cell_to_rank, meshsize, is_periodic, mpi);
  auto myrank = mpi.GetCommRank();
  auto print = [&](auto& task, std::string name) {
    auto& out = std::cout;
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
  parser.AddVariable<int>("--dim", 3).Help("Dimension").Options({1, 2, 3, 4});
  parser.AddVariable<int>("--nx", 4).Help("Mesh size in x-direction");
  parser.AddVariable<int>("--ny", 4).Help("Mesh size in y-direction");
  parser.AddVariable<int>("--nz", 4).Help("Mesh size in z-direction");
  parser.AddVariable<int>("--nw", 4).Help("Mesh size in w-direction");
  parser.AddVariable<int>("--perx", 1).Help("Periodic in x-direction");
  parser.AddVariable<int>("--pery", 1).Help("Periodic in y-direction");
  parser.AddVariable<int>("--perz", 1).Help("Periodic in z-direction");
  parser.AddVariable<int>("--perw", 1).Help("Periodic in w-direction");
  parser.AddVariable<int>("--block", 2).Help("Block size").Options({1, 2, 8});
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
      throw std::runtime_error(util::Format("Unknown dim={}", dim));
  }
}
