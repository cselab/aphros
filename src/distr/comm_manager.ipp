// Created by Petr Karnakov on 02.01.2021
// Copyright 2021 ETH Zurich

#include <set>

#include "util/format.h"

#include "comm_manager.h"

template <size_t dim_>
struct CommManager<dim_>::Imp {
  using Rank = int;

  static bool IsValid(
      MIdx w, MIdx globalsize, generic::Vect<bool, dim> is_periodic) {
    for (size_t i = 0; i < dim; ++i) {
      if (!is_periodic[i] && (w[i] < 0 || w[i] >= globalsize[i])) {
        return false;
      }
    }
    return true;
  }
  static MIdx GetPeriodic(MIdx w, MIdx globalsize) {
    return w.mod_positive(globalsize);
  }
  // Returns list of cells to receive from each rank.
  // blocks: blocks on current rank
  // cell_to_rank: function that returns rank from global cell index
  // halos: number of layers of halo cells
  //        (1 for 3x3x3 stencil, 2 for 5x5x5 stencil)
  // globalsize: global mesh size
  // is_periodic: components are `true` for periodic directions
  static std::map<Rank, std::vector<LocalCell>> GetRecvCells(
      const std::vector<Block>& blocks, std::function<Rank(MIdx)> cell_to_rank,
      int halos, bool full, MIdx globalsize,
      generic::Vect<bool, dim> is_periodic) {
    std::map<Rank, std::vector<LocalCell>> res;

    // Returns the number of non-zero components
    auto nnz = [](MIdx dw) {
      int cnt = 0;
      for (size_t i = 0; i < MIdx::dim; ++i) {
        cnt += (dw[i] != 0);
      }
      return cnt;
    };

    for (size_t ib = 0; ib < blocks.size(); ++ib) {
      auto& b = blocks[ib];
      std::set<IdxCell> seen_cells;
      for (auto wc : *b.incells) {
        GBlock<size_t, dim> offsets(MIdx(-halos), MIdx(2 * halos + 1));
        for (auto dw : offsets) {
          if (full || nnz(dw) == 1) {
            const auto w = wc + dw;
            if (IsValid(w, globalsize, is_periodic) &&
                !b.incells->IsInside(w)) {
              seen_cells.insert(b.indexc->GetIdx(w));
            }
          }
        }
      }
      for (auto c : seen_cells) {
        const auto wp = GetPeriodic(b.indexc->GetMIdx(c), globalsize);
        res[cell_to_rank(wp)].push_back({ib, c});
      }
    }
    return res;
  }
  // Returns list of cells to send in each rank.
  // blocks: blocks on current rank
  // cell_to_rank: function that returns rank from global cell index
  // halos: number of layers of halo cells
  //        (1 for 3x3x3 stencil, 2 for 5x5x5 stencil)
  // globalsize: global mesh size
  // is_periodic: components are `true` for periodic directions
  static std::map<Rank, std::vector<LocalCell>> GetSendCells(
      const std::vector<Block>& blocks, const MpiWrapper& mpi, MIdx globalsize,
      const std::map<Rank, std::vector<LocalCell>>& recv) {
    const Rank myrank = mpi.GetCommRank();
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
    std::map<Rank, std::vector<MIdx>> rank_to_midx;
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

    std::map<Rank, std::vector<LocalCell>> res;

    // Receive cell multi-indices.
    for (int i = 0; i < msg_count[myrank]; ++i) {
      MPI_Status status;
      MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, mpi.GetComm(), &status);
      int count;
      MPI_Get_count(&status, MPI_CHAR, &count);
      fassert_equal(count % sizeof(MIdx), 0);
      const Rank rank = status.MPI_SOURCE;
      fassert(
          !res.count(rank),
          util::Format(
              "Rank {} got more than one message from rank {}. "
              "Expected new message of size {}, buffer already contains {} "
              "cells",
              myrank, rank, count / sizeof(MIdx), res[rank].size()));
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
    MPI_Barrier(mpi.GetComm());
    return res;
  }
};

// blocks: blocks owned by current rank
template <size_t dim_>
auto CommManager<dim_>::GetTasks(
    const std::vector<Block>& blocks, std::function<int(MIdx)> cell_to_rank,
    MIdx globalsize, generic::Vect<bool, dim> is_periodic,
    const MpiWrapper& mpi) -> Tasks {
  Tasks res;
  res.full_two.recv =
      Imp::GetRecvCells(blocks, cell_to_rank, 2, true, globalsize, is_periodic);
  res.direct_two.recv = Imp::GetRecvCells(
      blocks, cell_to_rank, 2, false, globalsize, is_periodic);
  res.full_one.recv =
      Imp::GetRecvCells(blocks, cell_to_rank, 1, true, globalsize, is_periodic);
  res.direct_one.recv = Imp::GetRecvCells(
      blocks, cell_to_rank, 1, false, globalsize, is_periodic);

  for (auto t : {
           &res.full_two,
           &res.direct_two,
           &res.full_one,
           &res.direct_one,
       }) {
    t->send = Imp::GetSendCells(blocks, mpi, globalsize, t->recv);
  }

  return res;
}
