// Created by Petr Karnakov on 02.01.2021
// Copyright 2021 ETH Zurich

#include <set>

#include "util/format.h"

#include "comm_manager.h"

template <size_t dim_>
struct CommManager<dim_>::Imp {
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
    return (w + globalsize * 16) % globalsize;
  }
  // Returns list of cells to receive.
  // blocks: blocks on current rank
  // halos: number of layers of halo cells
  //        (1 for 3x3x3 stencil, 2 for 5x5x5 stencil)
  // globalsize: global mesh size
  // is_periodic: components are `true` for periodic directions
  static std::vector<LocalCell> GetRecvCells(
      const std::vector<Block>& blocks, int halos, bool full, MIdx globalsize,
      generic::Vect<bool, dim> is_periodic) {
    std::vector<LocalCell> res;

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
        res.push_back({ib, c});
      }
    }
    return res;
  }
  static std::vector<LocalCell> GetSendCells(
      const std::vector<Block>& blocks, MIdx globalsize,
      std::vector<LocalCell>& recv) {
    std::vector<MIdx> recv_midx; // multi-indices of inner cells to receive
    for (auto bc : recv) {
      const auto w = blocks[bc.block].indexc->GetMIdx(bc.cell);
      recv_midx.push_back(GetPeriodic(w, globalsize));
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
      fassert(cnt != 0, util::Format("Cell {} not found", w));
      fassert(cnt == 1, util::Format("Cell {} belongs to {} blocks", w, cnt));
      return ib_last;
    };

    std::vector<LocalCell> res;
    for (auto w : recv_midx) {
      auto ib = cell_to_block(w);
      auto c = blocks[ib].indexc->GetIdx(w);
      res.push_back({ib, c});
    }
    return res;
  }
};

// blocks: blocks owned by current rank
template <size_t dim_>
auto CommManager<dim_>::GetTasks(
    const std::vector<Block>& blocks, std::function<int(MIdx)>, MIdx globalsize,
    generic::Vect<bool, dim> is_periodic, const MpiWrapper&) -> Tasks {
  const int rank = 0;
  Tasks res;
  res.full_two.recv[rank] =
      Imp::GetRecvCells(blocks, 2, true, globalsize, is_periodic);
  res.direct_two.recv[rank] =
      Imp::GetRecvCells(blocks, 2, false, globalsize, is_periodic);
  res.full_one.recv[rank] =
      Imp::GetRecvCells(blocks, 1, true, globalsize, is_periodic);
  res.direct_one.recv[rank] =
      Imp::GetRecvCells(blocks, 1, false, globalsize, is_periodic);

  for (auto t : {
           &res.full_two,
           &res.direct_two,
           &res.full_one,
           &res.direct_one,
       }) {
    t->send[rank] = Imp::GetSendCells(blocks, globalsize, t->recv[rank]);
  }

  return res;
}
