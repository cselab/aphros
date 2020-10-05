// Created by Petr Karnakov on 28.06.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <limits>

#include "solver/solver.h"
#include "solver/trackerm.h"
#include "util/vof.h"

template <class M_>
class Tracker {
  using M = M_;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using MIdx = typename M::MIdx;
  using TRM = Trackerm<M>;
  static constexpr size_t dim = M::dim;

 public:
  // fccl: initial color [a]
  Tracker(M& m_, const FieldCell<Scal>& fccl) : m(m_), layers(1), fccl_(fccl) {
    trm_ = std::unique_ptr<TRM>(new TRM(m, layers));
  }
  // Propagates color.
  // fcu: volume fraction [a]
  void Update(
      const FieldCell<Scal>& fcu, Scal th, Scal clfixed, Vect clfixed_x,
      const MapCondFace& mfc, bool unionfind, bool reduce, bool grid);
  // Returns color [a]
  const FieldCell<Scal>& GetColor() const {
    return fccl_;
  }
  MIdx GetImage(IdxCell c) const {
    return trm_->GetImage(0, c);
  }
  static constexpr Scal kClNone = -1.; // no color

 private:
  M& m;
  GRange<size_t> layers;
  FieldCell<Scal> fccl_; // color [a]
  size_t dm_;
  UVof<M> uvof_;
  std::unique_ptr<TRM> trm_;
};

template <class M_>
void Tracker<M_>::Update(
    const FieldCell<Scal>& fcu, Scal th, Scal clfixed, Vect clfixed_x,
    const MapCondFace& mfc, bool unionfind, bool reduce, bool grid) {
  auto sem = m.GetSem("upd");

  struct {
    FieldCell<Scal> fccl0;
  } * ctx(sem);
  auto& fccl0 = ctx->fccl0;

  if (sem("")) {
    fccl0 = fccl_;
    for (auto c : m.AllCells()) {
      fccl_[c] = (fcu[c] > th ? 0 : kClNone);
    }
  }
  if (sem.Nested()) {
    uvof_.Recolor(
        layers, &fcu, &fccl_, &fccl0, clfixed, clfixed_x, 1e10, mfc, false,
        unionfind, reduce, grid, m);
  }
  if (sem.Nested()) {
    trm_->Update(&fccl_, &fccl0);
  }
}
