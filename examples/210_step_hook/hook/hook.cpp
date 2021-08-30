// Created by Petr Karnakov on 29.08.2021
// Copyright 2021 ETH Zurich

#include <iostream>

#include <kernel/hydro.h>
#include <util/format.h>
#include <util/posthook.h>

static void* lmp;
static long natoms;

template <class M>
void InitHook(Hydro<M>* hydro) {
  auto& m = hydro->m;
  if (m.IsRoot()) {
    std::cout << util::Format("InitHook t={:}\n\n", hydro->fs_->GetTime());
  }
}

template <class M>
void StepHook(Hydro<M>* hydro) {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;

  auto& m = hydro->m;
  auto& var = hydro->var;
  auto sem = m.GetSem();

  // State persistent across all stages (sections `sem()`) of one call
  struct {
    Vect velocity_mean{0};
    Scal volume = 0;
  } * ctx(sem);
  auto& t = *ctx;

  if (sem()) {
    auto& fc = hydro->fs_->GetVelocity();
    for (auto c : m.CellsM()) {
      t.velocity_mean += fc[c] * c.volume();
      t.volume += c.volume();
    }
    for (auto d : m.dirs) {
      m.Reduce(&t.velocity_mean[d], Reduction::sum);
    }
    m.Reduce(&t.volume, Reduction::sum);
  }
  if (sem()) {
    t.velocity_mean /= t.volume;
    if (m.IsRoot()) {
      std::cout << util::Format(
          "StepHook: t={:} uservar={:} velocity_mean={:}\n\n",
          hydro->fs_->GetTime(), var.Int["uservar"], t.velocity_mean);
    }
  }
}

template <class M>
void FinalHook(Hydro<M>* hydro) {
  auto& m = hydro->m;
  if (m.IsRoot()) {
    std::cout << util::Format("FinalHook t={:}\n\n", hydro->fs_->GetTime());
  }
}

using M = MeshCartesian<double, 3>;
template void InitHook(Hydro<M>*);
template void StepHook(Hydro<M>*);
template void FinalHook(Hydro<M>*);
