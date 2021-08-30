// Created by Petr Karnakov on 29.08.2021
// Copyright 2021 ETH Zurich

#include <iostream>

#include <kernel/hydro.h>
#include <util/format.h>
#include <util/posthook.h>

// Data to be shared among all local blocks
struct SharedState {
  int foo;
};

template <class M>
void InitHook(Hydro<M>* hydro) {
  auto& m = hydro->m;
  auto sem = m.GetSem();
  if (sem()) {
    if (m.IsLead()) { // Executed only on the lead block, once in each rank.
                      // There is only one lead block in each rank.
      // Allocate memory for the shared state and save pointer in the solver.
      hydro->par_.ptr = new SharedState();
      auto* shared = static_cast<SharedState*>(hydro->par_.ptr);
      shared->foo = 1000;
    }
    if (m.IsRoot()) { // Executed only on the root block, once over all ranks.
                      // There is only one root block over all ranks.
      auto* shared = static_cast<SharedState*>(hydro->par_.ptr);
      std::cout << util::Format(
          "InitHook t={:} foo={:}\n\n", hydro->fs_->GetTime(), shared->foo);
    }
  }
}

template <class M>
void StepHook(Hydro<M>* hydro) {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;

  auto& m = hydro->m;
  auto& var = hydro->var;
  auto sem = m.GetSem();
  auto* shared = static_cast<SharedState*>(hydro->par_.ptr);

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
    shared->foo += 1; // Increment once in each block.
  }
  if (sem()) {
    t.velocity_mean /= t.volume;
    if (m.IsRoot()) {
      std::cout << util::Format(
          "StepHook: t={:} uservar={:} foo={:} velocity_mean={:}\n\n",
          hydro->fs_->GetTime(), var.Int["uservar"], shared->foo,
          t.velocity_mean);
    }
  }
}

template <class M>
void FinalHook(Hydro<M>* hydro) {
  auto& m = hydro->m;
  auto sem = m.GetSem();
  auto* shared = static_cast<SharedState*>(hydro->par_.ptr);
  if (sem()) {
    if (m.IsRoot()) {
      std::cout << util::Format(
          "FinalHook t={:} foo={:}\n\n", hydro->fs_->GetTime(), shared->foo);
    }
    if (m.IsLead()) {
      delete shared;
    }
  }
}

using M = MeshCartesian<double, 3>;
template void InitHook(Hydro<M>*);
template void StepHook(Hydro<M>*);
template void FinalHook(Hydro<M>*);
