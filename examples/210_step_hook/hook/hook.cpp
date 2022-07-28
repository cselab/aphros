// Created by Petr Karnakov on 29.08.2021
// Copyright 2021 ETH Zurich

#include <iostream>

#include <kernel/hydro.h>
#include <util/format.h>
#include <util/posthook.h>

// Shared context, one for each processor.
struct HookContextShared {
  int foo;
};

// Local context, one for each block.
struct HookContext {
  int bar;
};


template <class M>
void InitHook(Hydro<M>* hydro) {
  auto& m = hydro->m;
  auto sem = m.GetSem();
  if (sem()) {
    if (m.IsLead()) { // Executed only on the lead block, once in each rank.
                      // There is only one lead block in each rank.
      // Create shared context.
      hydro->par_.hook_context_shared =
          std::make_unique<Holder<HookContextShared>>(new HookContextShared());
      auto* hsctx = dynamic_cast<Holder<HookContextShared>*>(
                        hydro->par_.hook_context_shared.get())
                        ->Get();
      hsctx->foo = 1000 * (m.GetId() + 1);
    }
    if (m.IsRoot()) { // Executed only on the root block, once over all ranks.
                      // There is only one root block over all ranks.
      auto* hsctx = dynamic_cast<Holder<HookContextShared>*>(
                        hydro->par_.hook_context_shared.get())
                        ->Get();
      std::cout << util::Format(
          "InitHook t={:} foo={:}\n\n", hydro->fs_->GetTime(), hsctx->foo);
    }
    // Create local context.
    hydro->hook_context_ =
        std::make_unique<Holder<HookContext>>(new HookContext());
    auto* hctx =
        dynamic_cast<Holder<HookContext>*>(hydro->hook_context_.get())->Get();
    hctx->bar = 1000 * (m.GetId() + 1);
  }
  if (sem()) {
    auto* hsctx = dynamic_cast<Holder<HookContextShared>*>(
                      hydro->par_.hook_context_shared.get())
                      ->Get();
    std::cout << util::Format(
        "InitHook t={:} foo={:}\n\n", hydro->fs_->GetTime(), hsctx->foo);
  }
}

template <class M>
void StepHook(Hydro<M>* hydro) {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;

  auto& m = hydro->m;
  auto& var = hydro->var;
  auto sem = m.GetSem();
  auto* hsctx = dynamic_cast<Holder<HookContextShared>*>(
                    hydro->par_.hook_context_shared.get())
                    ->Get();
  auto* hctx =
      dynamic_cast<Holder<HookContext>*>(hydro->hook_context_.get())->Get();

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
    hsctx->foo += 1; // Increment once in each block.
    hctx->bar += 1; // Increment once in each block.
  }
  if (sem()) {
    t.velocity_mean /= t.volume;
    if (m.IsRoot()) {
      std::cout << util::Format(
          "StepHook: t={:} uservar={:} velocity_mean={:}\n",
          hydro->fs_->GetTime(), var.Int["uservar"], t.velocity_mean);
    }
    std::cout << util::Format(
        "HookContext: id={:} foo={:} bar={:}\n", m.GetId(), hsctx->foo,
        hctx->bar);
  }
}

template <class M>
void FinalHook(Hydro<M>* hydro) {
  auto& m = hydro->m;
  auto sem = m.GetSem();
  auto* hsctx = dynamic_cast<Holder<HookContextShared>*>(
                    hydro->par_.hook_context_shared.get())
                    ->Get();
  auto* hctx =
      dynamic_cast<Holder<HookContext>*>(hydro->hook_context_.get())->Get();
  if (sem()) {
    if (m.IsRoot()) {
      std::cout << util::Format(
          "FinalHook t={:} foo={:} bar={:}\n\n", hydro->fs_->GetTime(), hsctx->foo,
          hctx->bar);
    }
  }
}

using M = MeshCartesian<double, 3>;
template void InitHook(Hydro<M>*);
template void StepHook(Hydro<M>*);
template void FinalHook(Hydro<M>*);
