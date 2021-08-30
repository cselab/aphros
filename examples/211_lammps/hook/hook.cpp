// Created by Petr Karnakov on 29.08.2021
// Copyright 2021 ETH Zurich

#include <iostream>

#include <kernel/hydro.h>
#include <library.h>
#include <util/format.h>
#include <util/posthook.h>

struct Lammps {
  void* lmp;
  long natoms;
  double* x;
};

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
    std::cout << util::Format(
        "StepHook: t={:} uservar={:} velocity_mean={:}\n\n",
        hydro->fs_->GetTime(), var.Int["uservar"], t.velocity_mean);
    Lammps* plmp = (Lammps*)hydro->par_.ptr;
    lammps_gather_atoms(plmp->lmp, (char*)"x", 1, 3, plmp->x);
    double** f = lammps_fix_external_get_force(plmp->lmp, "aphros_array");
    for (long i = 0; i < plmp->natoms; i++) {
      f[i][0] = 10 * plmp->x[3 * i];
      f[i][1] = 10 * plmp->x[3 * i + 1];
      f[i][2] = 10 * plmp->x[3 * i + 2];
    }
    lammps_command(plmp->lmp, "run 1");
  }
}

template <class M>
void InitHook(Hydro<M>* hydro) {
  static Lammps lmp;
  lmp.lmp = lammps_open_no_mpi(0, NULL, NULL);
  lammps_file(lmp.lmp, "in.lj");
  lmp.natoms = lammps_get_natoms(lmp.lmp);
  lmp.x = (double*)malloc(3 * lmp.natoms * sizeof *lmp.x);
  hydro->par_.ptr = (void*)&lmp;
}

using M = MeshCartesian<double, 3>;
template void StepHook(Hydro<M>*);
template void InitHook(Hydro<M>*);
