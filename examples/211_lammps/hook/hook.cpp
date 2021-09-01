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
  auto* lmp = static_cast<Lammps*>(hydro->par_.ptr);

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
    lammps_gather_atoms(lmp->lmp, (char*)"x", 1, 3, lmp->x);
    double** f = lammps_fix_external_get_force(lmp->lmp, "aphros_array");
    for (long i = 0; i < lmp->natoms; i++) {
      f[i][0] = 10 * lmp->x[3 * i];
      f[i][1] = 10 * lmp->x[3 * i + 1];
      f[i][2] = 10 * lmp->x[3 * i + 2];
    }
    lammps_command(lmp->lmp, "run 1");
  }
}

template <class M>
void InitHook(Hydro<M>* hydro) {
    auto& m = hydro->m;
  auto sem = m.GetSem();
  if (sem()) {
    if (m.IsLead()) { // Executed only on the lead block, once in each rank.
		      // There is only one lead block in each rank.
      // Allocate memory for the shared state and save pointer in the solver.
      int argc;
      int dash;
      const char **argv;
      const char *path;
      const char **a;

      hydro->par_.ptr = new Lammps();
      auto* lmp = static_cast<Lammps*>(hydro->par_.ptr);
      dash = sysinfo::misc.arg_after_double_dash;
      if (dash == 0) {
        const char *argv_none[] = {"aphros", NULL};
        argc = 0;
        argv = argv_none;
      } else {
        argc = sysinfo::misc.argc - dash - 1;
        argv = sysinfo::misc.argv + dash - 1;
      }
      path = "in.lj";
      for (a = argv; *++a != NULL;)
        if (a[0][0] == '-') switch (a[0][1]) {
          case 'i':
            a++;
            if ((path = *a) == NULL) {
              fprintf(stderr, "hook: -i needs an argumetn\n");
              fassert(0);
            }
            break;
          case 'h':
            setenv("PAGER", "cat", 1);
            const char *argv_help[] = {"aphros", "-h", NULL};
            lammps_open_no_mpi(2, const_cast<char**>(argv_help), NULL);
            fassert(0);
            break;
          }
      lmp->lmp = lammps_open_no_mpi(argc, const_cast<char**>(argv), NULL);
      lammps_file(lmp->lmp, path);
      lmp->natoms = lammps_get_natoms(lmp->lmp);
      lmp->x = (double*)malloc(3 * lmp->natoms * sizeof(double));
      if (lmp->x == NULL)
	fassert(false, "fail to allocate memory for data");
    }
    if (m.IsRoot()) { // Executed only on the root block, once over all ranks.
		      // There is only one root block over all ranks.
      //auto* lmp = static_cast<Lammps*>(hydro->par_.ptr);
    }
  }
}

template <class M>
void FinalHook(Hydro<M>* hydro) {
  auto& m = hydro->m;
  auto sem = m.GetSem();
  auto* lmp = static_cast<Lammps*>(hydro->par_.ptr);
  if (sem()) {
    if (m.IsRoot()) {
      std::cout << util::Format(
	  "FinalHook t={:} \n\n", hydro->fs_->GetTime());
    }
    if (m.IsLead()) {
      free(lmp->x);
      delete lmp;
    }
  }
}

using M = MeshCartesian<double, 3>;
template void StepHook(Hydro<M>*);
template void InitHook(Hydro<M>*);
template void FinalHook(Hydro<M>*);
