// Created by Petr Karnakov on 29.08.2021
// Copyright 2021 ETH Zurich

#include <iostream>

#include <kernel/hydro.h>
#include <library.h>
#include <util/format.h>
#include <util/posthook.h>
#include <func/init_vel.h>

struct Lammps {
  void* lmp;
  long natoms;
  double* x;
  double* v;
};

template <class M>
class CustomVelocity : public ModuleInitVelocity<M> {
 public:
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  CustomVelocity() : ModuleInitVelocity<M>("custom") {}
  void operator()(FieldCell<Vect>& fcv, const Vars&, const M& m) override {
    for (auto c : m.AllCellsM()) {
      const Vect x = c.center() * (2 * M_PI);
      Vect& v = fcv[c];
      v[0] = std::sin(x[0]) * std::cos(x[1]) * std::cos(x[2]);
      v[1] = -std::cos(x[0]) * std::sin(x[1]) * std::cos(x[2]);
      v[2] = 0;
    }
  }
};

template <class M>
void StepHook(Hydro<M>* hydro) {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;

  auto& m = hydro->m;
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
    double x, y, z, u, v, w;
    double** f;
    t.velocity_mean /= t.volume;
    lammps_gather_atoms(lmp->lmp, (char*)"x", 1, 3, lmp->x);
    lammps_gather_atoms(lmp->lmp, (char*)"v", 1, 3, lmp->v);
    f = lammps_fix_external_get_force(lmp->lmp, "aphros_array");
    const FieldCell<Vect>& fcv = hydro->fs_->GetVelocity();
    for (long i = 0; i < lmp->natoms; i++) {
      x = lmp->x[3 * i];
      y = lmp->x[3 * i + 1];
      z = lmp->x[3 * i + 2];
      u = lmp->v[3 * i];
      v = lmp->v[3 * i + 1];
      w = lmp->v[3 * i + 2];      
      const Vect& vel = fcv[m.GetCellFromPoint(Vect(x, y, z))];
      f[i][0] = (vel[0] - u)/64;
      f[i][1] = (vel[1] - v)/64;
      f[i][2] = (vel[2] - w)/64;
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
      path = "in.polymers";
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
      lmp->v = (double*)malloc(3 * lmp->natoms * sizeof(double));
      if (lmp->x == NULL || lmp->v == NULL)
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
      free(lmp->v);
      delete lmp;
    }
  }
}

using M = MeshCartesian<double, 3>;
template void StepHook(Hydro<M>*);
template void InitHook(Hydro<M>*);
template void FinalHook(Hydro<M>*);
bool reg[] = {
    RegisterModule<CustomVelocity<M>>(),
};
