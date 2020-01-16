// Created by Sergey Litvinov on 12.08.2019
// Copyright 2019 ETH Zurich

#undef NDEBUG
#include <cassert>
#include <cmath>
#include <functional>
#include <iostream>
#include <sstream>

#include "distr/distrsolver.h"
#include "geom/mesh.h"
#include "kernel/kernelmeshpar.h"
#include "parse/vars.h"
#include "solver/solver.h"
using namespace solver;

const int dim = 3;
using MIdx = GMIdx<dim>;
using Scal = double;
using Vect = GVect<Scal, dim>;
using Mesh = MeshStructured<Scal, dim>;

typedef void (*TFunc)(int, int, int, double*);
TFunc kFunc;

extern "C" {
int CMain(int, const char**, TFunc);
}

struct GPar {};

template <class M_>
class Simple : public KernelMeshPar<M_, GPar> {
 public:
  using P = KernelMeshPar<M_, GPar>; // parent
  using M = M_;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using MIdx = typename M::MIdx;
  using Par = GPar;
  static constexpr size_t dim = M::dim;

  using P::P;
  void Run() override;

 protected:
  using P::bi_;
  using P::m;
  using P::var;

 private:
  FieldCell<Scal> fcu_;
};

template <class M>
void Simple<M>::Run() {
  auto sem = m.GetSem();
  if (sem()) {
    fcu_.Reinit(m);
    for (auto c : m.AllCells()) {
      fcu_[c] = m.GetCenter(c).norm();
    }
    MIdx ws = m.GetIndexCells().GetSize();
    kFunc(ws[0], ws[1], ws[2], fcu_.data());
    m.Dump(&fcu_, "u");
  }
  if (sem()) {
  }
}

void Main(MPI_Comm comm, Vars& var) {
  using M = MeshStructured<double, 3>;
  using K = Simple<M>;
  using Par = typename K::Par;
  Par par;

  DistrSolver<M, K> ds(comm, var, par);
  ds.Run();
}

int CMain(int argc, const char** argv, TFunc f) {
  kFunc = f;
  return RunMpi(argc, argv, Main);
}
