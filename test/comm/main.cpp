#include <iostream>
#include <string>
#include <mpi.h>
#include <cassert>

#include "CubismDistr/Vars.h"
#include "CubismDistr/Interp.h"
#include "CubismDistr/Kernel.h"
#include "CubismDistr/Cubism.h"
#include "CubismDistr/Local.h"
#include "CubismDistr/Hydro.h"

#include "hydro/suspender.h"
#include "hydro/vect.hpp"
#include "hydro/mesh3d.hpp"
#include "hydro/solver.hpp"


template <class M>
class Simple : public Kernel {
 public:
  using Mesh = M;
  using Scal = double;
  using Vect = typename Mesh::Vect;
  using MIdx = typename Mesh::MIdx;
  static constexpr size_t dim = M::dim;

  Simple(Vars& par, const MyBlockInfo& bi);
  void Run() override;
  M& GetMesh() { return m; }

 private:
  Vars& par;
  std::string name_;
  MyBlockInfo bi_;
  M m;
  geom::FieldCell<Scal> fc_;
};

template <class _M>
class SimpleFactory : public KernelFactory {
 public:
  using M = _M;
  using K = Simple<M>;
  std::unique_ptr<Kernel> Make(Vars& par, const MyBlockInfo& bi) override {
    return std::unique_ptr<K>(new K(par, bi));
  }
};

template <class M>
Simple<M>::Simple(Vars& par, const MyBlockInfo& bi) 
  : par(par), bi_(bi), m(CreateMesh<M>(bi))
{
  name_ = 
      "[" + std::to_string(bi.index[0]) +
      "," + std::to_string(bi.index[1]) +
      "," + std::to_string(bi.index[2]) + "]";
}

bool Cmp(Scal a, Scal b) {
  return std::abs(a - b) < 1e-12;
}

template <class M>
void Simple<M>::Run() {
  auto sem = m.GetSem("run");
  auto f = [](Vect v) { 
    for (int i = 0; i < dim; ++i) {
      while (v[i] < 0.) {
        v[i] += 1.;
      }
      while (v[i] > 1.) {
        v[i] -= 1.;
      }
    }
    return std::sin(v[0]) * std::cos(v[1]) * std::exp(v[2]); 
  };
  auto& bc = m.GetBlockCells();
  if (sem("init")) {
    fc_.Reinit(m);
    for (auto i : m.Cells()) {
      fc_[i] = f(m.GetCenter(i));
    }
    m.Comm(&fc_);
  }
  if (sem("check")) {
    for (auto i : m.AllCells()) {
      auto x = m.GetCenter(i);
      if (!Cmp(fc_[i], f(m.GetCenter(i)))) {
        std::cerr 
          << bc.GetMIdx(i) << " " 
          << fc_[i] << " != " << f(x) << " "
          << std::endl;
        assert(false);
      }
    }
  }
}


void Main(MPI_Comm comm, bool loc, Vars& par) {
  // read config files, parse arguments, maybe init global fields
  using M = geom::MeshStructured<Scal, 3>;
  using K = Simple<M>;
  using KF = SimpleFactory<M>;
  using D = Distr;
  
  KF kf;

  const int es = 8;
  const int hl = par.Int["hl"];
  const int bs = 16;
  
  // Initialize buffer mesh and make Simple for each block.
  std::unique_ptr<Distr> d;
  if (loc) {
    d.reset(new Local<KF>(comm, kf, bs, es, hl, par));
  } else {
    d.reset(new Cubism<KF>(comm, kf, bs, es, hl, par));
  }

  while (!d->IsDone()) {
    d->Step();
  }
}


int main (int argc, const char ** argv) {
  int prov;
  MPI_Init_thread(&argc, (char ***)&argv, MPI_THREAD_MULTIPLE, &prov);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  MPI_Comm comm;

  Vars par;
  Interp ip(par);

  std::ifstream f("a.conf");
  ip.RunAll(f);
  if (rank == 0) {
    ip.PrintAll();
  }

  bool loc = par.Int["loc"];

  if (loc) {
    MPI_Comm_split(MPI_COMM_WORLD, rank, rank, &comm);
    if (rank == 0) {
      Main(comm, loc, par);
    }
  } else {
    comm = MPI_COMM_WORLD;
    Main(comm, loc, par);
  }


  MPI_Finalize();	
  return 0;
}
