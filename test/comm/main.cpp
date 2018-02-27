#include <iostream>
#include <string>
#include <mpi.h>
#include <cassert>

#include "CubismDistr/Vars.h"
#include "CubismDistr/Interp.h"
#include "CubismDistr/Kernel.h"
#include "CubismDistr/Cubism.h"
#include "CubismDistr/Local.h"

#include "hydro/suspender.h"
#include "hydro/vect.hpp"
#include "hydro/mesh3d.hpp"
#include "hydro/solver.hpp"

template <class M>
M CreateMesh(const MyBlockInfo& bi) {
  using MIdx = typename M::MIdx;
  using B = MyBlock;
  using Vect = typename M::Vect;
  using Rect = geom::Rect<Vect>;
  B& b = *(B*)bi.ptrBlock;
  int hl = bi.hl;
  MIdx s(B::sx, B::sy, B::sz); // block size inner

  Scal h = bi.h_gridpoint;
  auto w = bi.index;   // block index
  auto c = bi.origin; 
  Vect d0(c[0], c[1], c[2]); // origin coord
  Vect d1 = d0 + Vect(s) * h;      // end coord
  Rect d(d0, d1);

  MIdx o(w[0] * s[0], w[1] * s[1], w[2] * s[2]); // origin index
  std::cout 
      << "o=" << o 
      << " dom=" << d0 << "," << d1 
      << " h=" << h
      << std::endl;
  
  return geom::InitUniformMesh<M>(d, o, s, hl);
}

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

template <class M>
void Simple<M>::Run() {
  auto sem = m.GetSem("run");
  auto f = [](Vect x) { return x[0]; };
  if (sem("init")) {
    auto a = m.GetSem("run2");
    if (a()) {
      fc_.Reinit(m);
      for (auto i : m.Cells()) {
        fc_[i] = f(m.GetCenter(i));
      }
    }
  }
  if (sem("check")) {
    fc_.Reinit(m);
    for (auto i : m.Cells()) {
      assert(fc_[i] == f(m.GetCenter(i)));
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
