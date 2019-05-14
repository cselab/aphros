#undef NDEBUG
#include <iostream>
#include <string>
#include <mpi.h>
#include <cassert>
#include <iomanip>
#include <fstream>
#include <functional>
#include <utility>
#include <tuple>
#include <sstream>
#include <memory>

#include "parse/vars.h"
#include "kernel/kernelmeshpar.h"
#include "distr/distrsolver.h"
#include "util/suspender.h"
#include "geom/vect.h"
#include "geom/mesh.h"
#include "solver/solver.h"
#include "solver/embed.h"
#include "solver/reconst.h"
#include "dump/output.h"
#include "dump/output_paraview.h"

struct GPar {};

template <class M_>
class KernelEmbed : public KernelMeshPar<M_, GPar> {
 public:
  using P = KernelMeshPar<M_, GPar>; // parent
  using M = M_;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using MIdx = typename M::MIdx;
  using Par = GPar;
  using R = Reconst<Scal>;
  static constexpr size_t dim = M::dim;

  using P::P;
  void Run() override;

 protected:
  using P::var;
  using P::m;

 private:
  using EM = solver::Embed<M>;
  std::unique_ptr<EM> em_;
  FieldCell<Scal> fct_;  // cell type
  FieldCell<Scal> fcu_;  // field
};


template <class M>
void KernelEmbed<M>::Run() {
  auto sem = m.GetSem("Run");
  if (sem("init")) {
    FieldNode<Scal> fnf(m, 0);
    for (auto n : m.Nodes()) {
      auto x = m.GetNode(n);
      fnf[n] = x.norm() *
        (1 + std::sin(x[0] * 5)*0.1 + std::sin(x[1] * 7)*0.15) - 1.;
      //fnf[n] *= -1;
    }
    em_ = std::unique_ptr<EM>(new EM(m, fnf));
    fct_.Reinit(m);
    for (auto c : m.Cells()) {
      fct_[c] = size_t(em_->GetCellType()[c]);
    }
    m.Dump(&fct_, "type");

    fcu_.Reinit(m, 0.);
    using Type = typename EM::Type;
    auto& ffs = em_->GetFaceArea();
    auto& fcs = em_->GetCellArea();
    auto& fct = em_->GetCellType();
    FieldFace<Scal> ffg(m, 0); // gradient on faces
    FieldCell<Scal> fcg(m, 0); // gradient on embedded boundary
    Scal a = 1.; // value on boundary
    // gradient on faces
    for (auto f : m.Faces()) {
      auto cm = m.GetNeighbourCell(f, 0);
      auto cp = m.GetNeighbourCell(f, 1);
      if (fct[cm] != Type::excluded && fct[cp] != Type::excluded) {
        ffg[f] = (fcu_[cp] - fcu_[cm]) * m.GetArea(f) / m.GetVolume(cp);
      } else if (fct[cm] != Type::excluded) {
        ffg[f] = (a - fcu_[cm]) * m.GetArea(f) / m.GetVolume(cm);
      } else if (fct[cp] != Type::excluded) {
        ffg[f] = (fcu_[cp] - a) * m.GetArea(f) / m.GetVolume(cp);
      }
    }
    // gradient on embedded boundary
    for (auto c : m.Cells()) {
      auto xn = R::GetNearest(
          Vect(0), em_->GetNormal()[c], em_->GetPlane()[c], m.GetCellSize());
      fcg[c] = (a - fcu_[c]) / xn.norm();
    }
    FieldCell<Scal> fcq(m, 0);
    for (auto c : m.Cells()) {
      if (fct[c] == Type::cut) {
        Scal vol = 0;
        for (auto q : m.Nci(c)) {
          IdxFace f = m.GetNeighbourFace(c, q);
          auto x = m.GetCenter(c);
          fcq[c] += ffg[f] * ffs[f] * m.GetOutwardFactor(c, q);
          vol += x[0] * m.GetOutwardSurface(c, q)[0];
        }
        auto n = em_->GetNormal()[c];
        n /= n.norm();
        auto xn = R::GetNearest(
            Vect(0), em_->GetNormal()[c], 
            em_->GetPlane()[c], m.GetCellSize());
        fcq[c] += fcg[c] * fcs[c];
        vol += xn[0] * n[0] * fcs[c];
        fcq[c] /= vol;
      } else if (fct[c] == Type::regular) {
        for (auto q : m.Nci(c)) {
          IdxFace f = m.GetNeighbourFace(c, q);
          fcq[c] += ffg[f] * m.GetArea(f) * m.GetOutwardFactor(c, q);
        }
        fcq[c] /= m.GetVolume(c);
      }
    }
    Scal dt = 0.01;
    for (auto c : m.Cells()) {
      if (fct[c] != Type::excluded) {
        fcu_[c] += dt * fcq[c];
      }
    }
    m.Dump(&fcu_, "u");

  }
  if (sem.Nested("dump")) {
    em_->DumpPoly();
  }
}

using Scal = double;
using M = MeshStructured<Scal, 3>;
using K = KernelEmbed<M>;

void Main(MPI_Comm comm, Vars& var0) {
  int rank;
  MPI_Comm_rank(comm, &rank);

  Vars var = var0;

  using Par = typename K::Par;
  Par par;

  DistrSolver<M, K> ds(comm, var, par);
  ds.Run();
}


int main (int argc, const char ** argv) {
  return RunMpi(argc, argv, Main);
}
