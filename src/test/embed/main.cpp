#undef NDEBUG
#include <mpi.h>
#include <cassert>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>

#include "distr/distrsolver.h"
#include "dump/output.h"
#include "dump/output_paraview.h"
#include "geom/mesh.h"
#include "geom/vect.h"
#include "kernel/kernelmeshpar.h"
#include "parse/vars.h"
#include "solver/embed.h"
#include "solver/reconst.h"
#include "solver/solver.h"
#include "util/suspender.h"

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
  using P::m;
  using P::var;

 private:
  using EB = Embed<M>;
  std::unique_ptr<EB> eb_;
  FieldCell<Scal> fcu_; // field
  FieldCell<Scal> fcum_; // field
  FieldCell<Scal> fct_; // field
};

template <class M>
void KernelEmbed<M>::Run() {
  auto sem = m.GetSem("Run");
  if (sem("init")) {
    FieldNode<Scal> fnf(m, 0);
    for (auto n : m.Nodes()) {
      auto x = m.GetNode(n);
      fnf[n] = 1.01 - Vect(x[0], x[1], x[2]).dot(Vect(1., 1., 0.));
    }
    eb_ = std::unique_ptr<EB>(new EB(m, fnf));
    fcu_.Reinit(m, 0);
    fct_.Reinit(m, 0);
    for (auto c : m.Cells()) {
      fct_[c] = size_t(eb_->GetType(c));
    }
  }
  if (sem.Nested("dumppoly")) {
    auto& eb = *eb_;
    eb.DumpPoly();
  }
  if (sem("dumpcsv")) {
    auto& eb = *eb_;
    using Type = typename EB::Type;
    std::ofstream out("eb.csv");
    out << "x,y,z\n";
    for (auto c : m.Cells()) {
      if (eb.GetType(c) == Type::cut) {
        auto x = eb.GetCellCenter(c);
        out << x[0] << "," << x[1] << "," << x[2] << "\n";
      }
    }
  }
  if (sem("dumpcsvface")) {
    auto& eb = *eb_;
    using Type = typename EB::Type;
    std::ofstream out("ebf.csv");
    out << "x,y,z\n";
    for (auto c : m.Cells()) {
      if (eb.GetType(c) == Type::cut) {
        auto x = eb.GetFaceCenter(c);
        out << x[0] << "," << x[1] << "," << x[2] << "\n";
      }
    }
    for (auto f : m.Faces()) {
      if (eb.GetType(f) == Type::cut) {
        auto x = eb.GetFaceCenter(f);
        out << x[0] << "," << x[1] << "," << x[2] << "\n";
      }
    }
  }
  for (size_t t = 0; t < 50; ++t) {
    auto& eb = *eb_;
    if (sem("step")) {
      using Type = typename EB::Type;
      const Scal a = 1.; // value on boundary
      const Vect vel(1., 0., 0.); // advection velocity
      const Scal dt = 0.1 * m.GetCellSize()[0] / vel.norm();
      fcum_ = fcu_;
      for (auto c : m.Cells()) {
        if (eb.GetType(c) == Type::cut) {
          Scal s = 0;
          for (auto q : m.Nci(c)) {
            const IdxFace f = m.GetFace(c, q);
            const Scal u = fcum_[m.GetCell(f, 0)];
            s += u * vel.dot(
                         m.GetNormal(f) * eb.GetArea(f) *
                         m.GetOutwardFactor(c, q));
          }
          s += a * (eb.GetNormal(c) * eb.GetArea(c)).dot(vel);
          fcu_[c] = fcum_[c] - s / eb.GetVolume(c) * dt;
        } else if (eb.GetType(c) == Type::regular) {
          Scal s = 0;
          for (auto q : m.Nci(c)) {
            const IdxFace f = m.GetFace(c, q);
            Scal u = fcum_[m.GetCell(f, 0)];
            s += u * vel.dot(m.GetOutwardSurface(c, q));
          }
          fcu_[c] = fcum_[c] - s / m.GetVolume(c) * dt;
        }
      }
      m.Dump(&fcu_, "u");
      m.Dump(&fct_, "type");

      FieldEmbed<Scal> feu(m, 0);
      /*
      for (auto f : m.Faces()) {
        feu[f] = eb.GetFaceCenter(f)[0];
      }
      for (auto c : m.Cells()) {
        feu[c] = eb.GetFaceCenter(c)[0];
      }
      */
      fcu_ = eb.Interpolate(feu);
    }
  }
  if (sem()) {
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

int main(int argc, const char** argv) {
  return RunMpi(argc, argv, Main);
}
