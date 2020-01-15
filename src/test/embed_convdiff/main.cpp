#undef NDEBUG
#include <cassert>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <utility>

#include "distr/distrbasic.h"
#include "solver/convdiffe_eb.h"
#include "solver/embed.h"
#include "solver/reconst.h"

using M = MeshStructured<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using EB = Embed<M>;
using Type = typename EB::Type;
using CD = ConvDiffScalExpEmbed<M>;

FieldNode<Scal> GetLevelSet(const M& m) {
  FieldNode<Scal> fnl(m, 0);
  for (auto n : m.AllNodes()) {
    const Vect x = m.GetNode(n);
    const Vect xc(0.5, 0.5, 0.5);
    const Vect s(1., 1., 1.);
    auto rot = [](Vect xx) {
      const Scal a = 3.14159265358979323 * 0;
      const Scal as = std::sin(a);
      const Scal ac = std::cos(a);
      const Scal x = xx[0];
      const Scal y = xx[1];
      const Scal z = xx[2];
      return Vect(x * ac - y * as, x * as + y * ac, z);
    };
    fnl[n] = (rot(x - xc) / s).norm() - 0.2;
  }
  return fnl;
}

void Run(M& m, Vars&) {
  auto sem = m.GetSem("Run");
  struct {
    std::unique_ptr<EB> eb;
    std::unique_ptr<CD> cd;
    FieldCell<Scal> fcu;
    MapCondFace mfc;
    FieldCell<Scal> fcr;
    FieldEmbed<Scal> fed;
    FieldCell<Scal> fcs;
    FieldFace<Scal> ffv;
    size_t frame = 0;
  } * ctx(sem);
  auto& eb = ctx->eb;
  auto& cd = ctx->cd;
  auto& fcu = ctx->fcu;
  auto& ffv = ctx->ffv;
  auto& mfc = ctx->mfc;
  auto& frame = ctx->frame;

  if (sem("init")) {
    eb.reset(new EB(m, GetLevelSet(m)));
    ctx->fcr.Reinit(m, 1);
    ctx->fed.Reinit(m, 1);
    ctx->fcs.Reinit(m, 0);
    ffv.Reinit(m, 0);
    //const Vect vel(1., 0.5, 0.25);
    const Vect vel(0);
    for (auto f : m.Faces()) {
      ffv[f] = vel.dot(m.GetSurface(f));
    }
    fcu.Reinit(m);
    for (auto c : eb->AllCells()) {
      const Scal a = 12;
      fcu[c] =
          std::sin(m.GetCenter(c)[0] * a) * std::sin(m.GetCenter(c)[1] * a);
    }
    const size_t bc = 0;
    const Scal bcu = 0;
    const Scal dt = 0.001;
    typename CD::Par par;
    cd.reset(new CD(
        m, *eb, fcu, mfc, bc, bcu, &ctx->fcr, &ctx->fed, &ctx->fcs, &ctx->ffv,
        0, dt, par));
  }
  if (sem.Nested("dumppoly")) {
    eb->DumpPoly();
  }
  const size_t maxt = 100;
  const size_t nfr = 100;
  for (size_t t = 0; t < maxt; ++t) {
    if (sem.Nested("convdiff-start")) {
      cd->StartStep();
    }
    if (sem.Nested("convdiff-iter")) {
      cd->MakeIteration();
    }
    if (sem.Nested("convdiff-finish")) {
      cd->FinishStep();
    }
    if (t % std::max<size_t>(1, maxt / nfr) != 0) {
      continue;
    }
    if (sem("dump")) {
      m.Dump(&cd->GetField(), "u");
      ++frame;
    }
  }
  if (sem()) {
  }
}

int main(int argc, const char** argv) {
  std::string conf = R"EOF(
set int bx 1
set int by 1
set int bz 1

set int bsx 16
set int bsy 16
set int bsz 16
)EOF";
  return RunMpiBasic<M>(argc, argv, Run, conf);
}
