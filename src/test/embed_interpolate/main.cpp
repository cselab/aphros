// Created by Petr Karnakov on 05.07.2020
// Copyright 2020 ETH Zurich

#undef NDEBUG
#include <cassert>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>

#include "distr/distrsolver.h"
#include "dump/hdf.h"
#include "kernel/kernelmeshpar.h"
#include "parse/vars.h"
#include "solver/approx.h"
#include "solver/approx_eb.h"
#include "solver/embed.h"

template <class M_>
struct GPar {};

template <class M_>
class EmbedInterpolate : public KernelMeshPar<M_, GPar<M_>> {
 public:
  using P = KernelMeshPar<M_, GPar<M_>>;
  using M = M_;
  using Mesh = M;
  using Par = typename P::Par;
  static constexpr size_t dim = M::dim;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using Sem = typename Mesh::Sem;
  using BlockInfoProxy = generic::BlockInfoProxy<dim>;

  EmbedInterpolate(Vars& var, const BlockInfoProxy& block, Par& par);
  void Run() override;

  using P::m;
  using P::var;
  std::unique_ptr<Embed<M>> eb_;
};

template <class M>
EmbedInterpolate<M>::EmbedInterpolate(
    Vars& var_, const BlockInfoProxy& block, Par& par)
    : KernelMeshPar<M, Par>(var_, block, par) {}

template <class M>
void EmbedInterpolate<M>::Run() {
  auto sem = m.GetSem();
  struct {
    FieldNode<Scal> fnl;
    FieldCell<Scal> fcvx;
    FieldCell<Scal> fcvy;
    FieldCell<Scal> fcvz;
    FieldCell<Scal> fcomz;

    std::vector<Vect> px;
    std::vector<Scal> pomz;
    std::vector<Scal> pomz_interp;
  } * ctx(sem);
  auto& px = ctx->px;
  auto& pomz = ctx->pomz;
  auto& pomz_interp = ctx->pomz_interp;
  if (sem("ctor")) {
    eb_.reset(new Embed<M>(m, var.Double["embed_gradlim"]));
  }
  if (sem.Nested("levelset")) {
    UEmbed<M>::InitLevelSet(ctx->fnl, m, var, m.IsRoot());
  }
  if (sem.Nested("smoothen")) {
    SmoothenNode(ctx->fnl, m, var.Int["embed_smoothen_iters"]);
  }
  if (sem.Nested("init")) {
    eb_->Init(ctx->fnl);
  }
  if (sem.Nested()) {
    eb_->DumpPoly(var.Int["vtkbin"], var.Int["vtkmerge"]);
  }
  if (sem.Nested()) {
    dump::Hdf<M>::Read(ctx->fcvx, var.String["input_hdf_vx"], m);
  }
  if (sem.Nested()) {
    dump::Hdf<M>::Read(ctx->fcvy, var.String["input_hdf_vy"], m);
  }
  if (sem.Nested()) {
    dump::Hdf<M>::Read(ctx->fcomz, var.String["input_hdf_omz"], m);
  }
  if (sem("reduce_x")) {
    auto& eb = *eb_;
    for (auto c : eb.CFaces()) {
      px.push_back(eb.GetFaceCenter(c));
    }
    m.Reduce(&px, Reduction::concat);
  }
  if (sem("calc_omz")) {
    auto& eb = *eb_;
    MapEmbed<BCond<Scal>> mebc;
    for (auto c : eb.CFaces()) {
      mebc[c] = BCond<Scal>(BCondType::dirichlet, 0);
    }
    const auto fegx = UEmbed<M>::Gradient(ctx->fcvx, mebc, eb);
    const auto fegy = UEmbed<M>::Gradient(ctx->fcvy, mebc, eb);
    FieldEmbed<Scal> feomz(m, 0);
    for (auto c : eb.CFaces()) {
      const Vect n = -eb.GetNormal(c);
      const Vect t{n[1], -n[0], 0.};
      feomz[c] = fegx[c] * t[0] + fegy[c] * t[1];
    }
    for (auto c : eb.CFaces()) {
      pomz.push_back(feomz[c]);
    }
    m.Reduce(&pomz, Reduction::concat);
  }
  if (sem("interpolate_omz")) {
    auto& eb = *eb_;
    MapEmbed<BCond<Scal>> mebc;
    for (auto c : eb.CFaces()) {
      mebc[c] = BCond<Scal>(BCondType::neumann, 0);
    }
    const auto feomz = UEmbed<M>::Interpolate(ctx->fcomz, mebc, eb);
    for (auto c : eb.CFaces()) {
      pomz_interp.push_back(feomz[c]);
    }
    m.Reduce(&pomz_interp, Reduction::concat);
  }
  if (sem("write") && m.IsRoot()) {
    fassert_equal(px.size(), pomz.size());
    fassert_equal(px.size(), pomz_interp.size());
    std::ofstream o(var.String["output_csv_omz"]);
    o.precision(16);
    // header
    o << "x,y,z,omz,omz_interp";
    o << std::endl;
    // content
    for (size_t i = 0; i < px.size(); ++i) {
      o << px[i][0] << ',' << px[i][1] << ',' << px[i][2];
      o << ',' << pomz[i];
      o << ',' << pomz_interp[i];
      o << "\n";
    }
  }
  if (sem()) { // XXX empty stage
  }
}

void Main(MPI_Comm comm, Vars& var) {
  using M = MeshCartesian<double, 3>;
  using K = EmbedInterpolate<M>;
  using Par = typename K::Par;
  Par par;

  DistrSolver<M, K> ds(comm, var, par);
  ds.Run();
}

int main(int argc, const char** argv) {
  return RunMpi(argc, argv, Main);
}
