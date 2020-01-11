#undef NDEBUG
#include <cassert>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <utility>

#include "distr/distrbasic.h"
#include "solver/embed.h"
#include "solver/reconst.h"

using M = MeshStructured<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using R = Reconst<Scal>;
using EB = Embed<M>;
using Type = typename EB::Type;

void PrintStat(const FieldCell<Vect>& fcg, const M& m, const EB& eb) {
  Vect mean(1e10);
  Vect max(-1e10);
  Vect min(0);
  Scal w = 0;
  for (auto c : eb.Cells()) {
    auto& g = fcg[c];
    mean += g;
    max = max.max(g);
    min = min.min(g);
    w += 1;
  }
  mean /= w;
  std::cout << "mean=" << mean << std::endl;
  std::cout << "min=" << min << std::endl;
  std::cout << "max=" << max << std::endl;
}

void Run(M& m, Vars& var) {
  auto sem = m.GetSem("Run");
  struct {
    std::unique_ptr<EB> eb;
  } * ctx(sem);

  if (sem("init")) {
    auto fnl = InitEmbed(m, var);
    ctx->eb.reset(new EB(m, fnl));
  }
  if (sem.Nested("dumppoly")) {
    auto& eb = *ctx->eb;
    eb.DumpPoly();
  }
  auto& eb = *ctx->eb;
  if (sem("step")) {
    FieldEmbed<Scal> feu(m, 0);
    const FieldCell<Vect> fcg = eb.Gradient(feu); // compact gradient
    PrintStat(fcg, m, eb);
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

set string eb_init box
set vect eb_box_c 0.5 0.5 0.5
#set vect eb_box_r 0.249 0.249 0.249
set vect eb_box_r 10 0.249 10

set string eb_init sphere
set vect eb_sphere_c 0.5 0.5 0.5
set vect eb_sphere_r 0.249 0.249 0.249
set double eb_sphere_angle 0
)EOF";
  return RunMpiBasic<M>(argc, argv, Run, conf);
}
