#undef NDEBUG
#include <cassert>
#include <fstream>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <utility>
#include <tuple>

#include "distr/distrbasic.h"
#include "solver/embed.h"
#include "solver/reconst.h"

using M = MeshStructured<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using R = Reconst<Scal>;
using EB = Embed<M>;
using Type = typename EB::Type;

std::tuple<Vect, Vect, Vect> GetStat(
    const FieldCell<Vect>& fcg, const EB& eb) {
  Vect mean(0);
  Vect max(-1e10);
  Vect min(1e10);
  Scal sumv = 0;
  for (auto c : eb.Cells()) {
    auto& g = fcg[c];
    const Scal v = eb.GetVolume(c);
    mean += g * v;
    sumv += v;
    max = max.max(g);
    min = min.min(g);
  }
  mean /= sumv;
  return {mean, min, max};
}

void PrintStat(const std::tuple<Vect, Vect, Vect>& s) {
  std::cout << "mean=" << std::get<0>(s) << " ";
  std::cout << "min=" << std::get<1>(s) << " ";
  std::cout << "max=" << std::get<2>(s) << " ";
  std::cout << std::endl;
}

FieldEmbed<Scal> GetNoise(const EB& eb, size_t seed) {
  FieldEmbed<Scal> fe(eb.GetMesh());
  std::default_random_engine gen(seed);
  std::uniform_real_distribution<double> dis(-1, 1);
  for (auto f : eb.GetMesh().AllFaces()) {
    fe[f] = dis(gen);
  }
  return fe;
}

FieldEmbed<Scal> Add(
    const FieldEmbed<Scal>& fea, const FieldEmbed<Scal>& feb, const EB& eb) {
  auto fe = fea;
  for (auto f : eb.Faces()) {
    fe[f] += feb[f];
  }
  for (auto c : eb.CFaces()) {
    fe[c] += feb[c];
  }
  return fe;
}
void Run(M& m, Vars& var) {
  auto sem = m.GetSem("Run");
  struct {
    std::unique_ptr<EB> eb;
    FieldCell<Vect> fcg;
  } * ctx(sem);
  auto& fcg = ctx->fcg;

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
    auto func = [](Vect x) { return x[0]; };
    for (auto f : eb.Faces()) {
      feu[f] = func(eb.GetFaceCenter(f));
    }
    for (auto c : eb.CFaces()) {
      feu[c] = func(eb.GetFaceCenter(c));
    }
    feu = Add(feu, GetNoise(eb, 0), eb);
    fcg = eb.Gradient(feu); // compact gradient
    std::cout << "1/h=" << 1 / m.GetCellSize()[0] << std::endl;
    PrintStat(GetStat(fcg, eb));
    fcg = eb.AverageCutCells(fcg);
    PrintStat(GetStat(fcg, eb));

    m.Dump(&fcg, 0, "gx");
    m.Dump(&fcg, 1, "gy");
    m.Dump(&fcg, 2, "gz");
  }
  if (sem()) {}
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
set vect eb_box_r 0.251 0.251 0.251
#set vect eb_box_r 10 0.249 10

set string eb_init sphere
set vect eb_sphere_c 0.5 0.5 0.5
set vect eb_sphere_r 0.249 0.249 0.249
set double eb_sphere_angle 0

set int eb_init_inverse 1
)EOF";
    return RunMpiBasic<M>(argc, argv, Run, conf);
  }
