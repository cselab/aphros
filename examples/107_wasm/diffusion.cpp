// Created by Petr Karnakov on 04.03.2021
// Copyright 2021 ETH Zurich

#include <stdint.h>
#include <stdlib.h>
#include <emscripten.h>
#include <emscripten/html5.h>
#include <cmath>
#include <cstring>
#include <iostream>
#include <memory>
#include <sstream>

#include "distr/distrbasic.h"
#include "distr/distrsolver.h"
#include "func/init_u.h"
#include "geom/mesh.h"
#include "geom/rangemulti.h"
#include "solver/approx.h"
#include "solver/approx_eb.h"
#include "solver/cond.h"
#include "solver/solver.h"
#include "util/format.h"
#include "util/module.h"
#include "util/timer.h"

using M = MeshCartesian<double, 2>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using MIdx = typename M::MIdx;

void CopyToCanvas(uint32_t* ptr, int w, int h) {
  EM_ASM_(
      {
        let data = Module.HEAPU8.slice($0, $0 + $1 * $2 * 4);
        let context = Module['canvas'].getContext('2d');
        let imageData = context.getImageData(0, 0, $1, $2);
        imageData.data.set(data);
        context.putImageData(imageData, 0, 0);
      },
      ptr, w, h);
}

struct Par {};

class Solver : public KernelMeshPar<M, Par> {
 public:
  using P = KernelMeshPar<M, Par>; // parent
  using Par = typename P::Par;
  static constexpr size_t dim = M::dim;
  using Sem = typename M::Sem;
  using BlockInfoProxy = generic::BlockInfoProxy<dim>;

  Solver(Vars& var_, const BlockInfoProxy& proxy, Par& par)
      : KernelMeshPar<M, Par>(var_, proxy, par), fcu(m, 0), ff_flux(m, 0) {}
  void Run() override;

 protected:
  using P::m;
  using P::par_;
  using P::var;

 private:
  FieldCell<Scal> fcu;
  FieldEmbed<Scal> ff_flux;
};

struct Canvas {
  Canvas(MIdx size_) : size(size_), buf(size.prod()) {}
  MIdx size;
  std::vector<uint32_t> buf;
};

struct State {
  State(MPI_Comm comm, Vars& var) : distrsolver(comm, var, par) {}
  typename Solver::Par par;
  DistrSolver<M, Solver> distrsolver;
  int step = 0;
  Scal time = 0;
  Scal dt = 0;
  Scal diffusion = 0.01;
  Vect velocity{0, 1};
  bool pause = false;
  bool to_init_field = true;
  bool to_init_solver = true;
};

std::shared_ptr<State> g_state;
std::shared_ptr<Canvas> g_canvas;

static Scal Clamp(Scal f) {
  return f < 0 ? 0 : f > 1 ? 1 : f;
}

static void Init(FieldCell<Scal>& fcu, const M& m) {
  Vars par;
  par.String.Set("init_vf", "circlels");
  par.Vect.Set("circle_c", {0.5, 0.5});
  par.Double.Set("circle_r", 0.2);
  par.Int.Set("dim", 2);
  auto func = CreateInitU<M>(par, false);
  func(fcu, m);
}

static void Step(
    FieldCell<Scal>& fcu, Scal dt, Scal diffusion,
    const FieldFace<Scal>& ff_flux, const M& m) {
  using UEB = UEmbed<M>;
  const auto ffg = UEB::Gradient(fcu, {}, m);
  const auto fcg = UEB::AverageGradient(ffg, m);
  auto ffu = UEB::InterpolateUpwind(fcu, {}, ConvSc::superbee, fcg, ff_flux, m);

  for (auto c : m.CellsM()) {
    Scal sum = 0;
    for (auto q : m.Nci(c)) {
      auto f = c.face(q);
      sum += diffusion * ffg[f] * c.outward_factor(q) * f.area;
      sum -= ff_flux[f] * ffu[f] * c.outward_factor(q);
    }
    fcu[c] += dt * sum / c.volume;
  }
}

static void Render(Canvas& canvas, const FieldCell<Scal>& fcu, const M& m) {
  const auto msize = m.GetGlobalSize();
  for (auto c : m.CellsM()) {
    const MIdx w(c);
    const MIdx start = w * canvas.size / msize;
    const MIdx end = (w + MIdx(1)) * canvas.size / msize;
    unsigned char qr = 255 * Clamp(fcu[c]);
    unsigned char qg = 0;
    unsigned char qb = 0;
    uint32_t q = 0xff000000 | (qr << 0) | (qg << 8) | (qb << 16);
    for (int y = start[1]; y < end[1]; y++) {
      for (int x = start[0]; x < end[0]; x++) {
        canvas.buf[canvas.size[0] * (canvas.size[1] - y - 1) + x] = q;
      }
    }
  }
}

void Solver::Run() {
  auto sem = m.GetSem();
  auto state = g_state;
  auto& s = *state;
  if (sem("init") && s.to_init_field) {
    Init(fcu, m);
  }
  if (sem("flux")) {
    for (auto f : m.Faces()) {
      ff_flux[f] = s.velocity.dot(m.GetSurface(f));
    }
  }
  if (sem()) {
    s.to_init_field = false;
    const Scal h = m.GetCellSize()[0];
    s.dt = std::min(
        0.125 * h * h / s.diffusion, 0.25 * h / s.velocity.abs().max());
    Step(fcu, s.dt, s.diffusion, ff_flux, m);
    Render(*g_canvas, fcu, m);
    m.Comm(&fcu);
    if (m.IsRoot()) {
      s.time += s.dt;
      ++s.step;
    }
  }
}

static void main_loop() {
  auto state = g_state;
  auto& s = *state;
  if (s.pause) {
    return;
  }
  std::memset(g_canvas->buf.data(), 0, g_canvas->size.prod() * 4);
  SingleTimer timer;
  s.distrsolver.Run();
  const Scal tstep = timer.GetSeconds();

  if (s.step % 10 == 0) {
    std::cout << util::Format(
                     "step={:05} t={:.5f} dt={:.5f} diff={:.3f} vel={:.2f} "
                     "tstep={}ms",
                     s.step, s.time, s.dt, s.diffusion, s.velocity, tstep * 1e3)
              << std::endl;
  }
  CopyToCanvas(g_canvas->buf.data(), g_canvas->size[0], g_canvas->size[1]);
}

extern "C" {
Scal MulDiffusion(Scal factor) {
  auto& s = *g_state;
  s.diffusion *= factor;
  std::cout << util::Format("diff={}", s.diffusion) << std::endl;
  return s.diffusion;
}
Scal AddVelocityAngle(Scal add_deg) {
  auto& s = *g_state;
  const Vect v = s.velocity;
  Scal deg = std::atan2(v[1], v[0]) * 180. / M_PI;
  deg += add_deg;
  Scal rad = deg * M_PI / 180.;
  s.velocity = Vect(std::cos(rad), std::sin(rad)) * v.norm();
  std::cout << util::Format("vel={:.3f} angle={:.1f}", s.velocity, deg)
            << std::endl;
  return deg;
}
int Init() {
  auto& s = *g_state;
  s.to_init_field = true;
  return 0;
}
int TogglePause() {
  auto& s = *g_state;
  s.pause = !s.pause;
  return s.pause;
}

Vars var;

int SetMesh(int nx) {
  MPI_Comm comm = 0;
  std::stringstream conf;
  conf << GetDefaultConf();
  Subdomains<MIdx> sub(MIdx(nx), MIdx(nx), 1);
  conf << sub.GetConfig();
  Parser(var).ParseStream(conf);

  std::shared_ptr<State> new_state;
  try {
    new_state = std::make_shared<State>(comm, var);
  } catch (const std::exception& e) {
    std::cerr << FILELINE + "\nabort after throwing exception\n"
              << e.what() << '\n';
    std::terminate();
  } catch (...) {
    std::cerr << FILELINE + "\nabort after unknown exception\n";
    std::terminate();
  }

  if (g_state) {
    new_state->diffusion = g_state->diffusion;
    new_state->velocity = g_state->velocity;
  }
  g_state = new_state;
  std::cout << util::Format("mesh {}", MIdx(nx)) << std::endl;
  return 0;
}
int SetCanvas(int nx, int ny) {
  g_canvas = std::make_shared<Canvas>(MIdx(nx, ny));
  std::cout << util::Format("canvas {}", g_canvas->size) << std::endl;
  return 0;
}
}

int main() {
  FORCE_LINK(distr_local);
  FORCE_LINK(distr_native);

  SetCanvas(500, 500);
  SetMesh(32);
  emscripten_set_canvas_element_size(
      "#canvas", g_canvas->size[0], g_canvas->size[1]);
  emscripten_set_main_loop(main_loop, 30, 1);
}
