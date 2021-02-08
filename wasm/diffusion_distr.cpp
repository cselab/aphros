#include <cmath>
#include <cstring>
#include <emscripten.h>
#include <emscripten/html5.h>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdint.h>
#include <stdlib.h>

#include "distr/distrbasic.h"
#include "distr/distrsolver.h"
#include "geom/mesh.h"
#include "geom/rangemulti.h"
#include "solver/approx.h"
#include "solver/approx_eb.h"
#include "solver/cond.h"
#include "solver/solver.h"
#include "solver/vof.h"
#include "util/format.h"
#include "util/module.h"
#include "util/timer.h"

using M = MeshStructured<double, 2>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using MIdx = typename M::MIdx;

void CopyToCanvas(uint32_t *ptr, int w, int h) {
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

  Solver(Vars &var_, const BlockInfoProxy &proxy, Par &par)
      : KernelMeshPar<M, Par>(var_, proxy, par), fcu(m, 0), fc_src(m, 0),
        ff_flux(m, 0) {}
  void Run() override;

protected:
  using P::m;
  using P::par_;
  using P::var;

private:
  FieldCell<Scal> fcu;
  FieldFace<Scal> ff_flux;
  FieldCell<Scal> fc_src;
  MapEmbed<BCondAdvection<Scal>> bc_vof;
  std::unique_ptr<Vof<M>> vof;
};

struct Canvas {
  Canvas(MIdx size_) : size(size_), buf(size.prod()) {}
  MIdx size;
  std::vector<uint32_t> buf;
};

M GetMesh(MIdx size) {
  const Rect<Vect> dom(Vect(0), Vect(1));
  const MIdx begin(0);
  const int halos = 2;
  return InitUniformMesh<M>(dom, begin, size, halos, true, true, size, 0);
}

struct State {
  State(MPI_Comm comm, Vars &var)
      : m(GetMesh(MIdx(8))), distrsolver(comm, var, par) {}
  M m;
  typename Solver::Par par;
  DistrSolver<M, Solver> distrsolver;
  int step = 0;
  Scal time = 0;
  Scal diffusion = 0.01;
  Vect velocity{1., 0.5};
  bool pause = false;
  bool to_init_field = true;
  bool to_init_solver = true;
};

std::shared_ptr<State> g_state;
std::shared_ptr<Canvas> g_canvas;

static Scal Clamp(Scal f) { return f < 0 ? 0 : f > 1 ? 1 : f; }

static void Init(FieldCell<Scal> &fcu, const M &m) {
  for (auto c : m.AllCellsM()) {
    fcu[c] = (Vect(0.5).dist(c.center) < 0.2 ? 1 : 0);
  }
}

static void Step(FieldCell<Scal> &fcu, Scal dt, Scal diffusion, Vect velocity,
                 const M &m) {
  using UEB = UEmbed<M>;
  const auto ffg = UEB::Gradient(fcu, {}, m);
  const auto fcg = UEB::AverageGradient(ffg, m);
  FieldFace<Scal> ffv(m);
  for (auto f : m.Faces()) {
    ffv[f] = velocity.dot(m.GetSurface(f));
  }
  auto ffu = UEB::InterpolateUpwind(fcu, {}, ConvSc::superbee, fcg, ffv, m);

  for (auto c : m.CellsM()) {
    Scal sum = 0;
    for (auto q : m.Nci(c)) {
      auto f = c.face(q);
      sum += diffusion * ffg[f] * c.outward_factor(q) * f.area;
      sum -= ffv[f] * ffu[f] * c.outward_factor(q);
    }
    fcu[c] += dt * sum / c.volume;
  }
}

static void Render(Canvas &canvas, const FieldCell<Scal> &fcu, const M &m) {
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
  auto &s = *state;
  if (sem("init")) {
    s.to_init_solver = false;
    if (s.to_init_field) {
      Init(fcu, m);
    }
    for (auto f : m.Faces()) {
      ff_flux[f] = s.velocity.dot(m.GetSurface(f));
    }
  }
  if (sem("init_solver")) {
    auto p = ParsePar<Vof<M>>()(var);
    const FieldCell<Scal> fccl(m, 0);
    const Scal dt = 0.25 * h / s.velocity.abs().max();
    vof.reset(new Vof<M>(m, m, fcu, fccl, bc_vof, &ff_flux, &fc_src, 0, dt, p));
  }
  if (sem.Nested("start")) {
    vof->StartStep();
  }
  if (sem.Nested("iter")) {
    vof->MakeIteration();
  }
  if (sem.Nested("finish")) {
    vof->FinishStep();
  }
  if (sem()) {
    s.to_init_field = false;
    const Scal h = m.GetCellSize()[0];
    const Scal dt = std::min(0.125 * h * h / s.diffusion,
                             0.25 * h / s.velocity.abs().max());

    Step(fcu, dt, s.diffusion, s.velocity, m);
    //Render(*g_canvas, fcu, m);
    Render(*g_canvas, vof.GetField(), m);
    m.Comm(&fcu);

    if (m.IsRoot()) {
      s.time += dt;
      ++s.step;
    }
  }
}

static void main_loop() {
  auto state = g_state;
  auto &s = *state;
  if (s.pause) {
    return;
  }
  std::memset(g_canvas->buf.data(), 0, g_canvas->size.prod() * 4);
  SingleTimer timer;
  s.distrsolver.Run();
  const Scal tstep = timer.GetSeconds();

  if (s.step % 10 == 0) {
    auto dt = 0;
    std::cout << util::Format(
                     "step={:05} t={:.5f} dt={:.5f} diff={:.3f} vel={:.2f} "
                     "tstep={}ms",
                     s.step, s.time, dt, s.diffusion, s.velocity, tstep * 1e3)
              << std::endl;
  }
  CopyToCanvas(g_canvas->buf.data(), g_canvas->size[0], g_canvas->size[1]);
}

extern "C" {
Scal MulDiffusion(Scal factor) {
  auto &s = *g_state;
  s.diffusion *= factor;
  std::cout << util::Format("diff={}", s.diffusion) << std::endl;
  return s.diffusion;
}
Scal AddVelocityAngle(Scal add_deg) {
  auto &s = *g_state;
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
  auto &s = *g_state;
  s.to_init_field = true;
  return 0;
}
int TogglePause() {
  auto &s = *g_state;
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
  } catch (const std::exception &e) {
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

  SetCanvas(500, 500);
  SetMesh(32);
  emscripten_set_canvas_element_size("#canvas", g_canvas->size[0],
                                     g_canvas->size[1]);
  emscripten_set_main_loop(main_loop, 30, 1);
}
