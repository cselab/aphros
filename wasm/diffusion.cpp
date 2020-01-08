#include <stdint.h>
#include <stdlib.h>
#include <emscripten.h>
#include <emscripten/html5.h>
#include <cmath>
#include <cstring>
#include <iostream>
#include <memory>

#include "geom/mesh.h"
#include "geom/rangemulti.h"
#include "solver/approx.h"
#include "solver/approx_eb.h"
#include "solver/cond.h"
#include "solver/solver.h"
#include "util/format.h"

const int dim = 2;
using Scal = double;
using M = MeshStructured<Scal, dim>;
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

M GetMesh(MIdx size) {
  const Rect<Vect> dom(Vect(0), Vect(1));
  const MIdx begin(0);
  const int halos = 2;
  return InitUniformMesh<M>(dom, begin, size, halos, true, true, size, 0);
}

struct Canvas {
  Canvas(MIdx size_) : size(size_), buf(size.prod()) {}
  MIdx size;
  std::vector<uint32_t> buf;
};

struct State {
  State(M&& m_) : m(std::move(m_)), fcu(m) {}
  M m;
  FieldCell<Scal> fcu;
  int step = 0;
  Scal time = 0;
  Scal diffusion = 0.01;
  Vect velocity{1., 0.5};
  bool pause = false;
};

std::shared_ptr<State> g_state;
std::shared_ptr<Canvas> g_canvas;

static float Clamp(float f) {
  return f < 0 ? 0 : f > 1 ? 1 : f;
}

static void Init(FieldCell<Scal>& fcu, const M& m) {
  for (auto c : m.AllCellsM()) {
    fcu[c] = (Vect(0.5).dist(c.center) < 0.2 ? 1 : 0);
  }
}

static void Step(
    FieldCell<Scal>& fcu, Scal dt, Scal diffusion, Vect velocity, const M& m) {
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

static void main_loop() {
  auto state = g_state;
  auto& s = *state;
  if (s.pause) {
    return;
  }
  const Scal h = s.m.GetCellSize()[0];
  const Scal dt =
      std::min(0.125 * h * h / s.diffusion, 0.25 * h / s.velocity.abs().max());
  Step(s.fcu, dt, s.diffusion, s.velocity, s.m);
  s.time += dt;
  ++s.step;
  if (s.step % 10 == 0) {
    std::cout << util::Format(
                     "step={:05} t={:.5f} dt={:.5f} diff={:.3f} vel={:.2f}",
                     s.step, s.time, dt, s.diffusion, s.velocity)
              << std::endl;
  }
  std::memset(g_canvas->buf.data(), 0, g_canvas->size.prod() * 4);
  Render(*g_canvas, s.fcu, s.m);
  CopyToCanvas(g_canvas->buf.data(), g_canvas->size[0], g_canvas->size[1]);
}

extern "C" {
double MulDiffusion(double factor) {
  auto& s = *g_state;
  s.diffusion *= factor;
  std::cout << util::Format("diff={}", s.diffusion) << std::endl;
  return s.diffusion;
}
double AddVelocityAngle(double add_deg) {
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
  Init(s.fcu, s.m);
  return 0;
}
int TogglePause() {
  auto& s = *g_state;
  s.pause = !s.pause;
  return s.pause;
}
int SetMesh(int nx) {
  g_state = std::make_shared<State>(GetMesh(MIdx(nx)));
  auto& s = *g_state;
  Init(s.fcu, s.m);
  std::cout << util::Format("mesh {}", s.m.GetGlobalSize()) << std::endl;
  return 0;
}
int SetCanvas(int nx, int ny) {
  g_canvas = std::make_shared<Canvas>(MIdx(nx, ny));
  std::cout << util::Format("canvas {}", g_canvas->size) << std::endl;
  return 0;
}
}

int main() {
  SetCanvas(500, 500);
  SetMesh(32);
  emscripten_set_canvas_element_size(
      "#canvas", g_canvas->size[0], g_canvas->size[1]);
  emscripten_set_main_loop(main_loop, 20, 1);
}
