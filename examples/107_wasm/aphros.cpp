// Created by Petr Karnakov on 04.03.2021
// Copyright 2021 ETH Zurich

#include <stdint.h>
#include <stdlib.h>
#include <emscripten.h>
#include <emscripten/html5.h>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <memory>
#include <sstream>

#include "distr/distrbasic.h"
#include "distr/distrsolver.h"
#include "geom/rangemulti.h"
#include "kernel/hydro.h"
#include "kernel/kernelmeshpar.h"
#include "solver/approx_eb.ipp"
#include "util/timer.h"

static constexpr int kScale = 1;

using M = MeshStructured<double, 2>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using MIdx = typename M::MIdx;

#include "common.h"

struct State {
  State(MPI_Comm comm, Vars& var, typename Hydro<M>::Par par_)
      : par(par_), distrsolver(comm, var, par) {}
  typename Hydro<M>::Par par;
  DistrSolver<M, Hydro<M>> distrsolver;
  bool pause = false;
  std::vector<std::array<MIdx, 2>> lines; // interface lines
};

std::shared_ptr<State> g_state;
std::shared_ptr<Canvas> g_canvas;
std::string g_extra_config;
Vars g_var;

void StepCallback(void*, Hydro<M>* hydro) {
  if (!g_canvas) {
    return;
  }
  auto& m = hydro->m;
  auto& fs = hydro->fs_;

  auto& canvas = *g_canvas;

  const auto msize = m.GetGlobalSize();
  using Vect3 = generic::Vect<Scal, 3>;
  using MIdx3 = generic::Vect<unsigned char, 3>;
  const Scal vispressure = hydro->var.Double("vispressure", 0);

  auto get_color = [&](IdxCell c) -> Vect3 {
    auto& fcp = hydro->as_->GetField();
    Vect3 q(0);
    const auto p = (fcp[c] * vispressure) * 2 - 1;
    const auto f = Clamp(std::abs(p));
    if (p > 0) {
      q = interp::Linear(f, Vect3(1, 1, 1), Vect3(1, 0.12, 0.35));
    } else {
      q = interp::Linear(f, Vect3(1, 1, 1), Vect3(0, 0.6, 0.87));
    }
    return q;
  };

  for (auto c : m.CellsM()) {
    const MIdx w(c);
    const MIdx start = w * canvas.size / msize;
    const MIdx end = (w + MIdx(1)) * canvas.size / msize;
    auto q = get_color(c);
    const MIdx3 mq(q * 255);
    const uint32_t v =
        0xff000000 | (mq[0] << 0) | (mq[1] << 8) | (mq[2] << 16);
    for (int y = start[1]; y < end[1]; y++) {
      for (int x = start[0]; x < end[0]; x++) {
        canvas.buf[canvas.size[0] * (canvas.size[1] - y - 1) + x] = v;
      }
    }
  }

  if (m.IsRoot()) {
    Scal ekin = 0;
    for (auto c : m.CellsM()) {
      ekin += fs->GetVelocity()[c].sqrnorm() * c.volume;
    }
    std::cerr << util::Format("ekin={}\n", ekin);
  }
}

static void main_loop() {
  if (!g_state) {
    return;
  }
  auto state = g_state;
  auto& s = *state;
  if (s.pause) {
    return;
  }

  std::memset(g_canvas->buf.data(), 0, g_canvas->size.prod() * 4);

  s.distrsolver.Run();

  CopyToCanvas(g_canvas->buf.data(), g_canvas->size[0], g_canvas->size[1]);
  EM_ASM_({ Draw(); });
}

static std::string GetBaseConfig() {
  return R"EOF(
set string linsolver_symm conjugate
set double hypre_symm_tol 1e-3
set int hypre_symm_maxiter 100
set int hypre_periodic_x 0
set int hypre_periodic_y 0
set double dtmax 0.1
set double rho1 1
set double rho2 1
set double mu1 0.001
set double mu2 0.001
set double sigma 0
set vect gravity 0 -1
set double topvel 0
set int sharpen 0
set int layers 3
set double coalth 1.5
set double sharpen_cfl 0.1
set int coal 1

set int nsteps 1
set int reportevery 10

set double cfl 1
set double cflvis 0.5
set double cflsurf 2

set double visvel 0
set double visvf 1
set double visvort 0
set int visinterp 0

set int part 1
set double part_h 4
set double part_relax 0.5
set int part_np 7
set double part_segcirc 0
set int part_dn 0
set int part_dump_fr 1
set int part_ns 1
set double part_tol 1e-4
set int part_itermax 10
set int part_verb 0
set int dim 2
set int vtkbin 1
set int vtkmerge 1

include conf/coal/a.conf

set int verbose_time 0
set int return_after_each_step 1


set double tmax 0.5
set int verbose_stages 0
set int output 0
)EOF";
}

extern "C" {
void Spawn(float x, float y, float r) {
  if (!g_state) {
    return;
  }
  auto state = g_state;
  auto& s = *state;
  std::cout << util::Format("action") << std::endl;
}
int TogglePause() {
  if (!g_state) {
    return 0;
  }
  auto& s = *g_state;
  s.pause = !s.pause;
  return s.pause;
}
void SetExtraConfig(const char* conf) {
  g_extra_config = conf;
}

void SetRuntimeConfig(const char* str) {
  if (!g_state) {
    return;
  }
  std::stringstream conf(str);
  Parser(g_var).ParseStream(conf);
}

void SetMesh(int nx) {
  MPI_Comm comm = 0;
  std::stringstream conf;
  conf << GetDefaultConf();
  Subdomains<MIdx> sub(MIdx(nx), MIdx(nx), 1);
  conf << GetBaseConfig();
  conf << g_extra_config << '\n';
  conf << sub.GetConfig();
  Parser(g_var).ParseStream(conf);

  std::shared_ptr<State> new_state;
  Hydro<M>::Par par;
  par.step_callback = StepCallback;
  new_state = std::make_shared<State>(comm, g_var, par);

  g_state = new_state;
  std::cout << util::Format("mesh {}", MIdx(nx)) << std::endl;
}
void SetCanvas(int nx, int ny) {
  g_canvas = std::make_shared<Canvas>(MIdx(nx, ny));
  std::cout << util::Format("canvas {}", g_canvas->size) << std::endl;
}
}

int main() {
  FORCE_LINK(distr_local);
  FORCE_LINK(distr_native);
  FORCE_LINK(init_contang);
  FORCE_LINK(init_vel);

  aphros_SetErrorHandler(ErrorHandler);

  SetCanvas(512, 512);
  emscripten_set_canvas_element_size(
      "#canvas", g_canvas->size[0] * kScale, g_canvas->size[1] * kScale);

  emscripten_set_main_loop(main_loop, 20, 1);
}
