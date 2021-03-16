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
#include "util/timer.h"
#include "util/visual.h"

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
  bool render = true;
  bool pause = false;
  std::vector<std::array<MIdx, 2>> lines; // interface lines
  std::vector<std::array<MIdx, 2>> lines_eb; // embed lines
};

std::shared_ptr<State> g_state;
std::shared_ptr<Canvas> g_canvas;
std::string g_extra_config;
Vars g_var;

void StepCallback(void*, Hydro<M>* hydro) {
  if (!g_canvas) {
    return;
  }
  auto state = g_state;
  auto& s = *state;
  auto& m = hydro->m;
  auto& var = hydro->var;

  if (hydro->st_.step % var.Int("report_step_every", 1) == 0) {
    auto names = GetWords(var.String("print_vars", ""));
    for (auto name : names) {
      auto type = var.GetTypeName(name);
      if (!type.empty()) {
        std::cout << name << '=' << var.GetStr(type, name) << ' ';
      }
    }
    if (!names.empty()) {
      std::cout << std::endl;
    }
    if (auto* str = var.String.Find("print_string")) {
      std::cout << *str << std::endl;
    }
  }

  if (s.render) {
    auto& canvas = *g_canvas;
    using U = util::Visual<M>;
    typename U::CanvasView view(
        canvas.size, MIdx(0), canvas.size, canvas.buf.data());
    using Float3 = typename U::Float3;
    FieldCell<Float3> fc_color(m, Float3(1));
    std::stringstream str_entries(var.String["visual"]);
    auto entries = U::ParseEntries(str_entries);
    auto get_field = [&](std::string name) -> FieldCell<Scal> {
      if (name == "p" || name == "pressure") {
        return hydro->fs_->GetPressure();
      }
      if (name == "omz" || name == "vorticity") {
        if (hydro->eb_) {
          return GetVortScal(
              hydro->fs_->GetVelocity(), hydro->fs_->GetVelocityCond(),
              *hydro->eb_);
        }
        return GetVortScal(
            hydro->fs_->GetVelocity(), hydro->fs_->GetVelocityCond(), m);
      }
      if (name == "ebvf" || name == "embed fraction") {
        FieldCell<Scal> fc(m, 1);
        if (hydro->eb_) {
          auto& eb = *hydro->eb_;
          for (auto c : m.SuCells()) {
            fc[c] = eb.GetVolumeFraction(c);
          }
        }
        return fc;
      }
      if (name == "vf" || name == "volume fraction") {
        return hydro->as_->GetField();
      }
      if (name == "vmagn" || name == "velocity magnitude") {
        FieldCell<Scal> fc(m, 0);
        for (auto c : m.SuCells()) {
          fc[c] = hydro->fs_->GetVelocity()[c].norm();
        }
        return fc;
      }
      if (name == "vx" || name == "velocity x") {
        FieldCell<Scal> fc(m, 0);
        for (auto c : m.SuCells()) {
          fc[c] = hydro->fs_->GetVelocity()[c][0];
        }
        return fc;
      }
      if (name == "vy" || name == "velocity y") {
        FieldCell<Scal> fc(m, 0);
        for (auto c : m.SuCells()) {
          fc[c] = hydro->fs_->GetVelocity()[c][1];
        }
        return fc;
      }
      if (m.IsRoot()) {
        std::cerr << "Unknown field '" + name + "'\n";
      }
      return FieldCell<Scal>(m, 0);
    };
    U::RenderEntriesToField(fc_color, entries, get_field, m);
    if (var.Int("visual_interpolate", 1) && m.GetGlobalSize() < view.size) {
      U::RenderToCanvasBilinear(view, fc_color, m);
    } else {
      U::RenderToCanvasNearest(view, fc_color, m);
    }

    // Render interface lines
    if (m.IsRoot()) {
      s.lines.clear();
    }
    if (var.Int("visual_lines", 1)) {
      auto h = m.GetCellSize();
      const auto& plic = hydro->as_->GetPlic();
      for (auto l : plic.layers) {
        const auto& fci = *plic.vfci[l];
        const auto& fcn = *plic.vfcn[l];
        const auto& fca = *plic.vfca[l];
        for (auto c : m.Cells()) {
          if (fci[c]) {
            const auto poly =
                Reconst<Scal>::GetCutPoly(m.GetCenter(c), fcn[c], fca[c], h);
            if (poly.size() == 2) {
              s.lines.push_back({
                  GetCanvasCoords(poly[0], *g_canvas, m),
                  GetCanvasCoords(poly[1], *g_canvas, m),
              });
            }
          }
        }
      }
    }
    // Render embed lines
    if (m.IsRoot()) {
      s.lines_eb.clear();
    }
    if (var.Int("visual_lines_eb", 1) && hydro->eb_) {
      auto& eb = *hydro->eb_;
      auto h = m.GetCellSize();
      for (auto c : eb.CFaces()) {
        const auto poly = Reconst<Scal>::GetCutPoly(
            m.GetCenter(c), eb.GetNormal(c), eb.GetAlpha(c), h);
        if (poly.size() == 2) {
          s.lines_eb.push_back({
              GetCanvasCoords(poly[0], *g_canvas, m),
              GetCanvasCoords(poly[1], *g_canvas, m),
          });
        }
      }
    }
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

  for (int i = 0; i < g_var.Int("steps_per_frame", 1); ++i) {
    s.render = (i == 0);
    s.distrsolver.Run();
  }

  CopyToCanvas(g_canvas->buf.data(), g_canvas->size[0], g_canvas->size[1]);
  EM_ASM_({ Draw(); });
}

static std::string GetBaseConfig() {
  return R"EOF(
include conf/base.conf
include conf/std.conf

set int steps_per_frame 1
set int return_after_each_step 1
set string visual
)EOF";
}

extern "C" {
void Spawn(float x, float y, float r) {
  if (!g_state) {
    return;
  }
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

  if (g_var.Int("verbose_runtime_config", 0)) {
    std::cout << conf.str() << std::endl;
  }

  Parser(g_var).ParseStream(conf);
}

const char* GetConfigString(const char* name) {
  auto* ptr = g_var.String.Find(name);
  return ptr ? ptr->c_str() : "";
}

double GetConfigDouble(const char* name) {
  auto* ptr = g_var.Double.Find(name);
  return ptr ? *ptr : GetNan<double>();
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
int GetLines(int embed, uint16_t* data, int max_size) {
  if (!g_state) {
    return 0;
  }
  auto state = g_state;
  auto& s = *state;
  int i = 0;
  const auto& lines = embed ? s.lines_eb : s.lines;
  for (auto p : lines) {
    if (i + 3 >= max_size) {
      break;
    }
    data[i] = p[0][0];
    data[i + 1] = p[0][1];
    data[i + 2] = p[1][0];
    data[i + 3] = p[1][1];
    i += 4;
  }
  return i;
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

  emscripten_set_main_loop(main_loop, 20, 0);
}
