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
#include "util/git.h"
#include "util/hydro_post.h"
#include "util/timer.h"
#include "util/visual.h"

static constexpr int kScale = 1;

using M = MeshCartesian<double, 2>;
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

  bool to_spawn = false;
  Vect spawn_c;
  Scal spawn_r;
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
  auto sem = m.GetSem();

  if (sem()) {
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
      util::CanvasView view(
          canvas.size, MIdx(0), canvas.size, canvas.buf.data());
      using util::Float3;
      FieldCell<Float3> fc_color(m, Float3(1));
      std::stringstream str_entries(var.String["visual"]);
      auto entries = util::ParseEntries(str_entries);
      auto get_field = [&](std::string name) -> FieldCell<Scal> {
        return HydroPost<M>::GetField(hydro, name, m);
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

  if (sem() && s.to_spawn) {
    auto spawn = [&](auto& fcu, auto& fccl, const auto& meb) {
      FieldCell<Scal> fc_add(m);
      GetCircle(fc_add, s.spawn_c, s.spawn_r, m);
      for (auto c : meb.Cells()) {
        if (fc_add[c] > fcu[c]) {
          fcu[c] = fc_add[c];
          fccl[c] = 0;
        }
      }
    };
    if (auto vof = dynamic_cast<Vof<M>*>(hydro->as_.get())) {
      vof->AddModifier(
          [spawn](FieldCell<Scal>& fcu, FieldCell<Scal>& fccl, const M& m) { //
            spawn(fcu, fccl, m);
          });
    }
    if (auto vof = dynamic_cast<Vofm<M>*>(hydro->as_.get())) {
      vof->AddModifier([spawn](
                           const Multi<FieldCell<Scal>*>& fcu,
                           const Multi<FieldCell<Scal>*>& fccl, GRange<size_t>,
                           const M& m) { //
        spawn(*fcu[0], *fccl[0], m);
      });
    }
    if (auto vof = dynamic_cast<Vof<Embed<M>>*>(hydro->as_.get())) {
      vof->AddModifier([spawn](
                           FieldCell<Scal>& fcu, FieldCell<Scal>& fccl,
                           const Embed<M>& eb) { //
        spawn(fcu, fccl, eb);
      });
    }
    if (auto vof = dynamic_cast<Vofm<Embed<M>>*>(hydro->as_.get())) {
      vof->AddModifier([spawn](
                           const Multi<FieldCell<Scal>*>& fcu,
                           const Multi<FieldCell<Scal>*>& fccl, GRange<size_t>,
                           const Embed<M>& eb) { //
        spawn(*fcu[0], *fccl[0], eb);
      });
    }
  }
  if (sem() && m.IsRoot()) {
    s.to_spawn = false;
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
  auto state = g_state;
  auto& s = *state;
  s.to_spawn = true;
  s.spawn_c = Vect(x, y);
  s.spawn_r = r;
  std::cout << util::Format(
                   "sphere {:.3f} {:.3f} 0 {:.3f}", s.spawn_c[0], s.spawn_c[1],
                   s.spawn_r)
            << std::endl;
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

  g_var.String.Set("git_rev", GetGitRev());
  g_var.String.Set("git_diff", GetGitDiff());
  g_var.String.Set("git_msg", GetGitMsg());
  g_var.String.Set("logo", GetLogo());

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
