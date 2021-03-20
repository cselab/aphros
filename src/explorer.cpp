// Created by Petr Karnakov on 20.03.2021
// Copyright 2021 ETH Zurich

#include <stdint.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <memory>
#include <sstream>

#include "distr/distrbasic.h"
#include "distr/distrsolver.h"
#include "geom/rangemulti.h"
#include "kernel/hydro.h"
#include "kernel/kernelmeshpar.h"
#include "parse/argparse.h"
#include "util/distr.h"
#include "util/git.h"
#include "util/hydro_post.h"
#include "util/timer.h"
#include "util/visual.h"

const char* g_base_conf =
#include "explorer.inc"
;

using M = MeshStructured<double, 2>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using MIdx = typename M::MIdx;
using util::Canvas;
using util::CanvasView;
using util::MIdx2;
using util::Float3;
using util::Byte3;
using util::Pixel;

template <class M>
static MIdx2 GetViewCoords(Vect x, const CanvasView& view, const M& m) {
  MIdx2 w =
      view.start + MIdx2(x * Vect(view.end - view.start) / m.GetGlobalLength());
  w = w.max(view.start).min(view.end - MIdx2(1));
  return w;
}

struct State {
  State(MPI_Comm comm, Vars& var, typename Hydro<M>::Par par_)
      : par(par_), distrsolver(comm, var, par) {}
  typename Hydro<M>::Par par;
  DistrSolver<M, Hydro<M>> distrsolver;
  bool render = true;
  int frame = 0;
  std::vector<std::array<MIdx, 2>> lines; // interface lines
  std::vector<std::array<MIdx, 2>> lines_eb; // embed lines
};

std::shared_ptr<State> g_state;
std::shared_ptr<Canvas> g_canvas;
std::shared_ptr<CanvasView> g_view;
std::string g_extra_config;
Vars g_var;
bool g_exit = false;
bool g_is_root = false;
bool g_verbose = false;

static void StepCallback(void*, Hydro<M>* hydro) {
  if (!g_canvas) {
    return;
  }
  auto state = g_state;
  auto& s = *state;
  auto& m = hydro->m;
  auto& var = hydro->var;
  auto sem = m.GetSem();

  if (sem()) {
    if (m.IsRoot() && hydro->st_.step % var.Int("report_step_every", 1) == 0) {
      auto names = GetWords(var.String("print_vars", ""));
      for (auto name : names) {
        auto type = var.GetTypeName(name);
        if (!type.empty()) {
          std::cerr << name << '=' << var.GetStr(type, name) << ' ';
        }
      }
      if (!names.empty()) {
        std::cerr << std::endl;
      }
      if (auto* str = var.String.Find("print_string")) {
        std::cerr << *str << std::endl;
      }
    }

    if (s.render) {
      auto& canvas = *g_canvas;
      using U = util::Visual<M>;
      util::CanvasView view(
          canvas.size, MIdx(0), canvas.size, canvas.buf.data());
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
                    GetViewCoords(poly[0], *g_view, m),
                    GetViewCoords(poly[1], *g_view, m),
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
                GetViewCoords(poly[0], *g_view, m),
                GetViewCoords(poly[1], *g_view, m),
            });
          }
        }
      }
    }
  }
}

static void DrawLines(
    CanvasView& view, const std::vector<std::array<MIdx, 2>>& lines,
    Float3 color, float width) {
  for (auto line : lines) {
    auto draw = [&](bool q) {
      if (line[1][q] < line[0][q]) {
        std::swap(line[0], line[1]);
      }
      const int dx = line[1][q] - line[0][q];
      const int dy = std::abs(line[1][!q] - line[0][!q]);
      const int sy = line[1][!q] > line[0][!q] ? 1 : -1;
      int y = line[0][!q];
      int err = 2 * std::abs(dy) - dx;
      for (int x = line[0][q]; x < line[1][q]; ++x) {
        if (err <= 0) {
          err += 2 * dy;
        } else {
          err += 2 * (dy - dx);
          y += sy;
        }
        for (int hy = -width; hy < (width + 1); ++hy) {
          if (x >= view.start[q] && x < view.end[q] && y >= view.start[!q] &&
              y < view.end[!q]) {
            Vect v;
            v[q] = x;
            v[!q] = y + hy;
            float f;
            f = width * 0.5 -
                Reconst<Scal>::GetNearest(v, Vect(line[0]), Vect(line[1]))
                    .dist(v);
            f = Reconst<Scal>::Clip(f * 0.75, 0, 1);

            const Pixel p0 = (!q ? view(x, y + hy) : view(y + hy, x));
            Byte3 b0(p0 & 0xFF, (p0 >> 8) & 0xFF, (p0 >> 16) & 0xFF);
            Float3 c0 = Float3(b0) / 255;
            Float3 c = color * f + c0 * (1 - f);
            util::Byte3 b(c * 255);
            util::Pixel p =
                0xFF000000 | (b[0] << 0) | (b[1] << 8) | (b[2] << 16);
            (!q ? view(x, y + hy) : view(y + hy, x)) = p;
          }
        }
      }
    };
    const bool q =
        std::abs(line[1][0] - line[0][0]) < std::abs(line[1][1] - line[0][1]);
    draw(q);
    draw(!q);
  }
}

static void main_loop() {
  if (!g_state) {
    return;
  }

  auto state = g_state;
  auto& s = *state;

  g_verbose = g_var.Int["VERBOSE"];
  g_view->Clear(0xFFFFFFFF);
  for (int i = 0; i < g_var.Int("steps_per_frame", 1); ++i) {
    s.render = (i == 0);
    s.distrsolver.Run();
  }
  if (g_is_root) {
    DrawLines(*g_view, s.lines, Float3(0), g_var.Int["line_width"]);
    DrawLines(*g_view, s.lines_eb, Float3(0), g_var.Int["line_width"]);
  }
  if (g_is_root) {
    if (g_var.Int["stdout"]) {
      util::WritePpm(std::cout, *g_view);
    } else {
      const std::string path = util::Format("a_{:06d}.ppm", s.frame);
      util::WritePpm(path, *g_view);
      if (g_verbose) {
        std::cerr << path << std::endl;
      }
    }
  }
  ++s.frame;
  const int num_frames = g_var.Int["num_frames"];
  if (num_frames >= 0 && s.frame > num_frames) {
    g_exit = true;
  }
}

static std::string GetBaseConfig() {
  return R"EOF(
set int num_frames -1

set int steps_per_frame 1
set int return_after_each_step 1
set string visual

set int line_width 5

set string Cred 1 0.12 0.35
set string Cgreen 0 0.8 0.42
set string Cblue 0 0.6 0.87
set string Cpurple 0.686 0.345 0.729
set string Cyellow 1 0.776 0.118
set string Corange 0.949 0.522 0.133
set string Cwhite 1 1 1
set string Cblack 0 0 0
set string Cgray 0.627 0.694 0.729
)EOF";
}

static void SetMesh(const MpiWrapper& mpi, int nx, int bx) {
  std::stringstream conf;
  conf << GetDefaultConf();
  const MIdx meshsize(nx);
  const MIdx blocksize(bx > 0 ? bx : nx);
  const Subdomains<MIdx> sub(meshsize, blocksize, mpi.GetCommSize());
  conf << g_base_conf;
  conf << GetBaseConfig();
  conf << g_extra_config << '\n';
  conf << sub.GetConfig();
  Parser(g_var).ParseStream(conf);

  std::shared_ptr<State> new_state;
  Hydro<M>::Par par;
  par.step_callback = StepCallback;
  new_state = std::make_shared<State>(mpi.GetComm(), g_var, par);

  g_state = new_state;
  if (g_is_root && g_verbose) {
    std::cerr << util::Format("mesh {}", MIdx(nx)) << std::endl;
  }
}

static void SetCanvas(int nx, int ny) {
  g_canvas = std::make_shared<Canvas>(MIdx(nx, ny));
  g_view = std::make_shared<CanvasView>(*g_canvas);
  if (g_is_root && g_verbose) {
    std::cerr << util::Format("canvas {}", g_canvas->size) << std::endl;
  }
}

void FindReplace(std::string& s, std::string find, std::string replace) {
  auto pos = s.find(find);
  while (pos != std::string::npos) {
    s.replace(pos, find.size(), replace);
    pos = s.find(find);
  }
}

int main(int argc, const char** argv) {
  FORCE_LINK(distr_native);
  FORCE_LINK(init_contang);
  FORCE_LINK(init_vel);

  MpiWrapper mpi(&argc, &argv);
  g_is_root = mpi.IsRoot();

  ArgumentParser parser(
      "Sharpens the image using PLIC advection", mpi.IsRoot());
  parser.AddSwitch({"--verbose", "-v"}).Help("Report steps");
  parser.AddSwitch("--stdout").Help("Output ppm to STDOUT");
  parser.AddVariable<std::string>("--extra", "")
      .Help("Extra configuration (commands 'set ... ')");
  parser.AddVariable<int>({"--steps_per_frame", "-s"}, 10)
      .Help("Number of steps per frame");
  parser.AddVariable<int>({"--frames", "-f"}, 10).Help("Number of frames");
  parser.AddVariable<int>({"--mesh", "-m"}, 64).Help("Mesh size");
  parser.AddVariable<int>({"--block", "-b"}, 0)
      .Help("Block size, defaults to full mesh");
  parser.AddVariable<int>({"--canvas", "-c"}, 512).Help("Canvas size");
  parser.AddSwitch("--logo").Help("Print logo");
  parser.AddSwitch("--exit").Help("Exit before reading the configuration");
  parser.AddVariable<std::string>("config", "-")
      .Help("Path to input configuration file, - to read from STDIN");

  auto args = parser.ParseArgs(argc, argv);
  if (const int* p = args.Int.Find("EXIT")) {
    return *p;
  }

  if (args.Int["logo"] && mpi.IsRoot()) {
    std::cerr << GetLogo();
  }
  if (args.Int["exit"]) {
    return 0;
  }

  g_verbose = args.Int["verbose"];

  g_extra_config += util::Format(
      R"EOF(
set int num_frames {}
set int steps_per_frame {}
set int VERBOSE {}
set int silent {}
set int stdout {}
{}
)EOF",
      args.Int["frames"], args.Int["steps_per_frame"], args.Int["verbose"],
      !g_verbose, args.Int["stdout"], args.String["extra"]);

  if (!args.Int["verbose"]) {
    g_extra_config += R"EOF(
set int report_step_every 1000000
  )EOF";
  }

  const std::string confpath = args.String["config"];
  std::stringstream buf;
  if (confpath == "-") {
    buf << std::cin.rdbuf();
  } else {
    std::ifstream f(confpath);
    fassert(f.good(), "Can't open file '" + confpath + "' for reading");
    buf << f.rdbuf();
  }
  std::string bufstr = buf.str();
  FindReplace(bufstr, "FROMSLIDER", "");
  g_extra_config += bufstr;

  SetCanvas(args.Int["canvas"], args.Int["canvas"]);
  SetMesh(mpi, args.Int["mesh"], args.Int["block"]);

  while (!g_exit) {
    main_loop();
  }
}
