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
#include "linear/linear.h"
#include "parse/curv.h"
#include "parse/vof.h"
#include "solver/approx.h"
#include "solver/approx_eb.h"
#include "solver/cond.h"
#include "solver/curv.h"
#include "solver/proj.h"
#include "solver/solver.h"
#include "solver/vof.h"
#include "util/format.h"
#include "util/hydro.ipp"
#include "util/linear.h"
#include "util/module.h"
#include "util/timer.h"

using M = MeshStructured<double, 2>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using MIdx = typename M::MIdx;

void CopyToCanvas(uint32_t* buf, int w, int h) {
  EM_ASM_(
      {
        let data = Module.HEAPU8.slice($0, $0 + $1 * $2 * 4);
        let ctx = Module['canvas'].getContext('2d');
        let image = ctx.getImageData(0, 0, $1, $2);
        image.data.set(data);
        ctx.putImageData(image, 0, 0);
      },
      buf, w, h);
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
      : KernelMeshPar<M, Par>(var_, proxy, par)
      , fc_src(m, 0)
      , fc_rho(m, 1)
      , fc_mu(m, 0)
      , fc_curv(m, 0)
      , fc_force(m, Vect(0))
      , ff_bforce(m, 0) {}
  void Run() override;

 protected:
  using P::m;
  using P::par_;
  using P::var;

 private:
  FieldCell<Scal> fc_src;
  FieldCell<Scal> fc_rho;
  FieldCell<Scal> fc_mu;
  FieldCell<Scal> fc_curv;
  FieldCell<Vect> fc_force;
  FieldEmbed<Scal> ff_bforce;
  MapEmbed<BCondAdvection<Scal>> bc_vof;
  MapEmbed<BCondFluid<Vect>> bc_fluid;
  MapCell<std::shared_ptr<CondCellFluid>> mc_velcond;
  std::unique_ptr<Vof<M>> vof;
  std::unique_ptr<Proj<M>> fluid;
  Scal maxvel;
  typename PartStrMeshM<M>::Par psm_par;
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
  Scal dt = 0.001;
  bool pause = false;
  bool to_init_field = true;
  bool to_init_solver = true;
  std::vector<std::array<MIdx, 2>> lines; // interface lines

  bool to_spawn = false;
  Vect spawn_c;
  Scal spawn_r;
};

std::shared_ptr<State> g_state;
std::shared_ptr<Canvas> g_canvas;
std::string g_extra_config;
Vars g_var;

static Scal Clamp(Scal f, Scal min, Scal max) {
  return f < min ? min : f > max ? max : f;
}

static Scal Clamp(Scal f) {
  return Clamp(f, 0, 1);
}

void GetCircle(FieldCell<Scal>& fcu, Vect c, Scal r,  const M& m) {
  Vars par;
  par.String.Set("init_vf", "circlels");
  par.Vect.Set("circle_c", c);
  par.Double.Set("circle_r", r);
  par.Int.Set("dim", 2);
  auto func = CreateInitU<M>(par, false);
  func(fcu, m);
}

static MIdx GetCanvasCoords(Vect x, const Canvas& canvas, const M& m) {
  MIdx w = MIdx(x * Vect(canvas.size) / m.GetGlobalLength());
  w = w.max(MIdx(0)).min(canvas.size - MIdx(1));
  w[1] = canvas.size[1] - w[1] - 1;
  return w;
}

static void RenderField(
    Canvas& canvas, const FieldCell<Scal>& fcu, const FieldCell<Vect>& fcvel,
    const FieldCell<Scal>& fcp, const M& m) {
  const auto msize = m.GetGlobalSize();
  const Scal visvel = g_var.Double("visvel", 0);
  const Scal visvf = g_var.Double("visvf", 0);
  for (auto c : m.CellsM()) {
    const MIdx w(c);
    const MIdx start = w * canvas.size / msize;
    const MIdx end = (w + MIdx(1)) * canvas.size / msize;
    unsigned char qr = 0;
    unsigned char qg = 0;
    unsigned char qb = 0;
    if (visvel) {
      qg = qb = Clamp(visvel * fcvel[c].norm()) * 255;
    }
    qr = ~qr;
    qg = ~qg;
    qb = ~qb;
    if (visvf) {
      const auto f = Clamp(1 - visvf * fcu[c]);
      qr *= f;
      qg *= f;
      qb *= f;
    }
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
  if (sem("init_solver") && s.to_init_solver) {
    if (m.IsRoot()) {
      std::cout << "backend=" << var.String["backend"] << std::endl;
    }
    {
      typename Proj<M>::Par p;
      p.conv = Conv::exp;
      std::shared_ptr<linear::Solver<M>> linsolver(
          ULinear<M>::MakeLinearSolver(var, "symm", m));
      FieldCell<Vect> fcvel(m, Vect(0));

      const auto topvel = var.Double["topvel"];
      for (auto f : m.AllFacesM()) {
        size_t nci;
        if (m.IsBoundary(f, nci)) {
          { // fluid
            auto& bc = bc_fluid[f];
            bc.nci = nci;
            bc.type = BCondFluidType::wall;
            if (f[1] == m.GetGlobalSize()[1]) {
              bc.velocity = Vect(topvel, 0);
            }
          }
          { // vof
            auto& bc = bc_vof[f];
            bc.nci = nci;
            bc.halo = BCondAdvection<Scal>::Halo::fill;
            bc.fill_vf = 0;
          }
        }
      }
      const ProjArgs<M> args{fcvel,     bc_fluid,   mc_velcond, &fc_rho, &fc_mu,
                             &fc_force, &ff_bforce, &fc_src,    &fc_src, 0,
                             s.dt,      linsolver,  p};
      fluid.reset(new Proj<M>(m, m, args));
    }
    {
      typename Vof<M>::Par p;
      p.dim = 2;
      p.sharpen = var.Int["sharpen"];
      p.sharpen_cfl = var.Double["sharpen_cfl"];
      FieldCell<Scal> fcu(m, 0);
      const FieldCell<Scal> fccl(m, 0);
      vof.reset(new Vof<M>(
          m, m, fcu, fccl, bc_vof, &fluid->GetVolumeFlux(), &fc_src, 0, s.dt,
          p));
    }

    auto ps = ParsePar<PartStr<Scal>>()(m.GetCellSize().norminf(), var);
    psm_par = ParsePar<PartStrMeshM<M>>()(ps, var);
  }
  if (sem.Nested("start")) {
    fluid->StartStep();
  }
  if (sem.Nested("iter")) {
    fluid->MakeIteration();
  }
  if (sem.Nested("finish")) {
    fluid->FinishStep();
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
  if (sem.Nested("post")) {
    vof->PostStep();
  }
  if (sem("maxvel")) {
    maxvel = 0;
    const auto& ffv = fluid->GetVolumeFlux();
    for (auto f : m.Faces()) {
      maxvel = std::max(maxvel, std::abs(ffv[f] / m.GetArea(f)));
    }
    m.Reduce(&maxvel, Reduction::max);
  }
  if (sem.Nested("curv")) {
    UCurv<M>::CalcCurvPart(vof.get(), psm_par, &fc_curv, m);
  }
  if (sem("flux")) {
    const Scal rho1 = var.Double["rho1"];
    const Scal rho2 = var.Double["rho2"];
    const Scal mu1 = var.Double["mu1"];
    const Scal mu2 = var.Double["mu2"];
    const Scal sigma = var.Double["sigma"];
    const Vect gravity(var.Vect["gravity"]);

    const Scal h = m.GetCellSize()[0];
    s.dt = var.Double["cfl"] * h / maxvel;
    if (mu1 != 0) {
      s.dt = std::min(s.dt, var.Double["cflvis"] * h * h / mu1);
    }
    if (sigma != 0) {
      s.dt = std::min(
          s.dt, var.Double["cflsurf"] * std::sqrt(
                                            std::pow(h, 3) * (rho1 + rho2) /
                                            (4 * M_PI * std::abs(sigma))));
    }

    s.dt = std::min(s.dt, var.Double["dtmax"]);
    fluid->SetTimeStep(s.dt);
    vof->SetTimeStep(s.dt);
    const auto& fcu = vof->GetField();
    for (auto c : m.AllCells()) {
      const Scal u = fcu[c];
      fc_rho[c] = (1 - u) * rho1 + u * rho2;
      fc_mu[c] = (1 - u) * mu1 + u * mu2;
    }
    for (auto f : m.SuFacesM()) {
      auto rho = (fc_rho[f.cm] + fc_rho[f.cp]) * 0.5;
      ff_bforce[f] = rho * gravity.dot(m.GetNormal(f));
    }

    const FieldFace<Scal> ff_sigma(m, var.Double["sigma"]);
    AppendSurfaceTension(m, ff_bforce, fcu, fc_curv, ff_sigma);
  }
  if (sem("init") && s.to_init_field) {
    vof->AddModifier([](FieldCell<Scal>& fcu, FieldCell<Scal>&, const M& m) { //
      fcu.Reinit(m, 0);
    });
  }
  if (sem("spawn") && s.to_spawn) {
    vof->AddModifier(
        [&](FieldCell<Scal>& fcu, FieldCell<Scal>&, const M& m) { //
          FieldCell<Scal> fc_add(m);
          GetCircle(fc_add, s.spawn_c, s.spawn_r, m);
          for (auto c : m.Cells()) {
            fcu[c] = std::max(fcu[c], fc_add[c]);
          }
        });
  }
  if (sem()) {
    s.to_init_solver = false;
    s.to_init_field = false;
    s.to_spawn = false;

    RenderField(
        *g_canvas, vof->GetField(), fluid->GetVelocity(), fluid->GetPressure(),
        m);

    // Render interafce lines
    if (m.IsRoot()) {
      s.lines.clear();
    }
    auto h = m.GetCellSize();
    const auto& fci = vof->GetMask();
    const auto& fcn = vof->GetNormal();
    const auto& fca = vof->GetAlpha();
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

    if (m.IsRoot()) {
      ++s.step;
      s.time += vof->GetTimeStep();
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

  const int nsteps = g_var.Int["nsteps"];
  for (int i = 0; i < nsteps; ++i) {
    s.distrsolver.Run();
  }

  const Scal tstep = timer.GetSeconds() / nsteps;

  if (s.step % 10 == 0) {
    std::cout << util::Format(
                     "step={:05} t={:.5f} dt={:.4f} g={:.1f} tstep={:.1f}ms",
                     s.step, s.time, s.dt, Vect(g_var.Vect["gravity"]),
                     tstep * 1e3)
              << std::endl;
  }
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
set double sharpen_cfl 0.1

set int nsteps 1

set double cfl 1
set double cflvis 0.5
set double cflsurf 2

set int part 1
set double part_h 4
set double part_relax 0.5
set int part_np 7
set double part_segcirc 0
set int part_dn 0
set int part_dump_fr 1
set int part_ns 3
set double part_tol 1e-5
set int part_itermax 20
set int part_verb 0
set int dim 2
set int vtkbin 1
set int vtkmerge 1
)EOF";
}

extern "C" {
Scal AddVelocityAngle(Scal add_deg) {
  auto& var_g = g_var.Vect["gravity"];
  Vect g(var_g);
  Scal deg = std::atan2(g[1], g[0]) * 180. / M_PI;
  deg += add_deg;
  const Scal rad = deg * M_PI / 180.;
  g = Vect(std::cos(rad), std::sin(rad)) * g.norm();
  std::cout << util::Format("g={:.3f} angle={:.1f}", g, deg) << std::endl;
  var_g = g;
  return deg;
}
void Init() {
  auto& s = *g_state;
  s.to_init_field = true;
}
void Spawn(float x, float y, float r) {
  auto& s = *g_state;
  s.to_spawn = true;
  s.spawn_c = Vect(x, y);
  s.spawn_r = r;
  std::cout << util::Format("spawn c={:.3f} r={:.3f}", s.spawn_c, s.spawn_r)
            << std::endl;
}
int TogglePause() {
  auto& s = *g_state;
  s.pause = !s.pause;
  return s.pause;
}
void SetExtraConfig(const char* extra) {
  g_extra_config = extra;
}

void SetMesh(int nx) {
  MPI_Comm comm = 0;
  std::stringstream conf;
  conf << GetDefaultConf();
  Subdomains<MIdx> sub(MIdx(nx), MIdx(nx), 1);
  conf << sub.GetConfig();
  conf << GetBaseConfig();
  conf << g_extra_config << '\n';
  Parser(g_var).ParseStream(conf);

  std::shared_ptr<State> new_state;
  new_state = std::make_shared<State>(comm, g_var);

  g_state = new_state;
  std::cout << util::Format("mesh {}", MIdx(nx)) << std::endl;
}
void SetCanvas(int nx, int ny) {
  g_canvas = std::make_shared<Canvas>(MIdx(nx, ny));
  std::cout << util::Format("canvas {}", g_canvas->size) << std::endl;
}
int GetLines(uint16_t* data, int max_size) {
  auto state = g_state;
  auto& s = *state;
  int i = 0;
  for (auto p : s.lines) {
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

  SetCanvas(500, 500);
  emscripten_set_canvas_element_size(
      "#canvas", g_canvas->size[0], g_canvas->size[1]);
  emscripten_set_main_loop(main_loop, 30, 1);
}
