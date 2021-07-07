// Created by Petr Karnakov on 28.12.2019
// Copyright 2019 ETH Zurich

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>

#include <debug/isnan.h>
#include <distr/distrbasic.h>
#include <dump/dumper.h>
#include <dump/vtk.h>
#include <func/init.h>
#include <func/init_u.h>
#include <linear/linear.h>
#include <parse/argparse.h>
#include <solver/approx_eb.h>
#include <solver/curv.h>
#include <solver/reconst.h>
#include <solver/vofm.h>
#include <util/hydro.h>
#include <util/linear.h>
#include "func/init_bc.h"

using M = MeshCartesian<double, 2>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using MIdx = typename M::MIdx;
using R = Reconst<Scal>;
constexpr Scal kClNone = Vofm<M>::kClNone;

// Returns unique colors found in cells `cc`
template <class M>
std::vector<Scal> GetUniqueColors(
    const std::vector<IdxCell>& cc, const GRange<size_t>& layers,
    const Multi<const FieldCell<Scal>*> fccl, const M&) {
  std::vector<Scal> colors;
  for (auto c : cc) {
    for (auto l : layers) {
      const Scal cl = (*fccl[l])[c];
      if (cl != kClNone &&
          std::find(colors.begin(), colors.end(), cl) == colors.end()) {
        colors.push_back(cl);
      }
    }
  }
  return colors;
}

// Returns the number of unique colors found in cells `cc`
template <class M>
size_t GetNumColors(
    const std::vector<IdxCell>& cc, const GRange<size_t>& layers,
    const Multi<const FieldCell<Scal>*> fccl, const M& m) {
  return GetUniqueColors(cc, layers, fccl, m).size();
}

// Modifies volume flux `ff_flux` to make it divergence-free
// in cells with non-zero `fc_mask` with source term `fc_src`
// Solves equation:
//   ∇ · u = `src`
//   u = u* - ∇p
// where `u` is velocity after correction, `u*` is velocity before correction,
// volume flux is `flux = u* · Sf` with surface element `Sf`
template <class M>
void ProjectVolumeFluxMasked(
    FieldFace<typename M::Scal>& ff_flux, const FieldCell<Scal>& fc_mask,
    const FieldCell<Scal>& fc_src, std::shared_ptr<linear::Solver<M>> linsolver,
    M& m) {
  using Scal = typename M::Scal;
  using ExprFace = generic::Vect<Scal, 3>;
  using Expr = generic::Vect<Scal, M::dim * 2 + 2>;

  auto sem = m.GetSem();
  struct {
    FieldFace<ExprFace> ff_corr; // expression for corrected volume flux [i]
    FieldCell<Expr> fc_system; // linear system for pressure [i]
    FieldCell<Scal> fcp; // pressure (up to a constant)
  } * ctx(sem);
  auto& t = *ctx;

  if (sem("init")) {
    t.ff_corr.Reinit(m);
    for (auto f : m.FacesM()) {
      auto& e = t.ff_corr[f];
      size_t nci;
      const Scal h = m.GetCellSize()[0];
      if (m.IsBoundary(f, nci)) {
        e[0] = 0;
        e[1] = 0;
      } else {
        const Scal a = -f.area / h;
        e[0] = -a;
        e[1] = a;
      }
      e[2] = ff_flux[f];
    }

    t.fc_system.Reinit(m, Expr(0));
    for (auto c : m.CellsM()) {
      auto& e = t.fc_system[c];
      if (fc_mask[c]) {
        for (auto q : m.Nci(c)) {
          m.AppendExpr(e, t.ff_corr[c.face(q)] * c.outward_factor(q), q);
        }
        e /= c.volume;
        e.back() -= fc_src[c];
      } else {
        e[0] = 1;
        e.back() = 0;
      }
    }
  }
  if (sem.Nested("solve")) {
    linsolver->Solve(t.fc_system, nullptr, t.fcp, m);
  }
  if (sem("apply")) {
    for (auto f : m.FacesM()) {
      const auto& e = t.ff_corr[f];
      ff_flux[f] = e[0] * t.fcp[f.cm] + e[1] * t.fcp[f.cp] + e[2];
    }
  }
}

// Computes flux proportional to curvature.
//   fcu: volume fraction
//   fccl: color
//   fcn: normal
//   divfree: output divergence-free velocity
//   divfree_masked: apply divergence-free constraint only near interfaces
//   voidpenal: factor for penalization of voids
// Output:
//   ff_flux: result
//   fc_mask: mask for divergence-free constraint
template <class M>
void CalcMeanCurvatureFlowFlux(
    const GRange<size_t>& layers, const Multi<const FieldCell<Scal>*> fcu,
    const Multi<const FieldCell<Scal>*> fccl,
    const Multi<const FieldCell<Vect>*> fcn,
    const Multi<const FieldCell<Scal>*> fck,
    const MapEmbed<BCondFluid<Vect>>& mebc, Scal gamma, bool divfree,
    bool divfree_masked, Scal voidpenal, bool anchored_boundaries,
    std::shared_ptr<linear::Solver<M>> linsolver, M& m,
    /*out*/
    FieldFace<Scal>& ff_flux, FieldCell<Scal>* fc_mask) {
  (void)mebc;
  auto sem = m.GetSem();
  struct {
    FieldCell<Scal> fc_num_interfaces; // number of interfaces in cell
    FieldCell<Scal> fc_sum; // sum of volume fractions from all layers
    FieldCell<Scal> fc_mask; // non-zero in cells in which to solve Poisson
                             // 2: void cell
                             // 1: interface cell
                             // 0: otherwise
    FieldCell<Scal> fc_src; // source term for flux projection
  } * ctx(sem);
  auto& t = *ctx;

  if (sem("calc")) {
    t.fc_num_interfaces.Reinit(m, 0);
    for (auto c : m.Cells()) {
      for (auto l : layers) {
        if ((*fccl[l])[c] != kClNone) {
          const Scal k = (*fck[l])[c];
          const Vect n = (*fcn[l])[c];
          if (!IsNan(k) && !IsNan(n) && n.sqrnorm() > 0) {
            t.fc_num_interfaces[c] += 1;
          }
        }
      }
    }

    t.fc_mask.Reinit(m, 1);
    if (divfree_masked) {
      // Non-zero mask in cells with voids or interfaces
      for (auto c : m.Cells()) {
        t.fc_mask[c] = GetNumColors({c}, layers, fccl, m) == 0 ? 2
                       : t.fc_num_interfaces[c] > 0            ? 1
                                                               : 0;
      }
    }
    if (fc_mask) {
      (*fc_mask) = t.fc_mask;
    }

    t.fc_sum.Reinit(m, 0);
    for (auto c : m.AllCells()) {
      for (auto l : layers) {
        if ((*fccl[l])[c] != kClNone) {
          t.fc_sum[c] += (*fcu[l])[c];
        }
      }
    }

    t.fc_src.Reinit(m, 0);
    // Source term penalizing voids and overlaps
    for (auto c : m.CellsM()) {
      if (GetNumColors({c}, layers, fccl, m) <= 2) {
        const Scal h = m.GetCellSize()[0];
        t.fc_src[c] = -std::abs(1 - t.fc_sum[c]) * voidpenal * gamma / sqr(h);
      }
    }

    ff_flux.Reinit(m, 0);
    for (auto f : m.Faces()) {
      const IdxCell cm = m.GetCell(f, 0);
      const IdxCell cp = m.GetCell(f, 1);
      const std::vector<Scal> colors =
          GetUniqueColors({cm, cp}, layers, fccl, m);
      for (auto cl : colors) {
        Scal um = 0;
        Scal up = 0;
        Scal km = GetNan<Scal>();
        Scal kp = GetNan<Scal>();
        Vect nm(0);
        Vect np(0);
        for (auto l : layers) {
          if ((*fccl[l])[cm] == cl) {
            um = (*fcu[l])[cm];
            km = (*fck[l])[cm];
            nm = (*fcn[l])[cm];
          }
          if ((*fccl[l])[cp] == cl) {
            up = (*fcu[l])[cp];
            kp = (*fck[l])[cp];
            np = (*fcn[l])[cp];
          }
        }
        if (!IsNan(km) || !IsNan(kp)) {
          const Scal k = (std::abs(um - 0.5) < std::abs(up - 0.5) ? km : kp);
          const Vect n = (std::abs(um - 0.5) < std::abs(up - 0.5) ? nm : np);
          if (!IsNan(k) && !IsNan(n) && n.norm() > 0) {
            // Factor 0.5 since the force is computed for both components
            ff_flux[f] -= n.dot(m.GetSurface(f)) * (k * gamma / n.norm()) * 0.5;
          }
        }
      }
    }
    if (anchored_boundaries) {
      for (auto f : m.Faces()) {
        size_t nci;
        if (m.IsBoundary(f, nci)) {
          const IdxCell c = m.GetCell(f, nci);
          for (auto cn : m.Stencil5(c)) {
            for (auto q : m.Nci(cn)) {
              ff_flux[m.GetFace(cn, q)] = 0;
            }
          }
        }
      }
    }
  }
  if (divfree && sem.Nested("proj")) {
    ProjectVolumeFluxMasked(ff_flux, t.fc_mask, t.fc_src, linsolver, m);
  }
}

struct TrajEntry {
  Scal volume = 0;
  Vect center{Vect(0)};
  Scal k = 0; // average curvature
  Scal cells_k = 0; // cells with defined curvature
};

/** Returns trajectory entries for components of selected colors.
colors: colors to select

Returns:
array of trajectory entries for each color in `colors`
*/
void CalcTraj(
    const Vofm<M>::Plic& plic, const Multi<const FieldCell<Scal>*>& fck,
    const std::vector<Scal>& colors, M& m,
    /*out*/
    std::vector<TrajEntry>& entries) {
  auto sem = m.GetSem();
  struct {
    std::map<Scal, TrajEntry> map;
  } * ctx(sem);
  auto& t = *ctx;
  if (sem()) {
    // Create empty entries for all colors
    for (auto l : plic.layers) {
      for (auto c : m.Cells()) {
        const auto cl = (*plic.vfccl[l])[c];
        if (cl != kClNone) {
          t.map[cl];
        }
      }
    }
    for (auto l : plic.layers) {
      for (auto c : m.CellsM()) {
        const auto cl = (*plic.vfccl[l])[c];
        if (cl != kClNone) {
          auto& entry = t.map.at(cl);
          entry.volume += (*plic.vfcu[l])[c] * c.volume();
          entry.center += (*plic.vfcu[l])[c] * c.volume() * c.center();
          if (!IsNan((*fck[l])[c])) {
            entry.k += (*fck[l])[c];
            entry.cells_k += 1;
          }
        }
      }
    }
    // Copy selected entries from map to array
    entries.resize(colors.size());
    for (size_t i = 0; i < colors.size(); ++i) {
      auto it = t.map.find(colors[i]);
      if (it != t.map.end()) {
        entries[i] = it->second;
      }
      m.Reduce(&entries[i].volume, Reduction::sum);
      for (auto d : M::dirs) {
        m.Reduce(&entries[i].center[d], Reduction::sum);
      }
      m.Reduce(&entries[i].k, Reduction::sum);
      m.Reduce(&entries[i].cells_k, Reduction::sum);
    }
  }
  if (sem()) {
    for (auto& entry : entries) {
      if (entry.volume > 0) {
        entry.center /= entry.volume;
        entry.k /= entry.cells_k;
      }
    }
  }
}

void WriteTrajHeader(std::ostream& out, size_t num_colors) {
  bool first = true;
  auto delim = [&]() {
    if (first) {
      first = false;
    } else {
      out << ' ';
    }
  };
  delim();
  out << "t";
  delim();
  out << "frame";
  for (auto field : {"x", "y", "volume", "k"}) {
    for (size_t i = 0; i < num_colors; ++i) {
      delim();
      out << util::Format("{:}_{:}", field, i);
    }
  }
  out << '\n';
}

void WriteTrajEntry(
    std::ostream& out, Scal time, size_t frame,
    const std::vector<TrajEntry>& entries) {
  bool first = true;
  auto delim = [&]() {
    if (first) {
      first = false;
    } else {
      out << ' ';
    }
  };
  delim();
  out << time;
  delim();
  out << frame;
  for (auto& entry : entries) {
    delim();
    out << entry.center[0];
  }
  for (auto& entry : entries) {
    delim();
    out << entry.center[1];
  }
  for (auto& entry : entries) {
    delim();
    out << entry.volume;
  }
  for (auto& entry : entries) {
    delim();
    out << entry.k;
  }
  out << '\n';
}

Scal GetTimeStep(M& m, Vars& var) {
  const Scal cflst = var.Double["cflst"];
  const Scal gamma = var.Double["gamma"];
  return cflst * sqr(m.GetCellSize()[0]) * std::abs(gamma);
}

template <class M>
void SetZeroBoundaryFlux(FieldFace<Scal>& fcu, const M& m) {
  for (auto fb : m.Faces()) {
    size_t nci;
    if (m.IsBoundary(fb, nci)) {
      const IdxCell c = m.GetCell(fb, nci);
      for (auto q : m.Nci(c)) {
        const IdxFace f = m.GetFace(c, q);
        fcu[f] = 0;
      }
    }
  }
}

void Run(M& m, Vars& var) {
  auto sem = m.GetSem();
  struct {
    std::unique_ptr<Vofm<M>> as; // advection solver
    Multi<FieldCell<Scal>> fcu0; // initial volume fraction
    Multi<FieldCell<Scal>> fccl0; // initial color
    FieldCell<Scal> fc_src; // volume source
    FieldEmbed<Scal> fe_flux; // volume flux
    FieldCell<Vect> fc_vel; // velocity
    FieldCell<Scal> fc_mask; // mask for divergence-free constraint
    Multi<FieldCell<Scal>> fck; // curvature
    FieldCell<Scal> fck_single; // curvature from single layer
    MapEmbed<BCondFluid<Vect>> mebc_fluid; // face conditions
    MapEmbed<BCondAdvection<Scal>> mebc_adv; // face conditions
    GRange<size_t> layers;
    std::unique_ptr<curvature::Estimator<M>> curv;
    Scal dt = 0;
    size_t step = 0;
    size_t frame = 0;
    std::unique_ptr<Dumper> dumper;
    std::shared_ptr<linear::Solver<M>> linsolver;
    // `traj[frame][cl]` is trajectory entry of component with color `cl`
    std::vector<std::vector<TrajEntry>> traj;
    std::ofstream trajfile;
    bool dumptraj;
    size_t dumptraj_colors;
  } * ctx(sem);
  auto& layers = ctx->layers;
  auto& t = *ctx;

  auto& as = ctx->as;
  const Scal tmax = var.Double["tmax"];
  const Scal dtmax = var.Double["dtmax"];
  if (sem("init")) {
    t.dumper = std::make_unique<Dumper>(var, "dump_");
    t.fc_src.Reinit(m, 0);
    t.fe_flux.Reinit(m, 0);
    t.dumptraj = var.Int["dumptraj"];
    t.dumptraj_colors = var.Int["dumptraj_colors"];

    t.linsolver = ULinear<M>::MakeLinearSolver(var, "symm", m);

    m.flags.is_periodic[0] = var.Int["hypre_periodic_x"];
    m.flags.is_periodic[1] = var.Int["hypre_periodic_y"];
    m.flags.linreport = var.Int["linreport"];
    m.flags.edim = 2;

    {
      auto p = InitBc(var, m, {}, {});
      ctx->mebc_fluid = std::get<0>(p);
      ctx->mebc_adv = std::get<1>(p);
    }

    typename Vofm<M>::Par parvof;
    parvof.sharpen = var.Int["sharpen"];
    parvof.sharpen_cfl = var.Double["sharpen_cfl"];
    parvof.layers = var.Int["vofm_layers"];
    parvof.clipth = var.Double["clipth"];
    parvof.filterth = var.Double["filterth"];
    parvof.dim = 2;
    parvof.vtkbin = true;
    parvof.vtkpoly = false;
    layers = GRange<size_t>(parvof.layers);
    t.fcu0.resize(layers);
    t.fccl0.resize(layers);

    auto buf = ReadPrimList(var.String["list_path"], m.IsRoot());
    InitOverlappingComponents(buf, t.fcu0, t.fccl0, layers, m);
    as.reset(new Vofm<M>(
        m, m, t.fcu0, t.fccl0, ctx->mebc_adv, &t.fe_flux, &t.fc_src, 0, t.dt,
        parvof));
    auto mod = [&t](auto& fcu, auto& fccl, auto layers, auto&) {
      for (auto l : layers) {
        (*fcu[l]) = t.fcu0[l];
        (*fccl[l]) = t.fccl0[l];
      }
    };
    as->AddModifier(mod);

    if (t.dumptraj && m.IsRoot()) {
      t.trajfile.open("traj.dat");
      WriteTrajHeader(t.trajfile, t.dumptraj_colors);
    }

    t.fck.Reinit(layers, m, 0);
    t.curv = curvature::MakeEstimator(var, m, layers);

    if (m.IsRoot()) {
      std::cout << util::Format("meancurvflow dt={:}\n", GetTimeStep(m, var));
    }
  }
  sem.LoopBegin();
  if (sem.Nested("start")) {
    as->StartStep();
  }
  if (sem.Nested("iter")) {
    as->MakeIteration();
  }
  if (sem.Nested("finish")) {
    as->FinishStep();
  }
  if (sem.Nested("post")) {
    as->PostStep();
  }
  if (sem.Nested()) {
    t.curv->CalcCurvature(t.fck, as->GetPlic(), m, m);
  }
  if (sem.Nested("flux")) {
    CalcMeanCurvatureFlowFlux(
        layers, as->GetFieldM(), as->GetColor(), as->GetNormal(), t.fck,
        ctx->mebc_fluid, var.Double["gamma"], var.Int["divfree"],
        var.Int["divfree_masked"], var.Double["voidpenal"],
        var.Int["anchored_boundaries"], t.linsolver, m,
        t.fe_flux.GetFieldFace(), &t.fc_mask);
  }
  if (sem("dt")) {
    auto mod = [&var, &t](auto& fcu, auto& fccl, auto layers, auto& m) {
      if (var.Int["anchored_boundaries"]) {
        for (auto f : m.Faces()) {
          size_t nci;
          if (m.IsBoundary(f, nci)) {
            const IdxCell c = m.GetCell(f, nci);
            for (auto cn : m.Stencil5(c)) {
              for (auto l : layers) {
                (*fcu[l])[cn] = t.fcu0[l][cn];
                (*fccl[l])[cn] = t.fccl0[l][cn];
              }
            }
          }
        }
      }
    };
    as->AddModifier(mod);
    Scal maxv = 0;
    for (auto f : m.Faces()) {
      maxv = std::max(maxv, std::abs(t.fe_flux[f]));
    }
    const Scal cfl = var.Double["cfl"];
    t.dt = cfl * m.GetCellSize().prod() / maxv;
    t.dt = std::min(t.dt, GetTimeStep(m, var));
    m.Reduce(&t.dt, "min");
  }
  if (sem("mindt")) {
    t.dt = std::min(t.dt, dtmax);
    as->SetTimeStep(t.dt);
    if (m.IsRoot() && t.step % var.Int["report_every"] == 0) {
      std::cout << util::Format(
          "STEP={:04d} t={:.6f} dt={:}\n", t.step, as->GetTime(), t.dt);
    }
  }
  const bool dump = t.dumper->Try(as->GetTime(), as->GetTimeStep());
  if (sem.Nested() && var.Int["dumppoly"] && dump) {
    as->DumpInterface(GetDumpName("s", ".vtk", t.frame), {t.fck}, {"k"});
  }
  if (sem.Nested() && var.Int["dumppolymarch"] && dump) {
    as->DumpInterfaceMarch(GetDumpName("sm", ".vtk", t.frame));
  }
  if (sem.Nested() && var.Int["dumpaux"] && dump) {
    t.curv->DumpAux("heights", t.frame, m);
  }
  if (t.dumptraj) {
    if (sem() && dump) {
      t.traj.emplace_back();
    }
    if (sem.Nested() && dump) {
      const auto plic = as->GetPlic();
      std::vector<Scal> colors(t.dumptraj_colors);
      std::iota(colors.begin(), colors.end(), 0);
      CalcTraj(plic, t.fck, colors, m, t.traj.back());
    }
    if (sem() && dump) {
      if (m.IsRoot()) {
        WriteTrajEntry(t.trajfile, as->GetTime(), t.frame, t.traj.back());
        t.trajfile.flush();
      }
    }
  }
  if (sem("checkloop")) {
    if (dump) {
      ++t.frame;
      if (var.Int["dumpfields"]) {
        // Compute velocity from flux
        t.fc_vel.Reinit(m, Vect(0));
        for (auto c : m.AllCells()) {
          for (auto q : m.Nci(c)) {
            const auto f = m.GetFace(c, q);
            t.fc_vel[c] += m.GetNormal(f) * (t.fe_flux[f] * 0.5 / m.GetArea(f));
          }
        }
        // Fill curvature from any layer
        t.fck_single.Reinit(m, GetNan<Scal>());
        const auto plic = as->GetPlic();
        for (auto c : m.Cells()) {
          for (auto l : layers) {
            if (!IsNan(t.fck[l][c])) {
              t.fck_single[c] = t.fck[l][c];
              break;
            }
          }
        }
        m.Dump(&t.fc_vel, 0, "vx");
        m.Dump(&t.fc_vel, 1, "vy");
        m.Dump(&t.fck_single, "k");
        m.Dump(&t.fc_mask, "mask");
        for (auto l : layers) {
          m.Dump(as->GetPlic().vfcu[l], "u" + std::to_string(l));
          m.Dump(as->GetPlic().vfccl[l], "cl" + std::to_string(l));
        }
      }
    }
    if (as->GetTime() >= tmax) {
      sem.LoopBreak();
    }
    ++t.step;
  }
  sem.LoopEnd();
}

int main(int argc, const char** argv) {
  ArgumentParser parser("Constrained mean curvature flow");
  parser.AddVariable<std::string>("config", "std.conf")
      .Help("Path to config file");
  parser.AddVariable<int>("--nx", 16).Help("Mesh size in x-direction");
  parser.AddVariable<std::string>("--extra", "")
      .Help("Extra configuration (commands 'set ... ')");
  parser.AddVariable<int>("--ny", 0).Help(
      "Mesh size in y-direction. Defaults to NX");
  parser.AddVariable<int>("--bs", 0).Help(
      "Block size in all directions. Defaults to min(64, NX)");
  const auto args = parser.ParseArgs(argc, argv);
  if (const int* p = args.Int.Find("EXIT")) {
    return *p;
  }

  const auto confpath = args.String["config"];
  std::ifstream fconf(confpath);

  if (!fconf.good()) {
    throw std::runtime_error("Can't open config '" + confpath + "'");
  }

  std::stringstream conf;
  conf << fconf.rdbuf();

  const int nx = args.Int["nx"];
  const int ny = args.Int["ny"] ? args.Int["ny"] : args.Int["nx"];
  MIdx mesh_size(1);
  mesh_size[0] = nx;
  mesh_size[1] = ny;
  MIdx block_size(1);
  block_size[0] = args.Int["bs"] ? args.Int["bs"] : std::min(64, nx);
  block_size[1] = args.Int["bs"] ? args.Int["bs"] : std::min(64, ny);

  MpiWrapper mpi(&argc, &argv);
  Subdomains<MIdx> sub(mesh_size, block_size, mpi.GetCommSize());
  conf << sub.GetConfig() << '\n';
  conf << args.String["extra"] << '\n';

  return RunMpiBasicString<M>(mpi, Run, conf.str());
}
