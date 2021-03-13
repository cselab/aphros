// Created by Petr Karnakov on 24.12.2020
// Copyright 2020 ETH Zurich

#include <cassert>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>

#include <distr/distrbasic.h>
#include <dump/dump.h>
#include <dump/raw.h>
#include <geom/mesh.h>
#include <parse/argparse.h>
#include <parse/parser.h>
#include <parse/util.h>
#include <solver/approx_eb.h>
#include <solver/cond.h>
#include <solver/embed.h>
#include <util/filesystem.h>
#include <util/format.h>

using M = MeshCartesian<double, 4>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using MIdx = typename M::MIdx;

void DumpCsv(
    std::string path, FieldCell<Scal>& fcu, M& m, bool verbose = false) {
  auto sem = m.GetSem();
  struct {
    std::vector<std::string> csv_names;
    std::vector<std::vector<Scal>> csv_data;
  } * ctx(sem);
  auto& t = *ctx;
  if (sem()) {
    t.csv_names.clear();
    t.csv_data.clear();
    for (auto d : m.dirs) {
      std::vector<Scal> data;
      for (auto c : m.CellsM()) {
        data.push_back(c.center[d]);
      }
      const std::string name = std::string() + M::Dir(d).letter();
      t.csv_names.push_back(name);
      t.csv_data.push_back(data);
    }

    {
      std::vector<Scal> data;
      for (auto c : m.Cells()) {
        data.push_back(fcu[c]);
      }
      const std::string name = "u";
      t.csv_names.push_back(name);
      t.csv_data.push_back(data);
    }
    for (auto& v : t.csv_data) {
      m.Reduce(&v, Reduction::concat);
    }
  }
  if (sem()) {
    if (m.IsRoot()) {
      std::vector<std::pair<std::string, std::vector<Scal>>> csv;
      for (size_t i = 0; i < t.csv_names.size(); ++i) {
        csv.emplace_back(t.csv_names[i], t.csv_data[i]);
      }
      if (verbose) {
        std::cout << path << std::endl;
      }
      DumpCsv<Scal>(csv, path);
    }
  }
}

void DumpRaw(
    std::string path, FieldCell<Scal>& fcu, M& m, bool verbose = false) {
  auto sem = m.GetSem();
  using Raw = dump::Raw<M>;
  using Xmf = dump::Xmf<Vect>;
  struct {
    Raw::Meta meta;
  } * ctx(sem);
  auto& t = *ctx;
  if (sem("writexmf")) {
    t.meta = Xmf::GetMeta(MIdx(0), MIdx(1), m);
    t.meta.name = "u";
    t.meta.binpath = path;
    if (m.IsRoot()) {
      const auto xmfpath = util::SplitExt(path)[0] + ".xmf";
      Xmf::WriteXmf(xmfpath, t.meta);
      if (verbose) {
        std::cout << path << ' ' << xmfpath << std::endl;
      }
    }
  }
  if (sem.Nested("write")) {
    Raw::Write(fcu, t.meta, path, m);
  }
}

void Run(M& m, Vars& var) {
  using Scal = typename M::Scal;
  auto sem = m.GetSem();
  struct {
    FieldCell<Scal> fcu; // field
    MapEmbed<BCond<Scal>> mebc; // face conditions
    int step = 0;
    Scal sum0;
    Scal sum2;
    Scal max;
    std::set<std::string> dumpformats;
  } * ctx(sem);
  auto& t = *ctx;
  if (sem()) {
    t.dumpformats = GetWords(var.String["dumpformats"]);
    if (m.IsRoot()) {
      std::ofstream out("out.conf");
      Parser::PrintVars(var, out);
    }
  }
  if (sem("initial")) {
    t.fcu.Reinit(m, 0);
    for (auto c : m.AllCellsM()) {
      const Vect xc = m.GetGlobalLength() * 0.5;
      t.fcu[c] = xc.sqrdist(c.center) < 0.2 ? 1 : 0;
    }
  }
  sem.LoopBegin();
  if (sem() && m.IsRoot()) {
    if (t.step % var.Int["reportevery"] == 0) {
      std::cout << util::Format(
          "{}: sum0={:.3f} sum2={:.3f} max={:.3f}\n", t.step, t.sum0, t.sum2,
          t.max);
    }
  }
  if (sem.Nested() && t.step % var.Int["dumpevery"] == 0 &&
      t.dumpformats.count("csv")) {
    DumpCsv(util::Format("u_{:04d}.csv", t.step), t.fcu, m, true);
  }
  if (sem.Nested() && t.step % var.Int["dumpevery"] == 0 &&
      t.dumpformats.count("raw")) {
    DumpRaw(util::Format("u_{:04d}.raw", t.step), t.fcu, m, true);
  }
  if (sem("step")) {
    const Scal dt = var.Double["dt"];
    const auto ffg = UEmbed<M>::Gradient(t.fcu, t.mebc, m);
    FieldFace<Scal> ff_flux(m);
    for (auto f : m.FacesM()) {
      ff_flux[f] = ffg[f] * f.area;
    }

    for (auto c : m.CellsM()) {
      Scal sum = 0;
      for (auto q : m.Nci(c)) {
        sum += ff_flux[c.face(q)] * c.outward_factor(q);
      }
      t.fcu[c] += sum * dt / c.volume;
    }
    m.Comm(&t.fcu);

    t.sum0 = 0;
    t.sum2 = 0;
    t.max = 0;
    for (auto c : m.CellsM()) {
      const Vect xc = m.GetGlobalLength() * 0.5;
      t.sum0 += t.fcu[c] * c.volume;
      t.sum2 += t.fcu[c] * xc.sqrdist(c.center) * c.volume;
      t.max = std::max(t.max, t.fcu[c]);
    }
    m.Reduce(&t.sum0, Reduction::sum);
    m.Reduce(&t.sum2, Reduction::sum);
    m.Reduce(&t.max, Reduction::max);
    ++t.step;
  }
  if (sem()) {
    if (t.step >= var.Int["nsteps"]) {
      sem.LoopBreak();
    }
  }
  sem.LoopEnd();
}

int main(int argc, const char** argv) {
  MpiWrapper mpi(&argc, &argv);

  ArgumentParser argparser("Diffusion solver in 4D", mpi.IsRoot());
  argparser.AddVariable<std::string>("config", "a.conf")
      .Help("Path to configuration file");
  argparser.AddVariable<std::string>("--extra", "")
      .Help("Extra configuration (commands 'set ... ')");

  auto args = argparser.ParseArgs(argc, argv);
  if (const int* p = args.Int.Find("EXIT")) {
    return *p;
  }

  Vars var;
  Parser parser(var);
  const std::string confpath = args.String["config"];
  std::ifstream fconf(confpath);
  fassert(fconf.good(), "Can't open file '" + confpath + "'");
  std::stringstream conf;
  conf << fconf.rdbuf();
  {
    std::stringstream s(conf.str());
    parser.ParseStream(s);
  }

  const MIdx mesh_size(
      var.Int["nx"], var.Int["ny"], var.Int["nz"], var.Int["nw"]);
  const MIdx block_size(
      var.Int["bsx"], var.Int["bsy"], var.Int["bsz"], var.Int["bsw"]);

  Subdomains<MIdx> sub(mesh_size, block_size, mpi.GetCommSize());
  conf << sub.GetConfig();

  conf << args.String["extra"] << '\n';

  return RunMpiBasicString<M>(mpi, Run, conf.str());
}
