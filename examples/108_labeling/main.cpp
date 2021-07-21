// Created by Petr Karnakov on 25.09.2020
// Copyright 2020 ETH Zurich

#include <cassert>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

#include <distr/distrbasic.h>
#include <dump/hdf.h>
#include <func/init.h>
#include <parse/argparse.h>
#include <parse/vars.h>
#include <util/distr.h>
#include <util/filesystem.h>
#include <util/format.h>
#include <util/vof.h>

using M = MeshCartesian<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using MIdx = typename M::MIdx;

template <class M>
class LabelingCustom : public Labeling<M> {
 public:
  using Base = Labeling<M>;
  using Conf = typename Base::Conf;
  using Scal = typename M::Scal;
  using Base::conf;
  LabelingCustom(const Conf& conf_, const M&) : Base(conf_) {}
  ~LabelingCustom() {}
  void Recolor(
      const GRange<size_t>& layers, const Multi<const FieldCell<Scal>*>& fcu,
      const Multi<FieldCell<Scal>*>& fccl,
      const Multi<const FieldCell<Scal>*>& fccl_stable,
      const MapEmbed<BCond<Scal>>& mebc_cl, M& m) override {
    (void)layers;
    (void)fcu;
    (void)fccl;
    (void)fccl_stable;
    (void)mebc_cl;
    if (m.IsRoot()) {
      std::cerr << util::Format(
          "Custom implementation of Labeling.\n"
          "Warning: not implemented.\n");
    }
  }
};

template <class M>
class ModuleLabelingCustom : public ModuleLabeling<M> {
 public:
  using Conf = typename Labeling<M>::Conf;
  ModuleLabelingCustom() : ModuleLabeling<M>("custom") {}
  std::unique_ptr<Labeling<M>> Make(const Conf& conf, const M& m) override {
    return std::make_unique<LabelingCustom<M>>(conf, m);
  }
};

bool kReg_custom[] = {RegisterModule<ModuleLabelingCustom<M>>()};

void Init(
    GRange<size_t> layers, Multi<FieldCell<Scal>>& fcu, std::string prefix,
    const Vars& var, M& m) {
  auto sem = m.GetSem();
  struct {
    Vars var;
  } * ctx(sem);
  auto& t = *ctx;
  if (sem()) {
    fcu.Reinit(layers, m, 0);
  }
  for (auto l : layers) {
    const auto sl = std::to_string(l);
    if (sem()) {
      t.var.String.Set("init_vf", "list");
      t.var.String.Set("list_path", var.String["init_" + prefix + sl]);
      t.var.Int.Set("dim", m.GetGlobalSize()[2] == 1 ? 2 : 3);
      t.var.Int.Set("list_ls", 1);
    }
    if (sem.Nested()) {
      InitVf(fcu[l], t.var, m, false);
    }
  }
}

void Dump(
    const FieldCell<Scal>& fcu, std::string fieldname, std::string filename,
    M& m) {
  auto sem = m.GetSem();
  const auto path = filename + ".h5";
  if (sem.Nested()) {
    Hdf<M>::Write(fcu, path, m);
  }
  if (sem() && m.IsRoot()) {
    Hdf<M>::WriteXmf(util::SplitExt(path)[0] + ".xmf", fieldname, path, m);
  }
}

void Dump(
    GRange<size_t> layers, const Multi<FieldCell<Scal>>& fcu,
    std::string fieldname, std::string filename, M& m) {
  auto sem = m.GetSem();
  for (auto l : layers) {
    const auto sl = std::to_string(l);
    if (sem.Nested()) {
      Dump(fcu[l], fieldname, filename + sl, m);
    }
  }
}

void Run(M& m, Vars& var) {
  auto sem = m.GetSem();
  struct {
    GRange<size_t> layers;
    Multi<FieldCell<Scal>> fcvf; // volume fraction
    Multi<FieldCell<Scal>> fccl; // color
    MapEmbed<BCond<Scal>> mebc; // boundary conditions for volume fraction
    FieldCell<Scal> fcvf_sum; // sum of volume fractions from all layers
    FieldCell<Scal> fccl_sum; // colors from all layers
    std::unique_ptr<Labeling<M>> labeling; // labeling algorithm
  } * ctx(sem);
  constexpr Scal kClNone = -1;
  auto& t = *ctx;
  if (sem()) {
    t.layers = GRange<size_t>(var.Int["layers"]);
    t.fccl.Reinit(t.layers, m, kClNone);

    const auto modname = var.String["labeling"];
    if (auto* mod = ModuleLabeling<M>::GetInstance(modname)) {
      Labeling<M>::Conf conf;
      conf.verbose = true;
      t.labeling = mod->Make(conf, m);
    } else {
      fassert(
          false,
          util::Format(
              "Connected component labeling module '{}' not found", modname));
    }
  }
  if (sem.Nested()) {
    Init(t.layers, t.fcvf, "vf", var, m);
  }
  if (sem()) {
    // Clear color in cells with zero volume fraction
    for (auto l : t.layers) {
      for (auto c : m.AllCells()) {
        t.fccl[l][c] = (t.fcvf[l][c] > 0 ? l : kClNone);
      }
    }
  }
  if (sem.Nested()) {
    t.labeling->Recolor(t.layers, t.fcvf, t.fccl, t.fccl, t.mebc, m);
  }
  if (sem.Nested()) {
    Dump(t.layers, t.fcvf, "vf", "vf", m);
  }
  if (sem.Nested()) {
    Dump(t.layers, t.fccl, "cl", "cl", m);
  }
  if (sem()) {
    // Collect volume fractions and colors from all layers in one field.
    t.fcvf_sum.Reinit(m, 0);
    t.fccl_sum.Reinit(m, kClNone);
    for (auto c : m.Cells()) {
      for (auto l : t.layers) {
        t.fcvf_sum[c] += t.fcvf[l][c];
        if (t.fcvf[l][c] > 0) {
          t.fccl_sum[c] = t.fccl[l][c];
        }
      }
    }
  }
  if (sem.Nested()) {
    Dump(t.fcvf_sum, "vf", "vf_sum", m);
  }
  if (sem.Nested()) {
    Dump(t.fccl_sum, "cl", "cl_sum", m);
  }
}

int main(int argc, const char** argv) {
  MpiWrapper mpi(&argc, &argv);

  ArgumentParser parser(
      "Example for connected component labeling", mpi.IsRoot());
  parser.AddVariable<std::string>("--extra", "")
      .Help("Extra configuration (commands 'set ... ')");
  parser.AddVariable<std::string>("config", "a.conf")
      .Help("Path to configuration file");
  parser.AddVariable<int>("--nx", 128).Help("Mesh size in the x-direction");
  parser.AddVariable<int>("--ny", 128).Help("Mesh size in the y-direction");
  parser.AddVariable<int>("--nz", 1).Help("Mesh size in the z-direction");
  parser.AddVariable<int>("--bs", 32).Help("Block size");
  parser.AddVariable<int>("--layers", 4).Help("Number of layers");

  auto instances = []() {
    auto map = ModuleLabeling<M>::GetInstances();
    std::vector<std::string> res;
    for (auto p : map) {
      res.push_back(p.first);
    }
    return res;
  }();
  parser.AddVariable<std::string>("--labeling", "propagation")
      .Help("Algorithm for connected component labeling")
      .Options(instances);

  auto args = parser.ParseArgs(argc, argv);
  if (const int* p = args.Int.Find("EXIT")) {
    return *p;
  }

  std::stringstream conf;

  MIdx meshsize(args.Int["nx"], args.Int["ny"], args.Int["nz"]);
  const int bs = args.Int["bs"];
  MIdx blocksize(bs, bs, args.Int["nz"] == 1 ? 1 : bs);
  Subdomains<MIdx> sub(meshsize, blocksize, mpi.GetCommSize());
  conf << sub.GetConfig() << '\n';

  const auto configpath = args.String["config"];
  if (!configpath.empty()) {
    conf << "include " + configpath + '\n';
  }
  conf << "set string backend native\n";
  conf << "set double extent 1\n";
  conf << util::Format("set int layers {}\n", args.Int["layers"]);
  conf << util::Format("set string labeling {}\n", args.String["labeling"]);

  conf << args.String["extra"] << '\n';

  return RunMpiBasicString<M>(mpi, Run, conf.str());
}
