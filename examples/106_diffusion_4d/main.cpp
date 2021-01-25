// Created by Petr Karnakov on 24.12.2020
// Copyright 2020 ETH Zurich

#include <cassert>
#include <iostream>
#include <memory>
#include <string>

#include <dump/dump.h>
#include <geom/mesh.h>
#include <parse/argparse.h>
#include <parse/parser.h>
#include <solver/approx_eb.h>
#include <solver/cond.h>
#include <solver/embed.h>
#include <util/format.h>

using M = MeshCartesian<double, 4>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using MIdx = typename M::MIdx;

void Run(M& m, Vars& var) {
  FieldCell<Scal> fcu(m); // field
  MapEmbed<BCond<Scal>> mebc; // face conditions

  const Scal dt = var.Double["dt"];

  // initial
  fcu.Reinit(m, 0);
  for (auto c : m.AllCellsM()) {
    const Vect xc = m.GetGlobalLength() * 0.5;
    fcu[c] = xc.sqrdist(c.center) < 0.2 ? 1 : 0;
  }

  for (int step = 0; step < var.Int["nsteps"]; ++step) {
    const auto ffg = UEmbed<M>::Gradient(fcu, mebc, m);
    FieldFace<Scal> ff_flux(m);
    for (auto f : m.FacesM()) {
      ff_flux[f] = ffg[f] * f.area;
    }

    for (auto c : m.CellsM()) {
      Scal sum = 0;
      for (auto q : m.Nci(c)) {
        sum += ff_flux[c.face(q)] * c.outward_factor(q);
      }
      fcu[c] += sum * dt / c.volume;
    }

    Scal sum0 = 0;
    Scal sum2 = 0;
    for (auto c : m.CellsM()) {
      const Vect xc = m.GetGlobalLength() * 0.5;
      sum0 += fcu[c] * c.volume;
      sum2 += fcu[c] * xc.sqrdist(c.center) * c.volume;
    }
    std::cout << util::Format("{}) sum0={:.3f} sum2={:.3f}\n", step, sum0, sum2)
              << std::endl;

    std::vector<std::pair<std::string, std::vector<Scal>>> data;
    for (auto d : m.dirs) {
      std::vector<Scal> field;
      for (auto c : m.CellsM()) {
        field.push_back(c.center[d]);
      }
      const std::string name = std::string() + M::Dir(d).letter();
      data.emplace_back(name, field);
    }

    {
      std::vector<Scal> field;
      for (auto c : m.Cells()) {
        field.push_back(fcu[c]);
      }
      const std::string name = "u";
      data.emplace_back(name, field);
    }

    DumpCsv<Scal>(data, util::Format("u_{:04d}.csv", step));
  }
}

int main(int argc, const char** argv) {
  ArgumentParser argparser("Diffusion solver in 4D");
  argparser.AddVariable<std::string>("config", "a.conf")
      .Help("Path to configuration file");

  auto args = argparser.ParseArgs(argc, argv);
  if (const int* p = args.Int.Find("EXIT")) {
    return *p;
  }

  Vars var;
  Parser parser(var);
  parser.ParseFile(args.String["config"]);

  const Rect<Vect> dom(Vect(0), Vect(1));
  const MIdx begin(0);
  const MIdx size(var.Int["nx"], var.Int["ny"], var.Int["nz"], var.Int["nw"]);
  const int halos = 2;
  M m = InitUniformMesh<M>(dom, begin, size, halos, true, true, size, 0);
  Run(m, var);
}
