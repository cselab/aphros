// Created by Petr Karnakov on 22.03.2020
// Copyright 2020 ETH Zurich

#include <array>
#include <iostream>
#include <stdexcept>
#include <string>

#include <inside.h>
#include <util/posthook.h>

template <class M>
void InitEmbedHook(
    FieldNode<typename M::Scal>& fnl, const Vars& var, const M& m) {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  if (m.IsRoot()) {
    std::cout << "Embedded boundaries from lib/inside" << std::endl;
  }

  int nt;
  int nv;
  int* tri;
  double* ver;
  struct Inside* inside_state;

  const std::string path = var.String["inside_path"];
  if (inside_mesh_read(path.c_str(), &nt, &tri, &nv, &ver) != 0) {
    throw std::runtime_error(FILELINE + ": fail to read mesh '" + path + "'");
  }

  // normalize ver to [0,1]
  {
    // bounding box
    Vect box0(std::numeric_limits<Scal>::max());
    Vect box1(-std::numeric_limits<Scal>::max());
    for (int i = 0; i < nv * 3; ++i) {
      const int d = i % 3;
      box0[d] = std::min<Scal>(box0[d], ver[i]);
      box1[d] = std::max<Scal>(box1[d], ver[i]);
    }
    const Scal extent = (box1 - box0).abs().max();
    const Vect center = (box0 + box1) * 0.5;
    if (m.IsRoot()) {
      std::cout << "lib/inside"
                << " bounding box: " << box0 << " " << box1
                << " center: " << center << std::endl;
    }
    for (int i = 0; i < nv * 3; ++i) {
      const int d = i % 3;
      auto& x = ver[i];
      x = 0.5 + (x - center[d]) / extent;
    }
  }

  inside_ini(nt, tri, ver, &inside_state);
  auto distance = [&](Vect x) -> Scal {
    double p[3] = {x[0], x[1], x[2]};
    return inside_distance(inside_state, p);
  };
  auto inside = [&](Vect x) -> bool {
    double p[3] = {x[0], x[1], x[2]};
    return inside_inside(inside_state, p);
  };

  fnl.Reinit(m, GetNan<Scal>());
  for (auto c : m.AllCells()) {
    const size_t num_nodes = m.GetNumNodes(c);
    size_t num_inside = 0;
    for (size_t i = 0; i < num_nodes; ++i) {
      const IdxNode n = m.GetNode(c, i);
      if (inside(m.GetNode(n))) {
        ++num_inside;
      }
    }
    // levelset > 0 inside body
    const Scal inf = m.GetCellSize()[0] * 10; // large value
    if (num_inside == 0) { // whole cell outside
      for (size_t i = 0; i < num_nodes; ++i) {
        const IdxNode n = m.GetNode(c, i);
        if (IsNan(fnl[n])) {
          fnl[n] = -inf;
        }
      }
    } else if (num_inside == num_nodes) { // whole cell inside
      for (size_t i = 0; i < num_nodes; ++i) {
        const IdxNode n = m.GetNode(c, i);
        if (IsNan(fnl[n])) {
          fnl[n] = inf;
        }
      }
    } else { // cell crosses surface
      for (size_t i = 0; i < num_nodes; ++i) {
        const IdxNode n = m.GetNode(c, i);
        fnl[n] = -distance(m.GetNode(n));
      }
    }
  }
  if (var.Int["eb_init_inverse"]) {
    for (auto n : m.AllNodes()) {
      fnl[n] = -fnl[n];
    }
  }

  inside_fin(inside_state);
  inside_mesh_fin(tri, ver);
}

using M = MeshStructured<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;

template void InitEmbedHook(FieldNode<Scal>&, const Vars&, const M&);
