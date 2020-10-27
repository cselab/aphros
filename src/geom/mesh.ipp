// Created by Petr Karnakov on 27.10.2020
// Copyright 2020 ETH Zurich

#include "mesh.h"

using M = MeshStructured<double, 3>;

template <class _Scal, size_t _dim>
MeshStructured<_Scal, _dim>::MeshStructured(
    MIdx b, MIdx cs, Rect<Vect> domain, int halos, bool isroot, bool islead,
    MIdx gs, int id)
    // inner
    : blockci_(b, cs)
    , blockfi_(blockci_.GetBegin(), blockci_.GetSize())
    , blockni_(blockci_.GetBegin(), blockci_.GetSize() + MIdx(1))
    // all
    , blockca_(b - MIdx(halos), cs + MIdx(2 * halos))
    , blockfa_(blockca_.GetBegin(), blockca_.GetSize())
    , blockna_(blockca_.GetBegin(), blockca_.GetSize() + MIdx(1))
    // support
    , blockcs_(blockca_.GetBegin() + MIdx(1), blockca_.GetSize() - MIdx(2))
    , blockfs_(blockcs_.GetBegin(), blockcs_.GetSize())
    , blockns_(blockcs_.GetBegin(), blockcs_.GetSize() + MIdx(1))
    // index
    , indexc_(blockca_.GetBegin(), blockca_.GetSize() + MIdx(1))
    , indexf_(indexc_.GetBegin(), indexc_.GetSize())
    , indexn_(indexc_.GetBegin(), indexc_.GetSize())
    , isroot_(isroot)
    , islead_(islead)
    , incells_begin_(blockci_.GetBegin())
    , incells_end_(blockci_.GetEnd())
    , domain_(domain)
    , cell_size_(domain.GetDimensions() / Vect(blockci_.GetSize()))
    , half_cell_size_(cell_size_ * 0.5)
    , cell_volume_(cell_size_.prod())
    , global_size_(gs)
    , id_(id)
    , global_length_(Vect(gs) * cell_size_) {
  static_assert(dim == 3, "Not implemented for dim != 3");

  // mesh step

  // surface area
  face_area_[0] = cell_size_[1] * cell_size_[2];
  face_area_[1] = cell_size_[2] * cell_size_[0];
  face_area_[2] = cell_size_[0] * cell_size_[1];

  // surface vectors
  vs_[0] = Vect(face_area_[0], 0., 0.);
  vs_[1] = Vect(0., face_area_[1], 0.);
  vs_[2] = Vect(0., 0., face_area_[2]);

  // cell neighbour cell offset
  {
    MIdx w = indexc_.GetBegin(); // any cell
    IdxCell c = indexc_.GetIdx(w);
    for (auto q : Nci(c)) {
      MIdx wo;
      switch (q) {
        case 0:
          wo = MIdx(-1, 0, 0);
          break;
        case 1:
          wo = MIdx(1, 0, 0);
          break;
        case 2:
          wo = MIdx(0, -1, 0);
          break;
        case 3:
          wo = MIdx(0, 1, 0);
          break;
        case 4:
          wo = MIdx(0, 0, -1);
          break;
        default:
        case 5:
          wo = MIdx(0, 0, 1);
          break;
      };
      IdxCell cn = indexc_.GetIdx(w + wo);
      cell_cell_[q] = size_t(cn) - size_t(c);
    }
  }

  // cell neighbour face offset
  {
    MIdx w = indexc_.GetBegin(); // any cell
    IdxCell c = indexc_.GetIdx(w);
    for (auto q : Nci(c)) {
      MIdx wo;
      Dir d;
      switch (q) {
        case 0:
          wo = MIdx(0, 0, 0);
          d = Dir::i;
          break;
        case 1:
          wo = MIdx(1, 0, 0);
          d = Dir::i;
          break;
        case 2:
          wo = MIdx(0, 0, 0);
          d = Dir::j;
          break;
        case 3:
          wo = MIdx(0, 1, 0);
          d = Dir::j;
          break;
        case 4:
          wo = MIdx(0, 0, 0);
          d = Dir::k;
          break;
        default:
        case 5:
          wo = MIdx(0, 0, 1);
          d = Dir::k;
          break;
      };
      IdxFace fn = indexf_.GetIdx(std::make_pair(w + wo, d));
      cell_face_[q] = size_t(fn) - size_t(c);
    }
  }

  // cell outward factor
  {
    for (size_t q = 0; q < kCellNumNeighbourFaces; ++q) {
      cell_outward_[q] = (q % 2 == 0 ? -1. : 1.);
    }
  }

  // cell neighbour node offset
  {
    MIdx w = indexc_.GetBegin(); // any cell
    IdxCell c = indexc_.GetIdx(w);
    auto qm = kCellNumNeighbourNodes;
    for (size_t q = 0; q < qm; ++q) {
      MIdx wo;
      switch (q) {
        case 0:
          wo = MIdx(0, 0, 0);
          break;
        case 1:
          wo = MIdx(1, 0, 0);
          break;
        case 2:
          wo = MIdx(0, 1, 0);
          break;
        case 3:
          wo = MIdx(1, 1, 0);
          break;
        case 4:
          wo = MIdx(0, 0, 1);
          break;
        case 5:
          wo = MIdx(1, 0, 1);
          break;
        case 6:
          wo = MIdx(0, 1, 1);
          break;
        default:
        case 7:
          wo = MIdx(1, 1, 1);
          break;
      };
      IdxNode nn = indexn_.GetIdx(w + wo);
      cell_node_[q] = size_t(nn) - size_t(c);
    }
  }

  // face neighbour cell offset
  {
    const MIdx w = indexc_.GetBegin(); // any cell
    for (size_t d = 0; d < dim; ++d) {
      IdxFace f = indexf_.GetIdx(w, Dir(d));
      auto qm = kFaceNumNeighbourCells;
      for (size_t q = 0; q < qm; ++q) {
        MIdx wo(0);
        switch (q) {
          case 0:
            wo[d] = -1;
            break;
          default:
          case 1:
            wo[d] = 0;
            break;
        };
        IdxCell cn = indexc_.GetIdx(w + wo);
        face_cell_[d * qm + q] = size_t(cn) - size_t(f);
      }
    }
  }

  // face neighbour face offset
  {
    const MIdx w = indexc_.GetBegin(); // any cell
    for (size_t d = 0; d < dim; ++d) {
      const IdxFace f = indexf_.GetIdx(w, Dir(d));
      auto qm = kCellNumNeighbourFaces;
      for (size_t q = 0; q < qm; ++q) {
        MIdx wo(0);
        wo[q / 2] = (q % 2 == 0 ? -1 : 1);
        const IdxFace fn = indexf_.GetIdx(w + wo, Dir(d));
        face_face_[d * qm + q] = size_t(fn) - size_t(f);
      }
    }
  }

  // face neighbour cell offset
  {
    const MIdx w = indexc_.GetBegin(); // any cell
    for (size_t d = 0; d < dim; ++d) {
      IdxFace f = indexf_.GetIdx(w, Dir(d));
      auto qm = kFaceNumNeighbourNodes;
      for (size_t q = 0; q < qm; ++q) {
        MIdx wo;
        if (d == 0) {
          switch (q) {
            case 0:
              wo = MIdx(0, 0, 0);
              break;
            case 1:
              wo = MIdx(0, 1, 0);
              break;
            case 2:
              wo = MIdx(0, 1, 1);
              break;
            default:
            case 3:
              wo = MIdx(0, 0, 1);
              break;
          }
        } else if (d == 1) {
          switch (q) {
            case 0:
              wo = MIdx(0, 0, 0);
              break;
            case 1:
              wo = MIdx(0, 0, 1);
              break;
            case 2:
              wo = MIdx(1, 0, 1);
              break;
            default:
            case 3:
              wo = MIdx(1, 0, 0);
              break;
          }
        } else {
          switch (q) {
            case 0:
              wo = MIdx(0, 0, 0);
              break;
            case 1:
              wo = MIdx(1, 0, 0);
              break;
            case 2:
              wo = MIdx(1, 1, 0);
              break;
            default:
            case 3:
              wo = MIdx(0, 1, 0);
              break;
          }
        }
        IdxNode nn = indexn_.GetIdx(w + wo);
        face_node_[d * qm + q] = size_t(nn) - size_t(f);
      }
    }
  }
  // stencil 3x3x3 offsets
  {
    const IdxCell c = indexc_.GetIdx(blockci_.GetBegin());
    size_t i = 0;
    for (auto cn : StencilGeneral<1>(c)) {
      stencil_[i++] = size_t(cn) - size_t(c);
    }
    fassert_equal(i, kNumStencil);
  }
  // stencil 5x5x5 offsets
  {
    const IdxCell c = indexc_.GetIdx(blockci_.GetBegin());
    size_t i = 0;
    for (auto cn : StencilGeneral<2>(c)) {
      stencil5_[i++] = size_t(cn) - size_t(c);
    }
    fassert_equal(i, kNumStencil5);
  }
}

template <class M>
M InitUniformMesh(
    Rect<typename M::Vect> domain, typename M::MIdx begin, typename M::MIdx s,
    int halos, bool isroot, bool islead, typename M::MIdx gs, int id) {
  return M(begin, s, domain, halos, isroot, islead, gs, id);
}
