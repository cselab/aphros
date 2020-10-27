// Created by Petr Karnakov on 27.10.2020
// Copyright 2020 ETH Zurich

#include "mesh.h"

using M = MeshStructured<double, 3>;

template <class _Scal, size_t _dim>
MeshStructured<_Scal, _dim>::MeshStructured(
    MIdx b, MIdx cs, Rect<Vect> dom, int halos, bool isroot, bool islead, MIdx gs,
    int id)
    // inner
    : bci_(b, cs)
    , bfi_(bci_.GetBegin(), bci_.GetSize())
    , bni_(bci_.GetBegin(), bci_.GetSize() + MIdx(1))
    // all
    , bca_(b - MIdx(halos), cs + MIdx(2 * halos))
    , bfa_(bca_.GetBegin(), bca_.GetSize())
    , bna_(bca_.GetBegin(), bca_.GetSize() + MIdx(1))
    // support
    , bcs_(bca_.GetBegin() + MIdx(1), bca_.GetSize() - MIdx(2))
    , bfs_(bcs_.GetBegin(), bcs_.GetSize())
    , bns_(bcs_.GetBegin(), bcs_.GetSize() + MIdx(1))
    // raw
    , bcr_(bca_.GetBegin(), bca_.GetSize() + MIdx(1))
    , bfr_(bcr_.GetBegin(), bcr_.GetSize())
    , bnr_(bcr_.GetBegin(), bcr_.GetSize())
    , isroot_(isroot)
    , islead_(islead)
    , mb_(bci_.GetBegin())
    , me_(bci_.GetEnd())
    , dom_(dom)
    , h_(dom.GetDimensions() / Vect(bci_.GetSize()))
    , hh_(h_ * 0.5)
    , vol_(h_.prod())
    , gs_(gs)
    , id_(id)
    , gl_(Vect(gs) * h_)
    , checknan_(false)
    , edim_(dim) {
  static_assert(dim == 3, "Not implemented for dim != 3");

  // mesh step

  // surface area
  va_[0] = h_[1] * h_[2];
  va_[1] = h_[2] * h_[0];
  va_[2] = h_[0] * h_[1];

  // surface vectors
  vs_[0] = Vect(va_[0], 0., 0.);
  vs_[1] = Vect(0., va_[1], 0.);
  vs_[2] = Vect(0., 0., va_[2]);

  // cell neighbour cell offset
  {
    MIdx w = bcr_.GetBegin(); // any cell
    IdxCell c = bcr_.GetIdx(w);
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
      IdxCell cn = bcr_.GetIdx(w + wo);
      cnc_[q] = size_t(cn) - size_t(c);
    }
  }

  // cell neighbour face offset
  {
    MIdx w = bcr_.GetBegin(); // any cell
    IdxCell c = bcr_.GetIdx(w);
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
      IdxFace fn = bfr_.GetIdx(std::make_pair(w + wo, d));
      cnf_[q] = size_t(fn) - size_t(c);
    }
  }

  // cell outward factor
  {
    for (size_t q = 0; q < kCellNumNeighbourFaces; ++q) {
      co_[q] = (q % 2 == 0 ? -1. : 1.);
    }
  }

  // cell neighbour node offset
  {
    MIdx w = bcr_.GetBegin(); // any cell
    IdxCell c = bcr_.GetIdx(w);
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
      IdxNode nn = bnr_.GetIdx(w + wo);
      cnn_[q] = size_t(nn) - size_t(c);
    }
  }

  // face neighbour cell offset
  {
    const MIdx w = bcr_.GetBegin(); // any cell
    for (size_t d = 0; d < dim; ++d) {
      IdxFace f = bfr_.GetIdx(w, Dir(d));
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
        IdxCell cn = bcr_.GetIdx(w + wo);
        fnc_[d * qm + q] = size_t(cn) - size_t(f);
      }
    }
  }

  // face neighbour face offset
  {
    const MIdx w = bcr_.GetBegin(); // any cell
    for (size_t d = 0; d < dim; ++d) {
      const IdxFace f = bfr_.GetIdx(w, Dir(d));
      auto qm = kCellNumNeighbourFaces;
      for (size_t q = 0; q < qm; ++q) {
        MIdx wo(0);
        wo[q / 2] = (q % 2 == 0 ? -1 : 1);
        const IdxFace fn = bfr_.GetIdx(w + wo, Dir(d));
        fnf_[d * qm + q] = size_t(fn) - size_t(f);
      }
    }
  }

  // face neighbour cell offset
  {
    const MIdx w = bcr_.GetBegin(); // any cell
    for (size_t d = 0; d < dim; ++d) {
      IdxFace f = bfr_.GetIdx(w, Dir(d));
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
        IdxNode nn = bnr_.GetIdx(w + wo);
        fnn_[d * qm + q] = size_t(nn) - size_t(f);
      }
    }
  }
  // stencil 3x3x3 offsets
  {
    const IdxCell c = bcr_.GetIdx(bci_.GetBegin());
    size_t i = 0;
    for (auto cn : StencilGeneral<1>(c)) {
      stencil_[i++] = size_t(cn) - size_t(c);
    }
    fassert_equal(i, kNumStencil);
  }
  // stencil 5x5x5 offsets
  {
    const IdxCell c = bcr_.GetIdx(bci_.GetBegin());
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
