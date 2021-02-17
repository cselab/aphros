// Created by Petr Karnakov on 17.02.2021
// Copyright 2021 ETH Zurich

#include "idx.h"
#include "notation.h"

namespace geom {

template <class M>
struct Loop {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using MIdx = typename M::MIdx;
  static constexpr size_t dim = M::dim;

  struct Face;

  struct Cell {
    static IdxCell GetIdx(MIdx w, MIdx lead) {
      return IdxCell(w[0] + lead[0] * w[1] + lead[0] * lead[1] * w[2]);
    }
    Cell() = default;
    Cell(MIdx w_, MIdx lead_, Vect h_)
        : w(w_)
        , lead(lead_)
        , h(h_)
        , c(GetIdx(w, lead))
        , volume(h.prod())
        , center((Vect(w) + Vect(0.5)) * h) {}
    operator IdxCell() const {
      return c;
    }
    inline Face face(IdxNci q) const;
    inline Scal outward_factor(IdxNci q) const {
      return q.raw() % 2 == 0 ? -1 : 1;
    }

    MIdx w;
    MIdx lead;
    Vect h;
    IdxCell c;
    Scal volume;
    Vect center;
  };

  struct Face {
    static IdxFace GetIdx(size_t direction, MIdx w, MIdx lead) {
      return IdxFace(
          direction * lead.prod() + w[0] + lead[0] * w[1] +
          lead[0] * lead[1] * w[2]);
    }
    Face() = default;
    Face(size_t direction_, MIdx w_, MIdx lead_, Vect h_)
        : direction(direction_)
        , w(w_)
        , lead(lead_)
        , h(h_)
        , f(GetIdx(direction, w, lead))
        , area(h.prod() / h[direction]) {}
    operator IdxFace() const {
      return f;
    }
    inline Cell cell(Side s) const;

    size_t direction;
    MIdx w;
    MIdx lead;
    Vect h;
    IdxFace f;
    Scal area;
  };

  template <class Func>
  static void ForEachCell(const M& m, Func func) {
    const auto h = m.GetCellSize();
    const auto& indexer = m.GetIndexCells();
    const auto& block = m.GetInBlockCells();
    const MIdx lead = indexer.GetSize();
    const MIdx mstart = block.GetBegin();
    const MIdx mend = block.GetEnd();
    for (int iz = mstart[2]; iz < mend[2]; ++iz) {
      for (int iy = mstart[1]; iy < mend[1]; ++iy) {
        for (int ix = mstart[0]; ix < mend[0]; ++ix) {
          func(Cell(MIdx(ix, iy, iz), lead, h));
        }
      }
    }
  }

  template <class Func>
  static void ForEachFace(const M& m, Func func) {
    const auto h = m.GetCellSize();
    const auto& indexer = m.GetIndexCells();
    const auto& block = m.GetInBlockCells();
    const MIdx lead = indexer.GetSize();
    const MIdx mstart = block.GetBegin();
    const MIdx mend = block.GetEnd();
    for (size_t d = 0; d < dim; ++d) {
      for (int iz = mstart[2]; iz < mend[2] + (d == 2 ? 1 : 0); ++iz) {
        for (int iy = mstart[1]; iy < mend[1] + (d == 1 ? 1 : 0); ++iy) {
          for (int ix = mstart[0]; ix < mend[0] + (d == 0 ? 1 : 0); ++ix) {
            func(Face(d, MIdx(ix, iy, iz), lead, h));
          }
        }
      }
    }
  }
};

template <class M>
auto Loop<M>::Cell::face(IdxNci q) const -> Face {
  switch (q.raw()) {
    case 0:
      return Face(0, w, lead, h);
    case 1:
      return Face(0, w + MIdx::GetUnit(0), lead, h);
    case 2:
      return Face(1, w, lead, h);
    case 3:
      return Face(1, w + MIdx::GetUnit(1), lead, h);
    case 4:
      return Face(2, w, lead, h);
    case 5:
      return Face(2, w + MIdx::GetUnit(2), lead, h);
  }
  return {};
}

template <class M>
auto Loop<M>::Face::cell(Side s) const -> Cell {
  switch (s.raw()) {
    case 0:
      return Cell(w - MIdx::GetUnit(direction), lead, h);
    case 1:
      return Cell(w, lead, h);
  }
  return {};
}

} // namespace geom
