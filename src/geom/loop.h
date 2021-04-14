// Created by Petr Karnakov on 17.02.2021
// Copyright 2021 ETH Zurich

#pragma once

#include "idx.h"
#include "notation.h"

namespace geom {

template <class M>
struct Loop {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using MIdx = typename M::MIdx;
  using MIdx2 = generic::MIdx<2>;
  using MIdx3 = generic::MIdx<3>;
  using MIdx4 = generic::MIdx<4>;
  static constexpr size_t dim = M::dim;

  struct Face;

  struct Cell {
    static IdxCell GetIdx(MIdx w, MIdx lead) {
      return IdxCell(w.dot(lead.cumprod()));
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
    inline Face face(size_t q) const {
      return face(IdxNci(q));
    }
    inline Cell cell(IdxNci q) const;
    inline Cell cell(size_t q) const {
      return cell(IdxNci(q));
    }
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
      return IdxFace(direction * lead.prod() + w.dot(lead.cumprod()));
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
    inline Cell cell(size_t s) const {
      return cell(Side(s));
    }

    size_t direction;
    MIdx w;
    MIdx lead;
    Vect h;
    IdxFace f;
    Scal area;
  };

  template <class Func>
  static void ForEachCell(const M& m, Func func, MIdx3*) {
    const auto h = m.GetCellSize();
    const auto& indexer = m.GetIndexCells();
    const auto& block = m.GetInBlockCells();
    const MIdx lead = indexer.GetSize();
    const MIdx start = indexer.GetBegin();
    const MIdx bstart = block.GetBegin() - start;
    const MIdx bend = block.GetEnd() - start;
    for (int iz = bstart[2]; iz < bend[2]; ++iz) {
      for (int iy = bstart[1]; iy < bend[1]; ++iy) {
        for (int ix = bstart[0]; ix < bend[0]; ++ix) {
          func(Cell(MIdx(ix, iy, iz), lead, h));
        }
      }
    }
  }

  template <class Func>
  static void ForEachFace(const M& m, Func func, MIdx3*) {
    const auto h = m.GetCellSize();
    const auto& indexer = m.GetIndexCells();
    const auto& block = m.GetInBlockCells();
    const MIdx lead = indexer.GetSize();
    const MIdx start = indexer.GetBegin();
    const MIdx bstart = block.GetBegin() - start;
    const MIdx bend = block.GetEnd() - start;
    for (size_t d = 0; d < dim; ++d) {
      for (int iz = bstart[2]; iz < bend[2] + (d == 2 ? 1 : 0); ++iz) {
        for (int iy = bstart[1]; iy < bend[1] + (d == 1 ? 1 : 0); ++iy) {
          for (int ix = bstart[0]; ix < bend[0] + (d == 0 ? 1 : 0); ++ix) {
            func(Face(d, MIdx(ix, iy, iz), lead, h));
          }
        }
      }
    }
  }

  template <class Func>
  static void ForEachCell(const M& m, Func func, MIdx2*) {
    const auto h = m.GetCellSize();
    const auto& indexer = m.GetIndexCells();
    const auto& block = m.GetInBlockCells();
    const MIdx lead = indexer.GetSize();
    const MIdx start = indexer.GetBegin();
    const MIdx bstart = block.GetBegin() - start;
    const MIdx bend = block.GetEnd() - start;
    for (int iy = bstart[1]; iy < bend[1]; ++iy) {
      for (int ix = bstart[0]; ix < bend[0]; ++ix) {
        func(Cell(MIdx(ix, iy), lead, h));
      }
    }
  }

  template <class Func>
  static void ForEachFace(const M& m, Func func, MIdx2*) {
    const auto h = m.GetCellSize();
    const auto& indexer = m.GetIndexCells();
    const auto& block = m.GetInBlockCells();
    const MIdx lead = indexer.GetSize();
    const MIdx start = indexer.GetBegin();
    const MIdx bstart = block.GetBegin() - start;
    const MIdx bend = block.GetEnd() - start;
    for (size_t d = 0; d < dim; ++d) {
      for (int iy = bstart[1]; iy < bend[1] + (d == 1 ? 1 : 0); ++iy) {
        for (int ix = bstart[0]; ix < bend[0] + (d == 0 ? 1 : 0); ++ix) {
          func(Face(d, MIdx(ix, iy), lead, h));
        }
      }
    }
  }

  template <class Func>
  static void ForEachCell(const M& m, Func func) {
    ForEachCell(m, func, (MIdx*)nullptr);
  }
  template <class Func>
  static void ForEachFace(const M& m, Func func) {
    ForEachFace(m, func, (MIdx*)nullptr);
  }
};

template <class M>
auto Loop<M>::Cell::face(IdxNci q) const -> Face {
  const size_t d = q.raw() / 2;
  if (q.raw() % 2 == 0) {
    return Face(d, w, lead, h);
  } else {
    return Face(d, w + MIdx::GetUnit(d), lead, h);
  }
}

template <class M>
auto Loop<M>::Cell::cell(IdxNci q) const -> Cell {
  const size_t d = q.raw() / 2;
  if (q.raw() % 2 == 0) {
    return Cell(w - MIdx::GetUnit(d), lead, h);
  } else {
    return Cell(w + MIdx::GetUnit(d), lead, h);
  }
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
