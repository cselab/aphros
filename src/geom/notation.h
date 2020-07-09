#pragma once

#include "vect.h"

namespace generic {

template <class Scal, size_t dim>
class Direction {
 public:
  static constexpr size_t dim = dim_;
  using Vect = generic::Vect<Scal, dim>;

  explicit Direction(size_t index) : Direction(index, 1) {}
  Direction(size_t index, size_t forward)
      : index_(index % dim), forward_(forward) {}
  operator size_t() const {
    return index_;
  }
  Direction orient(const Vect& v) const {
    return Direction(index_, v[index_] > 0 ? 1 : 0);
  }
  Direction operator>>(size_t shift) const {
    return Direction((index_ + shift) % dim, forward_);
  }
  Direction operator-() const {
    return Direction(index_, 1 - forward_);
  }
  size_t GetIndex() const {
    return index_;
  }
  size_t GetForward() const {
    return forward_;
  }

 private:
  const size_t index_; // direction, [0, dim - 1]
  const size_t forward_; // 0: backward, 1: forward
};

} // namespace generic

template <class M>
class IdxCellMesh {
 public:
  static constexpr size_t dim = M::dim;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using Direction = generic::Direction<>;

  IdxCellMesh(IdxCell c, const M& m) : c_(c), m(m) {}
  IdxCellMesh operator+(Direction d) const {
    return IdxCellMesh(m.GetCell(c_, d.GetIndex() * 2 + d.GetForward()), m);
  }
  IdxCellMesh operator-(Direction d) const {
    return IdxCellMesh(m.GetCell(c_, d.GetIndex() * 2 + 1 - d.GetForward()), m);
  }

 private:
  IdxCell c_;
  const M& m;
};

template <class M>
class IdxFaceMesh {
 public:
  static constexpr size_t dim = M::dim;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using Direction = generic::Direction<>;

  IdxFaceMesh(IdxFace f, const M& m) : f_(f), m(m) {}
  IdxFaceMesh operator+(Direction d) const {
    return m.GetFace(
        m.GetCell(m.GetCell(f_, 0), d.GetIndex() * 2 + d.GetForward()), 1);
  }
  IdxFaceMesh operator-(Direction d) const {
    return IdxFaceMesh(m.GetCell(c_, d.GetIndex() * 2 + 1 - d.GetForward()), m);
  }

 private:
  IdxFace f_;
  size_t index_;
  const M& m;
};

