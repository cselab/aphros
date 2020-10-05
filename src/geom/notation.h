#pragma once

#include <iosfwd>
#include <string>

#include "vect.h"

namespace generic {

template <class Scal, size_t dim_>
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
  Scal sign() const {
    return forward_ == 1 ? 1 : -1;
  }
  Direction operator>>(int shift) const {
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
  bool operator==(const Direction& other) const {
    return other.index_ == index_ && other.forward_ == forward_;
  }
  bool operator!=(const Direction& other) const {
    return !(*this == other);
  }
  friend std::ostream& operator<<(std::ostream& out, const Direction& d) {
    out << d.index_ << ' ' << d.forward_;
    return out;
  }
  size_t nci() const {
    return index_ * 2 + forward_;
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
  using MIdx = typename M::MIdx;
  using Direction = generic::Direction<Scal, dim>;

  IdxCellMesh(IdxCell c, const M& m_) : center(*this), c_(c), m(m_) {}
  IdxCellMesh(const M& m_) : IdxCellMesh(IdxCell(0), m_) {}
  IdxCellMesh operator+(Direction d) const {
    return m(m.GetCell(c_, d.nci()));
  }
  IdxCellMesh operator-(Direction d) const {
    return m(m.GetCell(c_, (-d).nci()));
  }
  friend std::ostream& operator<<(std::ostream& out, const IdxCellMesh& c) {
    out << c.m.GetIndexCells().GetMIdx(c.c_);
    return out;
  }
  friend std::istream& operator>>(std::istream& in, IdxCellMesh& c) {
    MIdx w;
    in >> w;
    c.c_ = c.m.GetIndexCells().GetIdx(w);
    return in;
  }
  operator IdxCell() const {
    return c_;
  }
  explicit operator MIdx() const {
    return m.GetIndexCells().GetMIdx(c_);
  }
  IntIdx operator[](size_t i) const {
    return MIdx(*this)[i];
  }
  auto face(Direction d) const {
    return m(m.GetFace(c_, d.nci()));
  }

  class LazyCenter {
   public:
    LazyCenter(const IdxCellMesh& owner) : owner_(owner) {}
    operator Vect() const {
      return owner_.m.GetCenter(owner_.c_);
    }
    Vect operator()() const {
      return Vect(*this);
    }
    Scal operator[](size_t i) const {
      return Vect(*this)[i];
    }

   private:
    const IdxCellMesh& owner_;
  };
  LazyCenter center;

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
  using MIdx = typename M::MIdx;
  using Direction = generic::Direction<Scal, dim>;

  IdxFaceMesh(IdxFace f, const M& m_) : center(*this), cm(*this), cp(*this), f_(f), m(m_) {}
  IdxFaceMesh(const M& m_) : IdxFaceMesh(IdxFace(0), m_) {}
  IdxFaceMesh operator+(Direction d) const {
    return m(m.GetFace(f_, d.GetIndex() * 2 + d.GetForward()));
  }
  IdxFaceMesh operator-(Direction d) const {
    return m(m.GetFace(f_, d.GetIndex() * 2 + 1 - d.GetForward()));
  }
  friend std::ostream& operator<<(std::ostream& out, const IdxFaceMesh& f) {
    auto& index = f.m.GetIndexFaces();
    out << index.GetMIdx(f.f_) << index.GetDir(f.f_).GetLetter();
    return out;
  }
  friend std::istream& operator>>(std::istream& in, IdxFaceMesh& f) {
    MIdx w;
    std::string sdir;
    in >> w >> sdir;
    using Dir = typename M::Dir;
    Dir dir;
    if (sdir == "x") {
      dir = Dir(0);
    } else if (sdir == "y") {
      dir = Dir(1);
    } else if (sdir == "z") {
      dir = Dir(2);
    } else {
      in.setstate(std::ios_base::failbit);
      return in;
    }
    f.f_ = f.m.GetIndexFaces().GetIdx(w, dir);
    return in;
  }
  operator IdxFace() const {
    return f_;
  }
  explicit operator MIdx() const {
    return m.GetIndexFaces().GetMIdx(f_);
  }
  IntIdx operator[](size_t i) const {
    return MIdx(*this)[i];
  }
  operator Direction() const {
    return Direction(size_t(m.GetIndexFaces().GetDir(f_)));
  }
  Direction direction() const {
    return Direction(*this);
  }

  class LazyCenter {
   public:
    LazyCenter(const IdxFaceMesh& owner) : owner_(owner) {}
    operator Vect() const {
      return owner_.m.GetCenter(owner_.f_);
    }
    Vect operator()() const {
      return Vect(*this);
    }
    Scal operator[](size_t i) const {
      return Vect(*this)[i];
    }

   private:
    const IdxFaceMesh& owner_;
  };
  LazyCenter center;

  template <size_t id>
  class LazyCell {
   public:
    LazyCell(const IdxFaceMesh& owner) : owner_(owner) {}
    operator IdxCellMesh<M>() const {
      return owner_.m(owner_.m.GetCell(owner_.f_, id));
    }
    operator IdxCell() const {
      return owner_.m.GetCell(owner_.f_, id);
    }
    IdxCellMesh<M> operator()() const {
      return IdxCellMesh<M>(*this);
    }

   private:
    const IdxFaceMesh& owner_;
  };
  const LazyCell<0> cm;
  const LazyCell<1> cp;

 private:
  IdxFace f_;
  const M& m;
};
