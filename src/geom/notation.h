// Created by Petr Karnakov on 09.07.2020
// Copyright 2020 ETH Zurich

#pragma once

#include <iosfwd>
#include <string>

#include "vect.h"

template <class U, class T>
constexpr ptrdiff_t GetOffset(T U::*member) {
  return (char*)&((U*)nullptr->*member) - (char*)nullptr;
}

template <class U, class T>
constexpr U* GetOwner(T* ptr, T U::*member) {
  return (U*)((char*)ptr - GetOffset(member));
}

class Side {
 public:
  Side(size_t value) : value_(value) {}
  size_t raw() const {
    return value_;
  }
  bool operator==(const Side& other) const {
    return value_ == other.value_;
  }
  bool operator!=(const Side& other) const {
    return !((*this) == other);
  }
  Side opposite() const {
    return Side(1 - value_);
  }

 private:
  size_t value_; // 0: backward, 1: forward
};

namespace generic {

template <size_t dim_>
class Direction {
 public:
  static constexpr size_t dim = dim_;

  explicit Direction(size_t index) : Direction(index, 1) {}
  Direction(size_t index, Side side) : index_(index % dim), side_(side) {}
  operator size_t() const {
    return index_;
  }
  template <class T>
  Direction orient(const generic::Vect<T, dim>& v) const {
    return Direction(index_, v[index_] > 0 ? 1 : 0);
  }
  int sign() const {
    return -1 + side_.raw() * 2;
  }
  Direction next(int shift = 1) const {
    return Direction(index_ + shift, side_);
  }
  Direction operator-() const {
    return Direction(index_, 1 - side_.raw());
  }
  size_t index() const {
    return index_;
  }
  Side side() const {
    return side_;
  }
  bool operator==(const Direction& other) const {
    return other.index_ == index_ && other.side_ == side_;
  }
  bool operator!=(const Direction& other) const {
    return !(*this == other);
  }
  friend std::ostream& operator<<(std::ostream& out, const Direction& d) {
    out << d.index_ << ' ' << d.side_.raw();
    return out;
  }
  IdxNci nci() const {
    return IdxNci(index_ * 2 + side_.raw());
  }

 private:
  size_t index_; // direction, [0, dim - 1]
  Side side_; // 0: backward, 1: forward
};

} // namespace generic

template <class M>
class IdxCellMesh {
 public:
  static constexpr size_t dim = M::dim;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using MIdx = typename M::MIdx;
  using Direction = generic::Direction<dim>;

  IdxCellMesh(IdxCell c, const M& m_) : idxcell_(c), m(m_) {}
  IdxCellMesh operator+(Direction d) const {
    return m(m.GetCell(idxcell_, d.nci()));
  }
  IdxCellMesh operator-(Direction d) const {
    return m(m.GetCell(idxcell_, (-d).nci()));
  }
  friend std::ostream& operator<<(std::ostream& out, const IdxCellMesh& c) {
    out << c.m.GetIndexCells().GetMIdx(c.idxcell_);
    return out;
  }
  friend std::istream& operator>>(std::istream& in, IdxCellMesh& c) {
    MIdx w;
    in >> w;
    c.idxcell_ = c.m.GetIndexCells().GetIdx(w);
    return in;
  }
  operator IdxCell() const {
    return idxcell_;
  }
  explicit operator MIdx() const {
    return m.GetIndexCells().GetMIdx(idxcell_);
  }
  IntIdx operator[](size_t i) const {
    return MIdx(*this)[i];
  }
  auto face(Direction d) const {
    return m(m.GetFace(idxcell_, d.nci()));
  }
  auto face(IdxNci q) const {
    return m(m.GetFace(idxcell_, q));
  }
  int outward_factor(IdxNci q) const {
    return m.GetOutwardFactor(idxcell_, q);
  }
  Vect outward_surface(IdxNci q) const {
    return m.GetOutwardSurface(idxcell_, q);
  }

  class LazyCenter {
   public:
    operator Vect() const {
      auto owner = GetOwner(this, &IdxCellMesh::center);
      return owner->m.GetCenter(owner->idxcell_);
    }
    Vect operator()() const {
      return Vect(*this);
    }
    Scal operator[](size_t i) const {
      return Vect(*this)[i];
    }

   private:
    friend IdxCellMesh;
    // After construction of IdxCellMesh, LazyCenter cannot be
    // copied or moved since GetOwner() relies on the relative address.
    LazyCenter(const LazyCenter&) = default;
    LazyCenter& operator=(const LazyCenter&) = default;
  };
  const LazyCenter center{};

  class LazyVolume {
   public:
    operator Scal() const {
      auto owner = GetOwner(this, &IdxCellMesh::volume);
      return owner->m.GetVolume(owner->idxcell_);
    }
    Scal operator()() const {
      return Scal(*this);
    }

   private:
    friend IdxCellMesh;
    LazyVolume(const LazyVolume&) = default;
    LazyVolume& operator=(const LazyVolume&) = default;
  };
  const LazyVolume volume{};

 private:
  IdxCell idxcell_;
  const M& m;
};

template <class M>
class IdxFaceMesh {
 public:
  static constexpr size_t dim = M::dim;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using MIdx = typename M::MIdx;
  using Direction = generic::Direction<dim>;

  IdxFaceMesh(IdxFace f, const M& m_) : idxface_(f), m(m_) {}
  IdxFaceMesh operator+(Direction d) const {
    return m(m.GetFace(idxface_, d.nci()));
  }
  IdxFaceMesh operator-(Direction d) const {
    return m(m.GetFace(idxface_, (-d).nci()));
  }
  friend std::ostream& operator<<(std::ostream& out, const IdxFaceMesh& f) {
    auto& index = f.m.GetIndexFaces();
    out << index.GetMIdx(f.idxface_) << index.GetDir(f.idxface_).GetLetter();
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
    } else if (sdir == "w") {
      dir = Dir(3);
    } else {
      in.setstate(std::ios_base::failbit);
      return in;
    }
    f.idxface_ = f.m.GetIndexFaces().GetIdx(w, dir);
    return in;
  }
  operator IdxFace() const {
    return idxface_;
  }
  explicit operator MIdx() const {
    return m.GetIndexFaces().GetMIdx(idxface_);
  }
  IntIdx operator[](size_t i) const {
    return MIdx(*this)[i];
  }
  Direction direction() const {
    return Direction(size_t(m.GetIndexFaces().GetDir(idxface_)));
  }
  auto cell(Side s) const {
    return m(m.GetCell(idxface_, s));
  }

  class LazyCenter {
   public:
    operator Vect() const {
      auto owner = GetOwner(this, &IdxFaceMesh::center);
      return owner->m.GetCenter(owner->idxface_);
    }
    Vect operator()() const {
      return Vect(*this);
    }
    Scal operator[](size_t i) const {
      return Vect(*this)[i];
    }

   private:
    friend IdxFaceMesh;
    LazyCenter(const LazyCenter&) = default;
    LazyCenter& operator=(const LazyCenter&) = default;
  };
  const LazyCenter center{};

  class LazyCellCm {
   public:
    operator IdxCell() const {
      auto owner = GetOwner(this, &IdxFaceMesh::cm);
      return owner->m.GetCell(owner->idxface_, 0);
    }
    operator IdxCellMesh<M>() const {
      auto owner = GetOwner(this, &IdxFaceMesh::cm);
      return owner->m(owner->m.GetCell(owner->idxface_, 0));
    }
    IdxCellMesh<M> operator()() const {
      return IdxCellMesh<M>(*this);
    }

   private:
    friend IdxFaceMesh;
    LazyCellCm(const LazyCellCm&) = default;
    LazyCellCm& operator=(const LazyCellCm&) = default;
  };
  const LazyCellCm cm{};

  class LazyCellCp {
   public:
    operator IdxCell() const {
      auto owner = GetOwner(this, &IdxFaceMesh::cp);
      return owner->m.GetCell(owner->idxface_, 1);
    }
    operator IdxCellMesh<M>() const {
      auto owner = GetOwner(this, &IdxFaceMesh::cp);
      return owner->m(owner->m.GetCell(owner->idxface_, 1));
    }
    IdxCellMesh<M> operator()() const {
      return IdxCellMesh<M>(*this);
    }

   private:
    friend IdxFaceMesh;
    LazyCellCp(const LazyCellCp&) = default;
    LazyCellCp& operator=(const LazyCellCp&) = default;
  };
  const LazyCellCp cp{};

  class LazyArea {
   public:
    operator Scal() const {
      auto owner = GetOwner(this, &IdxFaceMesh::area);
      return owner->m.GetArea(owner->idxface_);
    }
    Scal operator()() const {
      return Scal(*this);
    }

   private:
    friend IdxFaceMesh;
    LazyArea(const LazyArea&) = default;
    LazyArea& operator=(const LazyArea&) = default;
  };
  const LazyArea area{};

  class LazySurface {
   public:
    operator Vect() const {
      auto owner = GetOwner(this, &IdxFaceMesh::surface);
      return owner->m.GetSurface(owner->idxface_);
    }
    Vect operator()() const {
      return Vect(*this);
    }
    Scal operator[](size_t i) const {
      return Vect(*this)[i];
    }

   private:
    friend IdxFaceMesh;
    LazySurface(const LazySurface&) = default;
    LazySurface& operator=(const LazySurface&) = default;
  };
  const LazySurface surface{};

 private:
  IdxFace idxface_;
  const M& m;
};
