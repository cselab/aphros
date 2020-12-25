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
  Direction next(int shift = 1) const {
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
  size_t index_; // direction, [0, dim - 1]
  size_t forward_; // 0: backward, 1: forward
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
  using Direction = generic::Direction<Scal, dim>;

  IdxFaceMesh(IdxFace f, const M& m_) : idxface_(f), m(m_) {}
  IdxFaceMesh operator+(Direction d) const {
    return m(m.GetFace(idxface_, d.GetIndex() * 2 + d.GetForward()));
  }
  IdxFaceMesh operator-(Direction d) const {
    return m(m.GetFace(idxface_, d.GetIndex() * 2 + 1 - d.GetForward()));
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

 private:
  IdxFace idxface_;
  const M& m;
};
