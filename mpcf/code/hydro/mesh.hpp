/*
 * mesh.hpp
 *
 *  Created on: Jan 25, 2016
 *      Author: Petr Karnakov (petr.karnakov@gmail.com)
 */

#pragma once

#include "vect.hpp"
#include <vector>
#include <array>
#include <cstddef>
#include <map>
#include <cassert>
#include <iostream>

namespace geom {

using IntIdx = std::ptrdiff_t;

template <size_t dim>
using MIdxGeneral = Vect<IntIdx, dim>;

template <int dummy>
class IdxGeneric {
  static constexpr int dummy_ = dummy; // make distinct classes
  size_t raw_;
  static constexpr size_t kNone = -1;
 public:
  IdxGeneric() {}
  explicit IdxGeneric(size_t raw)
      : raw_(raw)
  {}
  inline size_t GetRaw() const {
    return raw_;
  }
  inline void AddRaw(IntIdx add) {
    raw_ += add;
  }
  inline bool operator==(IdxGeneric other) const {
    return raw_ == other.raw_;
  }
  inline bool operator!=(IdxGeneric other) const {
    return !(*this == other);
  }
  inline static IdxGeneric None() {
    return IdxGeneric(-1);
  }
  inline bool IsNone() const {
    return *this == None();
  }
};

using IdxCell = IdxGeneric<0>;
using IdxFace = IdxGeneric<1>;
using IdxNode = IdxGeneric<2>;


template <class IdxType>
class Range {
  size_t pos_begin_, pos_end_;

 public:
  class iterator {
    size_t pos_;
   public:
    explicit iterator(size_t pos)
        : pos_(pos)
    {}
    iterator& operator++() {
      ++pos_;
      return *this;
    }
    iterator& operator--() {
      --pos_;
      return *this;
    }
    bool operator==(const iterator& other) const {
      return pos_ == other.pos_;
    }
    bool operator!=(const iterator& other) const {
      return !(*this == other);
    }
    IdxType operator*() const {
      return IdxType(pos_);
    }
  };

  Range()
      : pos_begin_(0)
      , pos_end_(0)
  {}
  Range(size_t pos_begin, size_t pos_end)
      : pos_begin_(pos_begin)
      , pos_end_(pos_end)
  {}
  iterator begin() const {
    return iterator(pos_begin_);
  }
  iterator end() const {
    return iterator(pos_end_);
  }
  size_t size() const {
    return pos_end_ - pos_begin_;
  }
};


template <class T, class Idx>
class FieldGeneric {
  using Vector = std::vector<T>;
  std::vector<T> data_;
  template <class OtherT, class OtherIdx>
  friend class FieldGeneric;

 public:
  using IdxType = Idx;
  using ValueType = T;
  FieldGeneric() {}
  template <class U>
  FieldGeneric(const FieldGeneric<U, Idx>& other)
      : data_(other.data_.begin(), other.data_.end()) {}
  explicit FieldGeneric(const Range<Idx>& range)
      : data_(range.size())
  {}
  FieldGeneric(const Range<Idx>& range, const T& value)
      : data_(range.size(), value)
  {}
  size_t size() const {
    return data_.size();
  }
  void Reinit(const Range<Idx>& range) {
    data_.resize(range.size());
  }
  void Reinit(const Range<Idx>& range, const T& value) {
    data_.assign(range.size(), value);
  }
  Range<IdxType> GetRange() const {
    return Range<IdxType>(0, size());
  }
  void resize(size_t size) {
    data_.resize(size);
  }
  bool empty() const {
    return data_.empty();
  }
  typename Vector::pointer data() {
    return data_.data();
  }
  typename Vector::const_pointer data() const {
    return data_.data();
  }
  void push_back(const T& value) {
    data_.push_back(value);
  }
  typename Vector::reference operator[](const Idx& idx) {
#ifdef __RANGE_CHECK
    assert(idx.GetRaw() >=0 && idx.GetRaw() < data_.size());
#endif
    return data_[idx.GetRaw()];
  }
  typename Vector::const_reference operator[](const Idx& idx) const {
#ifdef __RANGE_CHECK
    assert(idx.GetRaw() >=0 && idx.GetRaw() < data_.size());
#endif
    return data_[idx.GetRaw()];
  }
};

template <class T, class Idx>
FieldGeneric<typename T::value_type, Idx>
GetComponent(const FieldGeneric<T, Idx>& f_vect, size_t n) {
  FieldGeneric<typename T::value_type, Idx> f_scal(f_vect.GetRange());
  for (auto idx : f_vect.GetRange()) {
    f_scal[idx] = f_vect[idx][n];
  }
  return f_scal;
}

template <class T, class Idx>
void SetComponent(
    FieldGeneric<T, Idx>& f_vect,
    size_t n, const FieldGeneric<typename T::value_type, Idx>& f_scal) {
  for (auto idx : f_vect.GetRange()) {
    f_vect[idx][n] = f_scal[idx];
  }
}

template <class T, class Idx>
std::ostream& operator<<(std::ostream& out, const FieldGeneric<T, Idx>& field) {
  for (auto idx : field.GetRange()) {
    out << idx.GetRaw() << " " << field[idx] << "\n";
  }
  return out;
}

template <class T>
using FieldCell = FieldGeneric<T, IdxCell>;

template <class T>
using FieldFace = FieldGeneric<T, IdxFace>;

template <class T>
using FieldNode = FieldGeneric<T, IdxNode>;

template <class T, class Idx>
class MapGeneric {
  using Map = std::map<size_t, T>;
  Map data_;

 public:
  using IdxType = Idx;
  MapGeneric() {}
  explicit MapGeneric(const FieldGeneric<T, Idx>& field) {
    for (size_t i = 0; i < field.size(); ++i) {
      data_[i] = field[Idx(i)];
    }
  }
  size_t size() const {
    return data_.size();
  }
  T& operator[](const Idx& idx) {
    return data_[idx.GetRaw()];
  }
  const T& operator[](const Idx& idx) const {
    return data_[idx.GetRaw()];
  }
  T* find(const Idx& idx) {
    auto it = data_.find(idx.GetRaw());
    if (it != data_.end()) {
      return &it->second;
    }
    return nullptr;
  }
  const T * find(const Idx& idx) const {
    auto it = data_.find(idx.GetRaw());
    if (it != data_.end()) {
      return &it->second;
    }
    return nullptr;
  }
  void erase(const Idx& idx) {
    data_.erase(idx.GetRaw());
  }

  class iterator {
    class Proxy {
      typename Map::iterator it_;
     public:
      explicit Proxy(const typename Map::iterator& it)
          : it_(it)
      {}
      Idx GetIdx() const {
        return Idx(it_->first);
      }
      T& GetValue() const {
        return it_->second;
      }
    };

    typename Map::iterator it_;
    Proxy proxy_;

   public:
    explicit iterator(const typename Map::iterator& it)
        : it_(it)
        , proxy_(it_)
    {}
    iterator& operator++() {
      ++it_;
      proxy_ = Proxy(it_);
      return *this;
    }
    iterator& operator--() {
      --it_;
      proxy_ = Proxy(it_);
      return *this;
    }
    bool operator==(const iterator& other) const {
      return it_ == other.it_;
    }
    bool operator!=(const iterator& other) const {
      return !(*this == other);
    }
    const Proxy& operator*() {
      return proxy_;
    }
    Proxy const * operator->() {
      return &proxy_;
    }
  };
  iterator begin() {
    return iterator(data_.begin());
  }
  iterator end() {
    return iterator(data_.end());
  }
  class const_iterator {
    class Proxy {
      typename Map::const_iterator it_;
     public:
      explicit Proxy(const typename Map::const_iterator& it)
          : it_(it)
      {}
      Idx GetIdx() const {
        return Idx(it_->first);
      }
      const T& GetValue() const {
        return it_->second;
      }
    };

    typename Map::const_iterator it_;
    Proxy proxy_;

   public:
    explicit const_iterator(const typename Map::const_iterator& it)
        : it_(it)
        , proxy_(it_)
    {}
    const_iterator& operator++() {
      ++it_;
      proxy_ = Proxy(it_);
      return *this;
    }
    const_iterator& operator--() {
      --it_;
      proxy_ = Proxy(it_);
      return *this;
    }
    bool operator==(const const_iterator& other) const {
      return it_ == other.it_;
    }
    bool operator!=(const const_iterator& other) const {
      return !(*this == other);
    }
    const Proxy& operator*() {
      return proxy_;
    }
    Proxy const * operator->() {
      return &proxy_;
    }
  };
  const_iterator cbegin() const {
    return const_iterator(data_.cbegin());
  }
  const_iterator cend() const {
    return const_iterator(data_.cend());
  }
};

template <class T>
using MapCell = MapGeneric<T, IdxCell>;

template <class T>
using MapFace = MapGeneric<T, IdxFace>;

template <class T>
using MapNode = MapGeneric<T, IdxNode>;

// TODO: Neighbour faces iterator introducing (cell, face) pairs
template <class ScalArg, size_t dim, class K /*: Kernel*/>
class MeshGeneric {
 public:
  using Scal = ScalArg;
  using Vect = geom::Vect<Scal, dim>;
  MeshGeneric(K& kern) : kern_(kern) {}
  virtual ~MeshGeneric() {}
  virtual Vect GetCenter(IdxCell) const = 0;
  virtual Vect GetCenter(IdxFace) const = 0;
  virtual Vect GetSurface(IdxFace) const = 0;
  virtual Vect GetNode(IdxNode) const = 0;
  virtual Scal GetVolume(IdxCell) const = 0;
  virtual Scal GetArea(IdxFace) const = 0;
  virtual size_t GetNumCells() const = 0;
  virtual size_t GetNumFaces() const = 0;
  virtual size_t GetNumNodes() const = 0;
  virtual IdxCell GetNeighbourCell(IdxCell, size_t) const = 0;
  virtual IdxFace GetNeighbourFace(IdxCell, size_t) const = 0;
  virtual Scal GetOutwardFactor(IdxCell, size_t) const = 0;
  virtual Vect GetOutwardSurface(IdxCell, size_t) const = 0;
  virtual IdxNode GetNeighbourNode(IdxCell, size_t) const = 0;
  virtual IdxCell GetNeighbourCell(IdxFace, size_t) const = 0;
  virtual IdxNode GetNeighbourNode(IdxFace, size_t) const = 0;
  // For any cell, indices of its neighbour cells and faces are the same.
  // Thus, there's no GetNumNeighbourCells(IdxCell) function.
  virtual size_t GetNumNeighbourFaces(IdxCell) const = 0;
  virtual size_t GetNumNeighbourNodes(IdxCell) const = 0;
  virtual size_t GetNumNeighbourCells(IdxFace) const = 0;
  virtual size_t GetNumNeighbourNodes(IdxFace) const = 0;
  operator Range<IdxCell>() const {
    return Range<IdxCell>(0, GetNumCells());
  }
  operator Range<IdxFace>() const {
    return Range<IdxFace>(0, GetNumFaces());
  }
  operator Range<IdxNode>() const {
    return Range<IdxNode>(0, GetNumNodes());
  }
  Range<IdxCell> Cells() const {
    return Range<IdxCell>(*this);
  }
  Range<IdxFace> Faces() const {
    return Range<IdxFace>(*this);
  }
  Range<IdxNode> Nodes() const {
    return Range<IdxNode>(*this);
  }
  virtual bool IsInner(IdxCell) const = 0;
  virtual bool IsInner(IdxFace) const = 0;
  virtual bool IsInside(IdxCell, Vect) const = 0;
  virtual IdxCell FindNearestCell(Vect point) const {
    IdxCell idxcell_nearest = IdxCell(0);
    for (auto idxcell : Cells()) {
      if (point.dist(GetCenter(idxcell)) <
          point.dist(GetCenter(idxcell_nearest))) {
        idxcell_nearest = idxcell;
      }
    }
    return idxcell_nearest;
  }
  virtual size_t GetValidNeighbourCellId(IdxFace idxface) const {
    size_t id = 0;
    while (id < GetNumNeighbourCells(idxface) &&
        GetNeighbourCell(idxface, id).IsNone()) {
      ++id;
    }
    return id;
  }
  virtual Vect GetNormal(IdxFace idxface) const {
    return GetSurface(idxface) / GetArea(idxface);
  }

  using Sem = typename K::Sem;
  Sem GetSem(std::string name="") {
    return kern_.GetSem(name);
  }

 protected:
  K& kern_;
};


template <class IdxType, size_t dim>
class BlockGeneric {
  using MIdx = geom::MIdxGeneral<dim>;
  MIdx begin_, size_, end_;

 public:
  class iterator {
    const BlockGeneric* owner_;
    MIdx midx_;

   public:
    explicit iterator(const BlockGeneric* owner, MIdx midx)
        : owner_(owner)
        , midx_(midx)
    {}
    iterator& operator++() {
      for (size_t n = 0; n < dim; ++n) {
        ++midx_[n];
        if (midx_[n] == owner_->end_[n] && n < dim - 1) {
          midx_[n] = owner_->begin_[n];
        } else {
          break;
        }
      }
      return *this;
    }
    iterator& operator--() {
      for (size_t n = 0; n < dim; ++n) {
        if (midx_[n] == owner_->begin_[n]) {
          midx_[n] = owner_->end_[n] - 1;
        } else {
          --midx_[n];
          break;
        }
      }
      return *this;
    }
    bool operator==(const iterator& other) const {
      return midx_ == other.midx_;
    }
    bool operator!=(const iterator& other) const {
      return !(*this == other);
    }
    MIdx operator*() const {
      return midx_;
    }
  };

  BlockGeneric()
      : begin_(MIdx::kZero), size_(MIdx::kZero), end_(MIdx::kZero)
  {}
  BlockGeneric(MIdx size)
      : begin_(MIdx::kZero), size_(size), end_(begin_ + size_)
  {}
  BlockGeneric(MIdx begin, MIdx size)
      : begin_(begin), size_(size), end_(begin_ + size_)
  {}
  MIdx GetDimensions() const {
    return size_;
  }
  MIdx GetBegin() const {
    return begin_;
  }
  MIdx GetEnd() const {
    return end_;
  }
  size_t size() const {
    size_t res = 1;
    for (size_t i = 0; i < dim; ++i) {
      res *= size_[i];
    }
    return res;
  }
  operator Range<IdxType>() const {
    return Range<IdxType>(0, size());
  }
/*  void Advance(IdxType& idx, Direction dir, int n) const {
    switch (dir) {
      case Direction::i: {
        idx.AddRaw(n);
        break;
      }
      case Direction::j: {
        idx.AddRaw(n*size_[0]);
        break;
      }
      case Direction::k: {
        idx.AddRaw(n*size_[0]*size_[1]);
        break;
      }
    }
  }
  void Advance(IdxType& idx, MIdx midx) const {
    // TODO: unwrap multiple Advance()
    Advance(idx, Direction::i, midx[0]);
    Advance(idx, Direction::j, midx[1]);
    Advance(idx, Direction::k, midx[2]);
  }
  IdxType GetAdvanced(IdxType idx, MIdx midx) const {
    Advance(idx, midx);
    return idx;
  }
  void Inc(IdxType& idx, Direction dir) const {
    Advance(idx, dir, 1);
  }
  IdxType GetInc(IdxType idx, Direction dir) const {
    Inc(idx, dir);
    return idx;
  }
  void Dec(IdxType& idx, Direction dir) const {
    Advance(idx, dir, -1);
  }
  IdxType GetDec(IdxType idx, Direction dir) const {
    Dec(idx, dir);
    return idx;
  }*/
  IdxType GetIdx(MIdx midx) const {
    midx -= begin_;
    size_t res = 0;
    for (size_t i = dim; i != 0; ) {
      --i;
      res *= size_[i];
      res += midx[i];
    }
    return IdxType(res);
  }
  MIdx GetMIdx(IdxType idx) const {
    MIdx midx;
    size_t raw = idx.GetRaw();
    for (size_t i = 0; i < dim; ++i) {
      midx[i] = raw % size_[i];
      raw /= size_[i];
    }
    return midx;
  }
  bool IsInside(MIdx midx) const {
    return begin_ <= midx && midx < end_;
  }
  void Cut(MIdx& midx) const {
    for (size_t i = 0; i < midx.size(); ++i) {
      midx[i] = std::max(
          begin_[i],
          std::min(
              end_[i] - 1,
              midx[i]));
    }
  }
  iterator begin() const {
    return iterator(this, begin_);
  }
  iterator end() const {
    MIdx midx = begin_;
    midx[dim - 1] = end_[dim - 1];
    return iterator(this, midx);
  }
};

template <size_t dim>
using BlockCells = BlockGeneric<IdxCell, dim>;

template <size_t dim>
using BlockNodes = BlockGeneric<IdxNode, dim>;


template <size_t dim>
class Direction {
  using MIdx = Vect<IntIdx, dim>;
  size_t diridx;
 public:
  Direction() {}
  explicit Direction(size_t i)
    : diridx(i) {}
  char GetLetter() const {
    return std::string("xyz")[diridx];
  }
  operator size_t() const {
    return diridx;
  }
  operator MIdx() const {
    MIdx res = MIdx::kZero;
    ++res[diridx];
    return res;
  }
  bool operator==(const Direction& other) const {
    return diridx == other.diridx;
  }
  bool operator!=(const Direction& other) const {
    return !((*this) == other);
  }
  bool operator<(const Direction& other) const {
    return diridx < other.diridx;
  }
  static const Direction i;
  static const Direction j;
  static const Direction k;
};

template <size_t dim>
const Direction<dim> Direction<dim>::i(0);
template <size_t dim>
const Direction<dim> Direction<dim>::j(1);
template <size_t dim>
const Direction<dim> Direction<dim>::k(2);

template <size_t dim>
class BlockFaces {
  using MIdx = MIdxGeneral<dim>;
  using Direction = geom::Direction<dim>;
  MIdx block_cells_size_;
  size_t GetNumFaces(Direction dir) const {
    MIdx bcs = block_cells_size_;
    ++bcs[dir];
    size_t res = 1;
    for (size_t i = 0; i < dim; ++i) {
      res *= bcs[i];
    }
    return res;
  }
  size_t GetNumFaces(size_t i) const {
    return GetNumFaces(Direction(i));
  }
  size_t GetFlat(MIdx midx, Direction dir) const {
    MIdx bcs = block_cells_size_;
    ++bcs[dir];
    size_t res = 0;
    for (size_t i = dim; i != 0; ) {
      --i;
      res *= bcs[i];
      res += midx[i];
    }
    return res;
  }
  MIdx GetMIdxFromOffset(size_t raw, Direction dir) const {
    MIdx bcs = block_cells_size_;
    ++bcs[dir];
    MIdx midx;
    for (size_t i = 0; i < dim; ++i) {
      midx[i] = raw % bcs[i];
      raw /= bcs[i];
    }
    return midx;
  }
  std::pair<MIdx, Direction> GetMIdxDirection(IdxFace idxface) const {
    size_t raw = idxface.GetRaw();
    size_t diridx = 0;
    while (raw >= GetNumFaces(diridx) && diridx < dim) {
      raw -= GetNumFaces(diridx);
      ++diridx;
    }
    Direction dir(diridx);
    return {GetMIdxFromOffset(raw, dir), dir};
  }

 public:
  BlockFaces()
      : block_cells_size_(MIdx::kZero)
  {}
  BlockFaces(MIdx block_cells_size)
      : block_cells_size_(block_cells_size)
  {}
  size_t size() const {
    size_t res = 0;
    for (size_t i = 0; i < dim; ++i) {
      res += GetNumFaces(i);
    }
    return res;
  }
  operator Range<IdxFace>() const {
    return Range<IdxFace>(0, size());
  }
  IdxFace GetIdx(MIdx midx, Direction dir) const {
    size_t raw = 0;
    for (size_t i = 0; i < dir; ++i) {
      raw += GetNumFaces(i);
    }
    raw += GetFlat(midx, dir);
    return IdxFace(raw);
  }
  MIdx GetMIdx(IdxFace idx) const {
    return GetMIdxDirection(idx).first;
  }
  Direction GetDirection(IdxFace idx) const {
    return GetMIdxDirection(idx).second;
  }
};

template <class Mesh>
void InitUniformMesh(Mesh& mesh,
                     const Rect<typename Mesh::Vect>& domain,
                     typename Mesh::MIdx mesh_size) {
  using Vect = typename Mesh::Vect;
  using MIdx = typename Mesh::MIdx;
  typename Mesh::BlockNodes b_nodes(mesh_size + MIdx(1));
  FieldNode<Vect> fn_node(b_nodes);
  for (auto midx : b_nodes) {
    IdxNode idxnode = b_nodes.GetIdx(midx);
    Vect unit = Vect(midx) / Vect(mesh_size);
    fn_node[idxnode] = domain.lb + unit * domain.GetDimensions();
  }
  mesh = Mesh(b_nodes, fn_node);
}


template <class Mesh>
class SearchMesh {
  const Mesh& mesh;
  using Scal = typename Mesh::Scal;
  using Vect = typename Mesh::Vect;
  static constexpr size_t dim = Mesh::dim;
  using IdxSeg = geom::IdxGeneric<20160212>; // Cartesian segment idx
  using BlockSeg = geom::BlockGeneric<IdxSeg, dim>;
  using Rect = geom::Rect<Vect>;
  using MIdx = geom::MIdxGeneral<dim>;
  template <class T>
  using FieldSeg = geom::FieldGeneric<T, IdxSeg>;
  Rect bound_;
  BlockSeg block_;
  FieldSeg<std::vector<IdxCell>::iterator> fs_intersect_begin_;
  FieldSeg<std::vector<IdxCell>::iterator> fs_center_begin_;
  std::vector<IdxCell> content_intersect_;
  std::vector<IdxCell> content_center_;

  Rect GetBoundingRect(const Mesh& mesh) const {
    Rect res;
    res.lb = res.rt = mesh.GetNode(*mesh.Nodes().begin());
    for (auto idxnode : mesh.Nodes()) {
      const Vect node = mesh.GetNode(idxnode);
      for (size_t i = 0; i < dim; ++i) {
        res.lb[i] = std::min(res.lb[i], node[i]);
        res.rt[i] = std::max(res.rt[i], node[i]);
      }
    }
    return res;
  }
  Rect GetBoundingRect(IdxCell idxcell) const {
    Rect res;
    res.lb = res.rt = mesh.GetNode(mesh.GetNeighbourNode(idxcell, 0));
    for (size_t i = 0; i < mesh.GetNumNeighbourNodes(idxcell); ++i) {
      const Vect node = mesh.GetNode(mesh.GetNeighbourNode(idxcell, i));
      for (size_t i = 0; i < dim; ++i) {
        res.lb[i] = std::min(res.lb[i], node[i]);
        res.rt[i] = std::max(res.rt[i], node[i]);
      }
    }
    return res;
  }
  Rect GetRect(MIdx midx) const {
    const MIdx msize = block_.GetDimensions();
    const MIdx mlb = midx, mrt = midx + MIdx(1, 1);
    Rect res;
    for (size_t i = 0; i < dim; ++i) {
      res.lb[i] = (bound_.lb[i] * (block_.GetEnd()[i] - 1 - mlb[i]) +
          bound_.rt[i] * (mlb[i] - block_.GetBegin()[i])) / msize[i];
      res.rt[i] = (bound_.lb[i] * (block_.GetEnd()[i] - 1 - mrt[i]) +
          bound_.rt[i] * (mrt[i] - block_.GetBegin()[i])) / msize[i];
    }
    return res;
  }
  Rect GetRect(IdxSeg idxseg) const {
    return GetRect(block_.GetMIdx(idxseg));
  }
  MIdx GetMIdx(Vect vect) const {
    const Vect vsize = bound_.GetDimensions();

    MIdx res;
    for (size_t i = 0; i < dim; ++i) {
      res[i] = static_cast<int>((
          block_.GetBegin()[i] * (bound_.rt[i] - vect[i]) +
          (block_.GetEnd()[i] - 1) * (vect[i] - bound_.lb[i])) / vsize[i]);
    }

    block_.Cut(res);
    return res;
  }
  IdxSeg GetIdxSeg(Vect vect) const {
    return block_.GetIdx(GetMIdx(vect));
  }
  bool HaveCommonPoints(IdxCell, MIdx) const {
    return true;
  }
  void InitRefs(MIdx block_dimensions, size_t& num_intersect_seg_max,
                size_t& num_segs, size_t& num_center_seg_max) {
    block_ = BlockSeg(block_dimensions);

    // Find intersecting cells for each segment
    FieldSeg<std::vector<IdxCell>> intersect_chains_(block_);
    size_t intersect_total = 0;
    for (auto idxcell : mesh.Cells()) {
      Rect rect = GetBoundingRect(idxcell);
      MIdx mlb = GetMIdx(rect.lb);
      MIdx mrt = GetMIdx(rect.rt);
      for (auto midx : geom::BlockGeneric<size_t, dim>(
          mlb, mrt - mlb + MIdx(1))) {
        if (HaveCommonPoints(idxcell, midx)) {
          intersect_chains_[block_.GetIdx(midx)].push_back(idxcell);
          ++intersect_total;
        }
      }
    }

    num_intersect_seg_max = 0;
    fs_intersect_begin_.Reinit(geom::Range<IdxSeg>(0, block_.size() + 1));
    content_intersect_.resize(intersect_total);
    size_t intersect_offset = 0;
    for (size_t i = 0; i < block_.size(); ++i) {
      IdxSeg idxseg(i);
      fs_intersect_begin_[idxseg] =
          content_intersect_.begin() + intersect_offset;
      num_intersect_seg_max = std::max(num_intersect_seg_max,
                                    intersect_chains_[idxseg].size());
      for (auto idxcell : intersect_chains_[idxseg]) {
        content_intersect_[intersect_offset] = idxcell;
        ++intersect_offset;
      }
    }
    fs_intersect_begin_[IdxSeg(block_.size())] = content_intersect_.end();
    num_segs = block_.size();

    // Find cell centers lying in each segment
    FieldSeg<std::vector<IdxCell>> center_chains_(block_);
    for (auto idxcell : mesh.Cells()) {
      center_chains_[GetIdxSeg(mesh.GetCenter(idxcell))].push_back(idxcell);
    }

    num_center_seg_max = 0;
    fs_center_begin_.resize(block_.size() + 1);
    content_center_.resize(mesh.GetNumCells());
    size_t center_offset = 0;
    for (size_t i = 0; i < block_.size(); ++i) {
      IdxSeg idxseg(i);
      fs_center_begin_[idxseg] = content_center_.begin() + center_offset;
      num_center_seg_max = std::max(num_center_seg_max,
                                    center_chains_[idxseg].size());
      for (auto idxcell : center_chains_[idxseg]) {
        content_center_[center_offset] = idxcell;
        ++center_offset;
      }
    }
    fs_center_begin_[IdxSeg(block_.size())] = content_center_.end();
  }

 public:
  SearchMesh(const Mesh& mesh, size_t num_intersect_seg_limit,
             size_t num_segs_limit, size_t num_center_seg_limit)
      : mesh(mesh)
  {
    bound_ = GetBoundingRect(mesh);

    size_t num_intersect_seg_max;
    size_t num_center_seg_max;
    size_t num_segs;
    IntIdx block_dim = 1;
    do {
      InitRefs(MIdx(block_dim),
               num_intersect_seg_max, num_segs,
               num_center_seg_max);
      block_dim *= 2;
    } while (num_segs <= num_segs_limit &&
        (num_intersect_seg_max > num_intersect_seg_limit ||
            num_center_seg_max > num_center_seg_limit));

    std::cout << "Init SearchMesh with"
        << "\n num_intersect_seg_max = " << num_intersect_seg_max
        << "\n num_segs = " << num_segs
        << "\n num_center_seg_max = " << num_center_seg_max
        << std::endl;
  }
  IdxCell FindCell(Vect vect) const {
    IdxSeg idxseg = GetIdxSeg(vect);
    for (auto it = fs_intersect_begin_[idxseg];
        it != fs_intersect_begin_[IdxSeg(idxseg.GetRaw() + 1)]; ++it) {
      if (mesh.IsInside(*it, vect)) {
        return *it;
      }
    }
    return IdxCell::None();
  }
  std::vector<IdxCell> GetNearbyCells(Vect center, Scal radius) const {
    std::vector<IdxCell> res;
    Rect bound(center - Vect(radius), center + Vect(radius));
    MIdx mlb = GetMIdx(bound.lb);
    MIdx mrt = GetMIdx(bound.rt);
    for (auto midx : geom::BlockGeneric<size_t, dim>(
      mlb, mrt - mlb + MIdx(1))) {
      IdxSeg idxseg = block_.GetIdx(midx);
      for (auto it = fs_center_begin_[idxseg];
          it != fs_center_begin_[IdxSeg(idxseg.GetRaw() + 1)]; ++it) {
        if (center.dist(mesh.GetCenter(*it)) <= radius) {
          res.push_back(*it);
        }
      }
    }
    return res;
  }

};


} // namespace geom
