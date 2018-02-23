/*
 * mesh.hpp *
 *  Created on: Jan 25, 2016
 *      Author: Petr Karnakov (petr.karnakov@gmail.com)
 */

#pragma once

#include <vector>
#include <array>
#include <cstddef>
#include <map>
#include <cassert>
#include <iostream>

#include "suspender.h"
#include "vect.hpp"


namespace geom {

using IntIdx = std::ptrdiff_t;


template <size_t dim>
using GMIdx = GVect<IntIdx, dim>;


template <int dummy>
class GIdx {
  static constexpr int dummy_ = dummy; // make distinct classes
  size_t raw_;
  static constexpr size_t kNone = -1;
 public:
  GIdx() {}
  explicit GIdx(size_t raw)
      : raw_(raw)
  {}
  inline size_t GetRaw() const {
    return raw_;
  }
  inline void AddRaw(IntIdx add) {
    raw_ += add;
  }
  inline bool operator==(GIdx other) const {
    return raw_ == other.raw_;
  }
  inline bool operator!=(GIdx other) const {
    return !(*this == other);
  }
  inline static GIdx None() {
    return GIdx(-1);
  }
  inline bool IsNone() const {
    return *this == None();
  }
};


template <class _Idx>
class GRange {
  size_t pos_begin_, pos_end_;

 public:
  using Idx = _Idx;
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
    Idx operator*() const {
      return Idx(pos_);
    }
  };

  GRange()
      : pos_begin_(0)
      , pos_end_(0)
  {}
  GRange(size_t pos_begin, size_t pos_end)
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


template <class _Value, class _Idx>
class GField {
 public:
  using Idx = _Idx;
  using Value = _Value;
  using Range = GRange<Idx>;
  using Cont = std::vector<Value>;

  GField() {}
  template <class U>
  GField(const GField<U, Idx>& o /*other*/)
      : data_(o.data_.begin(), o.data_.end()) {}
  explicit GField(const Range& range)
      : data_(range.size()) {}
  GField(const Range& range, const Value& value)
      : data_(range.size(), value)
  {}
  size_t size() const {
    return data_.size();
  }
  void Reinit(const Range& range) {
    data_.resize(range.size());
  }
  void Reinit(const Range& range, const Value& value) {
    data_.assign(range.size(), value);
  }
  Range GetRange() const {
    return Range(0, size());
  }
  void resize(size_t size) {
    data_.resize(size);
  }
  bool empty() const {
    return data_.empty();
  }
  typename Cont::pointer data() {
    return data_.data();
  }
  typename Cont::const_pointer data() const {
    return data_.data();
  }
  void push_back(const Value& value) {
    data_.push_back(value);
  }
  typename Cont::reference operator[](const Idx& idx) {
#ifdef __RANGE_CHECK
    assert(idx.GetRaw() >=0 && idx.GetRaw() < data_.size());
#endif
    return data_[idx.GetRaw()];
  }
  typename Cont::const_reference operator[](const Idx& idx) const {
#ifdef __RANGE_CHECK
    assert(idx.GetRaw() >=0 && idx.GetRaw() < data_.size());
#endif
    return data_[idx.GetRaw()];
  }

 private:

  template <class OValue, class OIdx>
  friend class GField;

  std::vector<Value> data_;
};




using IdxCell = GIdx<0>;
using IdxFace = GIdx<1>;
using IdxNode = GIdx<2>;


// Extract component of field with values GVect<T,dim>
template <class Vect, class Idx>
GField<typename Vect::value_type, Idx>
GetComponent(const GField<Vect, Idx>& fv, size_t n) {
  GField<typename Vect::value_type, Idx> fs(fv.GetRange());
  for (auto idx : fv.GetRange()) {
    fs[idx] = fv[idx][n];
  }
  return fs;
}

// Set component of field with values GVect<T,dim>
template <class Vect, class Idx>
void SetComponent(
    GField<Vect, Idx>& fv,
    size_t n, 
    const GField<typename Vect::value_type, Idx>& fs) {
  for (auto idx : fv.GetRange()) {
    fv[idx][n] = fs[idx];
  }
}

template <class T, class Idx>
std::ostream& operator<<(std::ostream& out, const GField<T, Idx>& field) {
  for (auto idx : field.GetRange()) {
    out << idx.GetRaw() << " " << field[idx] << "\n";
  }
  return out;
}

template <class T>
using FieldCell = GField<T, IdxCell>;

template <class T>
using FieldFace = GField<T, IdxFace>;

template <class T>
using FieldNode = GField<T, IdxNode>;

template <class T, class _Idx>
class GMap {
 public:
  using Idx = _Idx;
  using Value = T;
  using Cont = std::map<size_t, Value>;

  GMap() {}
  explicit GMap(const GField<T, Idx>& field) {
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
      typename Cont::iterator it_;
     public:
      explicit Proxy(const typename Cont::iterator& it)
          : it_(it)
      {}
      Idx GetIdx() const {
        return Idx(it_->first);
      }
      T& GetValue() const {
        return it_->second;
      }
    };

    typename Cont::iterator it_;
    Proxy proxy_;

   public:
    explicit iterator(const typename Cont::iterator& it)
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
      typename Cont::const_iterator it_;
     public:
      explicit Proxy(const typename Cont::const_iterator& it)
          : it_(it)
      {}
      Idx GetIdx() const {
        return Idx(it_->first);
      }
      const T& GetValue() const {
        return it_->second;
      }
    };

    typename Cont::const_iterator it_;
    Proxy proxy_;

   public:
    explicit const_iterator(const typename Cont::const_iterator& it)
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
 
 private:
  Cont data_;
};

template <class T>
using MapCell = GMap<T, IdxCell>;

template <class T>
using MapFace = GMap<T, IdxFace>;

template <class T>
using MapNode = GMap<T, IdxNode>;

template <class _Idx, size_t _dim>
class GBlock {
 public:
  using Idx = _Idx;
  using MIdx = GMIdx<_dim>;

  static constexpr size_t dim = _dim;

  class iterator {
    const GBlock* owner_;
    MIdx midx_;

   public:
    explicit iterator(const GBlock* owner, MIdx midx)
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

  GBlock()
      : begin_(MIdx::kZero), size_(MIdx::kZero), end_(MIdx::kZero)
  {}
  GBlock(MIdx size)
      : begin_(MIdx::kZero), size_(size), end_(begin_ + size_)
  {}
  GBlock(MIdx begin, MIdx size)
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
  operator GRange<Idx>() const {
    return GRange<Idx>(0, size());
  }
  Idx GetIdx(MIdx midx) const {
    midx -= begin_;
    size_t res = 0;
    for (size_t i = dim; i != 0; ) {
      --i;
      res *= size_[i];
      res += midx[i];
    }
    return Idx(res);
  }
  MIdx GetMIdx(Idx idx) const {
    MIdx midx;
    size_t raw = idx.GetRaw();
    for (size_t i = 0; i < dim; ++i) {
      midx[i] = raw % size_[i];
      raw /= size_[i];
    }
    return begin_ + midx;
  }
  bool IsInside(MIdx midx) const {
    return begin_ <= midx && midx < end_;
  }
  // TODO: rename to clip
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

 private:
  MIdx begin_, size_, end_;
};

template <size_t _dim>
class GDir {
 public:
  using MIdx = GVect<IntIdx, _dim>;
  static constexpr size_t dim = _dim;

  GDir() {}
  explicit GDir(size_t d)
    : d_(d) {}
  char GetLetter() const {
    return std::string("xyz")[d_];
  }
  operator size_t() const {
    return d_;
  }
  operator MIdx() const {
    MIdx r = MIdx::kZero;
    ++r[d_];
    return r;
  }
  bool operator==(const GDir& o) const {
    return d_ == o.d_;
  }
  bool operator!=(const GDir& o) const {
    return !((*this) == o);
  }
  bool operator<(const GDir& o) const {
    return d_ < o.d_;
  }
  static const GDir i;
  static const GDir j;
  static const GDir k;

 private:
  size_t d_;
};

template <size_t dim>
const GDir<dim> GDir<dim>::i(0);
template <size_t dim>
const GDir<dim> GDir<dim>::j(1);
template <size_t dim>
const GDir<dim> GDir<dim>::k(2);


template <size_t _dim>
class GBlock<IdxFace, _dim> {
 public:
  using Idx = IdxFace;
  using MIdx = GMIdx<_dim>;
  using Dir = GDir<_dim>;
  static constexpr size_t dim = _dim;

  GBlock()
      : b_(MIdx::kZero), cs_(MIdx::kZero)
  {}
  GBlock(MIdx block_cells_size)
      : b_(MIdx::kZero), cs_(block_cells_size)
  {}
  GBlock(MIdx begin, MIdx block_cells_size)
      : b_(begin), cs_(block_cells_size)
  {}
  size_t size() const {
    size_t res = 0;
    for (size_t i = 0; i < dim; ++i) {
      res += GetNumFaces(i);
    }
    return res;
  }
  operator GRange<Idx>() const {
    return GRange<Idx>(0, size());
  }
  Idx GetIdx(MIdx midx, Dir dir) const {
    size_t raw = 0;
    for (size_t i = 0; i < dir; ++i) {
      raw += GetNumFaces(i);
    }
    raw += GetFlat(midx, dir);
    return Idx(raw);
  }
  MIdx GetMIdx(Idx idx) const {
    return GetMIdxDir(idx).first;
  }
  Dir GetDir(Idx idx) const {
    return GetMIdxDir(idx).second;
  }

 private:
  MIdx b_;
  MIdx cs_; // cells size
  size_t GetNumFaces(Dir dir) const {
    MIdx bcs = cs_;
    ++bcs[dir];
    size_t res = 1;
    for (size_t i = 0; i < dim; ++i) {
      res *= bcs[i];
    }
    return res;
  }
  size_t GetNumFaces(size_t i) const {
    return GetNumFaces(Dir(i));
  }
  size_t GetFlat(MIdx midx, Dir dir) const {
    midx -= b_;
    MIdx bcs = cs_;
    ++bcs[dir];
    size_t res = 0;
    for (size_t i = dim; i != 0; ) {
      --i;
      res *= bcs[i];
      res += midx[i];
    }
    return res;
  }
  MIdx GetMIdxFromOffset(size_t raw, Dir dir) const {
    MIdx bcs = cs_;
    ++bcs[dir];
    MIdx midx;
    for (size_t i = 0; i < dim; ++i) {
      midx[i] = raw % bcs[i];
      raw /= bcs[i];
    }
    return b_ + midx;
  }
  std::pair<MIdx, Dir> GetMIdxDir(Idx idxface) const {
    size_t raw = idxface.GetRaw();
    size_t diridx = 0;
    while (raw >= GetNumFaces(diridx) && diridx < dim) {
      raw -= GetNumFaces(diridx);
      ++diridx;
    }
    Dir dir(diridx);
    return {GetMIdxFromOffset(raw, dir), dir};
  }
};

template <size_t dim>
using GBlockCells = GBlock<IdxCell, dim>;

template <size_t dim>
using GBlockNodes = GBlock<IdxNode, dim>;

template <size_t dim>
using GBlockFaces = GBlock<IdxFace, dim>;


template <class _Idx, int _dim>
class GRangeIn {
  using Idx = _Idx;
  using B = GBlock<Idx, _dim>; // block 
  using I = typename B::iterator; // block iterator
  const B& ba_; // block all
  const B& bi_; // block inner

 public:
  class iterator {
    const B& ba_; // block all
    I i_;  // block iterator over bi
   public:
    explicit iterator(const B& ba, const I& i)
        : ba_(ba), i_(i)
    {}
    iterator& operator++() {
      ++i_;
      return *this;
    }
    iterator& operator--() {
      --i_;
      return *this;
    }
    bool operator==(const iterator& o /*other*/) const {
      return i_ == o.i_;
    }
    bool operator!=(const iterator& o /*other*/) const {
      return !(*this == o);
    }
    Idx operator*() const {
      return ba_.GetIdx(*i_);
    }
  };

  GRangeIn(const B& ba /*block all*/, const B& bi /*block inner*/)
      : ba_(ba), bi_(bi)
  {}
  iterator begin() const {
    return iterator(ba_, bi_.begin());
  }
  iterator end() const {
    return iterator(ba_, bi_.end());
  }
};


template <int dim>
class GRangeIn<IdxFace, dim> {
  using B = GBlockFaces<dim>; // block 
  using R = GRange<IdxFace>;  // range
  using I = typename R::iterator; // range iterator
  const B& ba_; // block all
  const B& bi_; // block inner
  R ri_; // range inner

 public:
  class iterator {
    const B& ba_; // block all
    const B& bi_; // block inner
    I i_;         // range iterator over bi
   public:
    explicit iterator(const B& ba, const B& bi, const I& i)
        : ba_(ba), bi_(bi), i_(i)
    {}
    iterator& operator++() {
      ++i_;
      return *this;
    }
    iterator& operator--() {
      --i_;
      return *this;
    }
    bool operator==(const iterator& o /*other*/) const {
      return i_ == o.i_;
    }
    bool operator!=(const iterator& o /*other*/) const {
      return !(*this == o);
    }
    IdxFace operator*() const {
      auto x = bi_.GetMIdx(*i_); 
      auto d = bi_.GetDir(*i_); 
      return ba_.GetIdx(x, d);
    }
  };

  GRangeIn(const B& ba /*block all*/, const B& bi /*block inner*/)
      : ba_(ba), bi_(bi), ri_(bi_)
  {}
  iterator begin() const {
    return iterator(ba_, bi_, ri_.begin());
  }
  iterator end() const {
    return iterator(ba_, bi_, ri_.end());
  }
};

/*
template <class ScalArg, size_t dim>
class MeshGeneric {
 public:
  using Scal = ScalArg;
  using Vect = GVect<Scal, dim>;
  MeshGeneric() = default;
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
  virtual bool IsInner(IdxCell) const = 0;
  virtual bool IsInner(IdxFace) const = 0;
  virtual bool IsInside(IdxCell, Vect) const = 0;
};
*/

// TODO: Neighbour faces iterator introducing (cell, face) pairs
// TODO: consider computing some on-the-fly to reduce memory access
// TODO: separate assert() for range check and flag to disable on release
//
// Flag _generic added so that specializations could
// inherit the generic class without recursion.
template <class _Scal, size_t _dim>
class MeshStructured {
 public:
  static constexpr size_t dim = _dim;
  using Scal = _Scal;
  using Vect = GVect<Scal, dim>;
  using Dir = GDir<dim>;
  using MIdx = GVect<IntIdx, dim>;
  using BlockNodes = GBlockNodes<dim>;
  using BlockCells = GBlockCells<dim>;
  using BlockFaces = GBlockFaces<dim>;
  static constexpr size_t kCellNumNeighbourFaces = 2 * dim;
  static constexpr size_t kCellNumNeighbourNodes = std::pow(2, dim);
  static constexpr size_t kFaceNumNeighbourNodes = std::pow(2, dim - 1);
  static constexpr size_t kFaceNumNeighbourCells = 2;

 private:
  // b:Block, fc:FieldCell, ff:FieldFace, fn:FieldNode
  BlockNodes b_nodes_;
  BlockCells b_cells_;
  BlockFaces b_faces_;
  // inner 
  BlockNodes b_innodes_;
  BlockCells b_incells_;
  BlockFaces b_infaces_;

  FieldNode<Vect> fn_node_;

  FieldCell<Vect> fc_center_;
  FieldCell<Scal> fc_volume_;
  FieldFace<Vect> ff_center_;
  FieldFace<Scal> ff_area_;
  FieldFace<Vect> ff_surface_;
  std::array<IntIdx, kCellNumNeighbourFaces> cell_neighbour_cell_offset_;
  FieldCell<std::array<IdxFace, kCellNumNeighbourFaces>> fc_neighbour_face_;
  FieldCell<IdxNode> fc_neighbour_node_base_;
  std::array<IntIdx, kCellNumNeighbourNodes> cell_neighbour_node_offset_;
  FieldFace<Dir> ff_direction_;
  FieldFace<std::array<IdxCell, kFaceNumNeighbourCells>> ff_neighbour_cell_;
  FieldFace<std::array<IdxNode, kFaceNumNeighbourNodes>> ff_neighbour_node_;
  FieldFace<std::array<Vect, kFaceNumNeighbourCells>> ff_to_cell_;
  FieldCell<bool> fc_is_inner_;
  FieldFace<bool> ff_is_inner_;

 public:
  MeshStructured(const BlockNodes& b_nodes, 
                 const FieldNode<Vect>& fn_node, 
                 int hl /*halos*/);
  const BlockCells& GetBlockCells() const {
    return b_cells_;
  }
  const BlockFaces& GetBlockFaces() const {
    return b_faces_;
  }
  const BlockNodes& GetBlockNodes() const {
    return b_nodes_;
  }
  const BlockCells& GetInBlockCells() const {
    return b_incells_;
  }
  const BlockFaces& GetInBlockFaces() const {
    return b_infaces_;
  }
  const BlockNodes& GetInBlockNodes() const {
    return b_innodes_;
  }
  Vect GetCenter(IdxCell idx) const {
    return fc_center_[idx];
  }
  Vect GetCenter(IdxFace idx) const {
    return ff_center_[idx];
  }
  Vect GetSurface(IdxFace idx) const {
    return ff_surface_[idx];
  }
  Vect GetNode(IdxNode idx) const {
    return fn_node_[idx];
  }
  Scal GetVolume(IdxCell idx) const {
    return fc_volume_[idx];
  }
  Scal GetArea(IdxFace idx) const {
    return ff_area_[idx];
  }
  size_t GetNumCells() const {
    return b_cells_.size();
  }
  size_t GetNumFaces() const {
    return b_faces_.size();
  }
  size_t GetNumNodes() const {
    return b_nodes_.size();
  }
  IdxCell GetNeighbourCell(IdxCell idx, size_t n) const {
    // TODO: 3d specific
    assert(n < kCellNumNeighbourFaces);
    if (fc_is_inner_[idx]) {
      idx.AddRaw(cell_neighbour_cell_offset_[n]);
    } else {
      MIdx base = b_cells_.GetMIdx(idx);
      MIdx offset = (std::vector<MIdx>{
        MIdx{-1, 0, 0}, MIdx{1, 0, 0},
        MIdx{0, -1, 0}, MIdx{0, 1, 0},
        MIdx{0, 0, -1}, MIdx{0, 0, 1}})[n];
      if (b_cells_.IsInside(base + offset)) {
        idx.AddRaw(cell_neighbour_cell_offset_[n]);
      } else {
        idx = IdxCell::None();
      }
    }
    return idx;
  }
  IdxFace GetNeighbourFace(IdxCell idxcell, size_t n) const {
    assert(n < kCellNumNeighbourFaces);
    return fc_neighbour_face_[idxcell][n];
  }
  Scal GetOutwardFactor(IdxCell, size_t n) const {
    // TODO: <= 3d specific, maybe replace with odd/even convention
    assert(n < kCellNumNeighbourFaces);
    Scal factor;
    switch (n) {
      case 0:
      case 2:
      case 4:
        factor = -1.;
        break;
      default: // n == 1, 3, 5
        factor = 1.;
    }
    return factor;
  }
  Vect GetOutwardSurface(IdxCell idxcell, size_t n) const {
    assert(n < kCellNumNeighbourFaces);
    return GetSurface(GetNeighbourFace(idxcell, n)) *
        GetOutwardFactor(idxcell, n);
  }
  IdxNode GetNeighbourNode(IdxCell idxcell, size_t n) const {
    assert(n < kCellNumNeighbourNodes);
    IdxNode idxnode = fc_neighbour_node_base_[idxcell];
    idxnode.AddRaw(cell_neighbour_node_offset_[n]);
    return idxnode;
  }
  Dir GetDir(IdxFace idx) const {
    return ff_direction_[idx];
  }
  IdxCell GetNeighbourCell(IdxFace idx, size_t n) const {
    assert(n < kFaceNumNeighbourCells);
    return ff_neighbour_cell_[idx][n];
  }
  Vect GetVectToCell(IdxFace idx, size_t n) const {
    assert(n < kFaceNumNeighbourCells);
    return ff_to_cell_[idx][n];
  }
  IdxNode GetNeighbourNode(IdxFace idx, size_t n) const {
    assert(n < kFaceNumNeighbourNodes);
    return ff_neighbour_node_[idx][n];
  }
  size_t GetNumNeighbourFaces(IdxCell) const {
    return kCellNumNeighbourFaces;
  }
  size_t GetNumNeighbourNodes(IdxCell) const {
    return kCellNumNeighbourNodes;
  }
  size_t GetNumNeighbourCells(IdxFace) const {
    return kFaceNumNeighbourCells;
  }
  size_t GetNumNeighbourNodes(IdxFace) const {
    return kFaceNumNeighbourNodes;
  }
  // inner means not halo
  bool IsInner(IdxCell idxcell) const {
    return fc_is_inner_[idxcell];
  }
  bool IsInner(IdxFace idxface) const {
    return ff_is_inner_[idxface];
  }

  // Returns range of all indices
  template <class Idx>
  GRange<Idx> GetAll() const { 
    return GetAll((Idx*)0);
  }
  GRange<IdxCell> GetAll(IdxCell*) const { 
    return GRange<IdxCell>(0, GetNumCells());
  }
  GRange<IdxFace> GetAll(IdxFace*) const { 
    return GRange<IdxFace>(0, GetNumFaces());
  }
  GRange<IdxNode> GetAll(IdxNode*) const { 
    return GRange<IdxNode>(0, GetNumNodes());
  }
  GRange<IdxCell> AllCells() const {
    return GetAll<IdxCell>();
  }
  GRange<IdxFace> AllFaces() const {
    return GetAll<IdxFace>();
  }
  GRange<IdxNode> AllNodes() const {
    return GetAll<IdxNode>();
  }

  // Returns range of inner indices
  template <class Idx>
  GRangeIn<Idx, dim> Get() const { 
    return Get((Idx*)0);
  }
  GRangeIn<IdxCell, dim> Get(IdxCell*) const {
    return GRangeIn<IdxCell, dim>(GetBlockCells(), GetInBlockCells());
  }
  GRangeIn<IdxFace, dim> Get(IdxFace*) const {
    return GRangeIn<IdxFace, dim>(GetBlockFaces(), GetInBlockFaces());
  }
  GRangeIn<IdxNode, dim> Get(IdxNode*) const {
    return GRangeIn<IdxNode, dim>(GetBlockNodes(), GetInBlockNodes());
  }
  GRangeIn<IdxCell, dim> Cells() const {
    return Get<IdxCell>();
  }
  GRangeIn<IdxFace, dim> Faces() const {
    return Get<IdxFace>();
  }
  GRangeIn<IdxNode, dim> Nodes() const {
    return Get<IdxNode>();
  }

  bool IsInside(IdxCell idxcell, Vect vect) const {
    for (size_t i = 0; i < GetNumNeighbourFaces(idxcell); ++i) {
      IdxFace idxface = GetNeighbourFace(idxcell, i);
      if (GetOutwardSurface(idxcell, i).dot(vect - GetCenter(idxface)) > 0) {
        return false;
      }
    }
    return true;
  }
  IdxCell FindNearestCell(Vect point) const {
    IdxCell idxcell_nearest = IdxCell(0);
    for (auto idxcell : AllCells()) {
      if (point.dist(GetCenter(idxcell)) <
          point.dist(GetCenter(idxcell_nearest))) {
        idxcell_nearest = idxcell;
      }
    }
    return idxcell_nearest;
  }
  size_t GetValidNeighbourCellId(IdxFace idxface) const {
    size_t id = 0;
    while (id < GetNumNeighbourCells(idxface) &&
        GetNeighbourCell(idxface, id).IsNone()) {
      ++id;
    }
    return id;
  }
  Vect GetNormal(IdxFace idxface) const {
    return GetSurface(idxface) / GetArea(idxface);
  }

 private:
   Vect CalcCenter(IdxCell idxcell) const {
     Vect res = Vect::kZero;
     for (size_t i = 0; i < GetNumNeighbourNodes(idxcell); ++i) {
       res += GetNode(GetNeighbourNode(idxcell, i));
     }
     return res / Scal(GetNumNeighbourNodes(idxcell));
   }
   // TODO: Face center using weighted sum
   Vect CalcCenter(IdxFace idxface) const {
     Vect res = Vect::kZero;
     for (size_t i = 0; i < GetNumNeighbourNodes(idxface); ++i) {
       res += GetNode(GetNeighbourNode(idxface, i));
     }
     return res / Scal(GetNumNeighbourNodes(idxface));
   }
   static Vect CalcTriangleSurface(Vect a, Vect b, Vect c) {
     return (b - a).cross(c - a) * 0.5;
   }
   static Vect CalcTriangleCenter(Vect a, Vect b, Vect c) {
     return (a + b + c) / 3.;
   }
   Vect CalcSurface(IdxFace idxface) const {
     Vect res = Vect::kZero;

     for (size_t i = 2; i < GetNumNeighbourNodes(idxface); ++i) {
       Vect a = GetNode(GetNeighbourNode(idxface, 0));
       Vect b = GetNode(GetNeighbourNode(idxface, i - 1));
       Vect c = GetNode(GetNeighbourNode(idxface, i));
       res += CalcTriangleSurface(a, b, c);
     }

     return res;
   }
   Scal CalcArea(IdxFace idxface) const {
     return CalcSurface(idxface).norm();
   }
   Scal CalcVolume(IdxCell idxcell) const {
     Scal res = 0.;
     for (size_t i = 0; i < GetNumNeighbourFaces(idxcell); ++i) {
       IdxFace idxface = GetNeighbourFace(idxcell, i);
       res += GetCenter(idxface)[0] * GetOutwardSurface(idxcell, i)[0];
     }
     return res;
   }

  // TODO: move to separate class: Sem, LS, Comm, Reduce, Solve
  // BEGIN DISTR
 public:
  using Sem = Suspender::Sem;
  Sem GetSem(std::string name="") {
    return susp_.GetSem(name);
  }
  bool Pending() const {
    return susp_.Pending();
  }

  std::string GetCurName() const {
    return susp_.GetCurName();
  }
  size_t GetDepth() const {
    return susp_.GetDepth();
  }
  struct LS { // linear system ax=b
    std::vector<MIdx> st; // stencil
    std::vector<Scal>* a;
    std::vector<Scal>* b; 
    std::vector<Scal>* x;
  };

  void Comm(FieldCell<Scal>* u) {
    vcm_.push_back(u);
  }
  void Reduce(Scal* u) {
    vrd_.push_back(u);
  }
  void Solve(const LS& ls) {
    vls_.push_back(ls);
  }
  const std::vector<FieldCell<Scal>*>& GetComm() const {
    return vcm_;
  }
  void ClearComm() {
    vcm_.clear();
  }
  const std::vector<Scal*>& GetReduce() const {
    return vrd_;
  }
  void ClearReduce() {
    vrd_.clear();
  }
  const std::vector<LS>& GetSolve() const {
    return vls_;
  }
  void ClearSolve() {
    vls_.clear();
  }

 private:
  Suspender susp_;
  std::vector<FieldCell<Scal>*> vcm_; // fields for [c]o[m]munication
  std::vector<Scal*> vrd_; // scalars for reduction
  std::vector<LS> vls_; // linear system
  // END DISTR
};


template <class _Scal, size_t _dim>
MeshStructured<_Scal, _dim>::MeshStructured(
    const BlockNodes& b_nodes,
    const geom::FieldNode<Vect>& fn_node,
    int hl /*halos*/)
    : b_nodes_(b_nodes)
    , fn_node_(fn_node)
    , b_cells_(b_nodes_.GetBegin(), b_nodes_.GetDimensions() - MIdx(1))
    , b_faces_(b_nodes_.GetBegin(), b_cells_.GetDimensions())
    , b_incells_(b_cells_.GetBegin() + MIdx(hl), 
                 b_cells_.GetDimensions() - MIdx(2*hl))
    , b_infaces_(b_cells_.GetBegin() + MIdx(hl), 
                 b_cells_.GetDimensions() - MIdx(2*hl))
    , b_innodes_(b_nodes_.GetBegin() + MIdx(hl), 
                 b_nodes_.GetDimensions() - MIdx(2*hl))
    , fc_center_(b_cells_)
    , fc_volume_(b_cells_)
    , ff_center_(b_faces_)
    , ff_area_(b_faces_)
    , ff_surface_(b_faces_)
    , fc_neighbour_face_(b_cells_)
    , fc_neighbour_node_base_(b_cells_)
    , ff_direction_(b_faces_)
    , ff_neighbour_cell_(b_faces_)
    , ff_neighbour_node_(b_faces_)
    , ff_to_cell_(b_faces_)
    , fc_is_inner_(b_cells_, false)
    , ff_is_inner_(b_faces_, false)
{
  static_assert(dim == 3, "Not implemented for dim != 3");
  MIdx mb = b_cells_.GetBegin(), me = b_cells_.GetEnd();

  // Base for cell neighbours
  for (auto midx : b_cells_) {
    IdxCell idxcell(b_cells_.GetIdx(midx));
    fc_neighbour_node_base_[idxcell] = b_nodes_.GetIdx(midx);

    auto nface = [this, midx](IntIdx i, IntIdx j, IntIdx k, Dir dir) {
      return b_faces_.GetIdx(midx + MIdx(i, j, k), dir);
    };
    fc_neighbour_face_[idxcell] = {{
        nface(0, 0, 0, Dir::i),
        nface(1, 0, 0, Dir::i),
        nface(0, 0, 0, Dir::j),
        nface(0, 1, 0, Dir::j),
        nface(0, 0, 0, Dir::k),
        nface(0, 0, 1, Dir::k)}};
  }

  // Mark inner cells
  for (auto i : this->Cells()) {
    fc_is_inner_[i] = true;
  }

  // Offset for cell neighbour cells
  {
    auto offset = [this](IntIdx i, IntIdx j, IntIdx k) {
      return static_cast<IntIdx>(b_cells_.GetIdx(MIdx(i, j, k)).GetRaw() -
          b_cells_.GetIdx(MIdx(0, 0, 0)).GetRaw());
    };

    cell_neighbour_cell_offset_ = {{offset(-1, 0, 0), offset(1, 0, 0),
                                  offset(0, -1, 0), offset(0, 1, 0),
                                  offset(0, 0, -1), offset(0, 0, 1)}};
  }

  // Offset for cell neighbour nodes
  {
    auto offset = [this](IntIdx i, IntIdx j, IntIdx k) {
      return static_cast<IntIdx>(b_nodes_.GetIdx(MIdx(i, j, k)).GetRaw() -
          b_nodes_.GetIdx(MIdx(0, 0, 0)).GetRaw());
    };

    cell_neighbour_node_offset_ = {{offset(0, 0, 0), offset(1, 0, 0),
                                  offset(0, 1, 0), offset(1, 1, 0),
                                  offset(0, 0, 1), offset(1, 0, 1),
                                  offset(0, 1, 1), offset(1, 1, 1)}};
  }

  // Base for i-face and j-face neighbours
  for (Dir dir : {Dir::i, Dir::j, Dir::k}) {
    for (auto midx : BlockCells(mb, me - mb + dir)) {
      IdxFace idxface = b_faces_.GetIdx(midx, dir);
      ff_direction_[idxface] = dir;
      ff_neighbour_cell_[idxface] = {{
          b_cells_.GetIdx(midx - dir),
          b_cells_.GetIdx(midx)}};

      auto l = [this, &midx](IntIdx i, IntIdx j, IntIdx k) {
        return b_nodes_.GetIdx(midx + MIdx(i, j, k));
      };
      if (dir == Dir::i) {
        ff_neighbour_node_[idxface] = {{
            l(0,0,0), l(0,1,0), l(0,1,1), l(0,0,1)}};
      } else if (dir == Dir::j) {
        ff_neighbour_node_[idxface] = {{
            l(0,0,0), l(0,0,1), l(1,0,1), l(1,0,0)}};
      } else {
        // dir == Dir::k
        ff_neighbour_node_[idxface] = {{
            l(0,0,0), l(1,0,0), l(1,1,0), l(0,1,0)}};
      }

      if (midx[dir] == mb[dir]) {
        ff_neighbour_cell_[idxface][0] = IdxCell::None();
      }
      if (midx[dir] == me[dir]) {
        ff_neighbour_cell_[idxface][1] = IdxCell::None();
      }
      #ifdef PERX
      // adhoc for periodic in x
      if (midx[0] == mb[0] && dir == Dir::i) {
        ff_neighbour_cell_[idxface][0] = 
            b_cells_.GetIdx(MIdx(me[0] - 1, midx[1], midx[2]));
      }
      if (midx[0] == me[0] && dir == Dir::i) {
        ff_neighbour_cell_[idxface][1] = 
            b_cells_.GetIdx(MIdx(mb[0], midx[1], midx[2]));
      }
      #endif
      #ifdef PERZ
      // adhoc for periodic in x
      if (midx[2] == mb[2] && dir == Dir::k) {
        ff_neighbour_cell_[idxface][0] = 
            b_cells_.GetIdx(MIdx(midx[0], midx[1], me[2] - 1));
      }
      if (midx[2] == me[2] && dir == Dir::k) {
        ff_neighbour_cell_[idxface][1] = 
            b_cells_.GetIdx(MIdx(midx[0], midx[1], mb[2]));
      }
      #endif
    }
  }

  // Face centers, area
  for (auto idx : this->AllFaces()) {
    ff_center_[idx] = CalcCenter(idx);
    ff_area_[idx] = CalcArea(idx);
    ff_surface_[idx] = CalcSurface(idx);
  }

  // Cell centers, volume
  for (auto idx : this->AllCells()) {
    fc_center_[idx] = CalcCenter(idx);
    fc_volume_[idx] = CalcVolume(idx);
  }

  // Vect to cell
  for (Dir dir : {Dir::i, Dir::j, Dir::k}) {
    for (auto midx : BlockCells(mb, me - mb + dir)) {
      IdxFace idxface = b_faces_.GetIdx(midx, dir);
      IdxCell c;
      c = GetNeighbourCell(idxface, 0);
      if (!c.IsNone()) {
        ff_to_cell_[idxface][0] = GetCenter(c) - GetCenter(idxface);
      }
      c = GetNeighbourCell(idxface, 1);
      if (!c.IsNone()) {
        ff_to_cell_[idxface][1] = GetCenter(c) - GetCenter(idxface);
      }
      #ifdef PERX
      // adhoc for periodic in x
      if (midx[0] == mb[0] && dir == Dir::i) {
        auto pf = b_faces_.GetIdx(MIdx(me[0], midx[1], midx[2]), dir);
        auto pc = b_cells_.GetIdx(MIdx(me[0] - 1, midx[1], midx[2]));
        ff_to_cell_[idxface][0] = GetCenter(pc) - GetCenter(pf);
      }
      if (midx[0] == me[0] && dir == Dir::i) {
        auto pf = b_faces_.GetIdx(MIdx(mb[0], midx[1], midx[2]), dir);
        auto pc = b_cells_.GetIdx(MIdx(mb[0], midx[1], midx[2]));
        ff_to_cell_[idxface][1] = GetCenter(pc) - GetCenter(pf);
      }
      #endif
      
      #ifdef PERZ
      // adhoc for periodic in z
      if (midx[2] == mb[2] && dir == Dir::k) {
        auto pf = b_faces_.GetIdx(MIdx(midx[0], midx[1], me[2]), dir);
        auto pc = b_cells_.GetIdx(MIdx(midx[0], midx[1], me[2] - 1));
        ff_to_cell_[idxface][0] = GetCenter(pc) - GetCenter(pf);
      }
      if (midx[2] == me[2] && dir == Dir::k) {
        auto pf = b_faces_.GetIdx(MIdx(midx[0], midx[1], mb[2]), dir);
        auto pc = b_cells_.GetIdx(MIdx(midx[0], midx[1], mb[2]));
        ff_to_cell_[idxface][1] = GetCenter(pc) - GetCenter(pf);
      }
      #endif
    }
  }

  // Mark inner faces
  for (auto idxface : this->Faces()) {
    ff_is_inner_[idxface] = true; 
  }
}

// Create uniform mesh
// domain    - rectangle covering inner cells
// begin     - index of lower inner cells
// s        - number of inner cells in each direction
// hl        - number of halo layers
template <class Mesh>
Mesh InitUniformMesh(Rect<typename Mesh::Vect> domain,
                     typename Mesh::MIdx begin,
                     typename Mesh::MIdx s, int hl) {
  using Vect = typename Mesh::Vect;
  using MIdx = typename Mesh::MIdx;
  typename Mesh::BlockNodes b_nodes(begin - MIdx(hl), s + MIdx(1 + 2*hl));
  FieldNode<Vect> fn_node(b_nodes);
  for (auto midx : b_nodes) {
    IdxNode idxnode = b_nodes.GetIdx(midx);
    Vect unit = Vect(midx - b_nodes.GetBegin() - MIdx(hl)) / Vect(s);
    fn_node[idxnode] = domain.lb + unit * domain.GetDimensions();
  }
  return Mesh(b_nodes, fn_node, hl);
}

template <class Mesh>
Mesh InitUniformMesh(const Rect<typename Mesh::Vect>& domain,
                     typename Mesh::MIdx mesh_size) {
  using MIdx = typename Mesh::MIdx;
  return InitUniformMesh<Mesh>(domain, MIdx(0), mesh_size);
}


template <class Mesh>
class SearchMesh {
  const Mesh& mesh;
  using Scal = typename Mesh::Scal;
  using Vect = typename Mesh::Vect;
  static constexpr size_t dim = Mesh::dim;
  using IdxSeg = GIdx<20160212>; // Cartesian segment idx
  using BlockSeg = GBlock<IdxSeg, dim>;
  using Rect = geom::Rect<Vect>;
  using MIdx = GMIdx<dim>;
  template <class T>
  using FieldSeg = GField<T, IdxSeg>;
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
      for (auto midx : GBlock<size_t, dim>(
          mlb, mrt - mlb + MIdx(1))) {
        if (HaveCommonPoints(idxcell, midx)) {
          intersect_chains_[block_.GetIdx(midx)].push_back(idxcell);
          ++intersect_total;
        }
      }
    }

    num_intersect_seg_max = 0;
    fs_intersect_begin_.Reinit(GRange<IdxSeg>(0, block_.size() + 1));
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
    for (auto midx : GBlock<size_t, dim>(
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


/*
template <class _Scal, size_t _dim>
MeshStructured<_Scal, 3>::MeshStructured(
    const BlockNodes& b_nodes, 
    const FieldNode<Vect>& fn_node, int hl) {
  assert(false && "MeshStructured not implemented");
}
*/

} // namespace geom
