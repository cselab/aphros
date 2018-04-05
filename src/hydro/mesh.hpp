#pragma once

#include <vector>
#include <array>
#include <cstddef>
#include <map>
#include <cassert>
#include <iostream>
#include <utility>

#include "suspender.h"
#include "vect.hpp"
#include "idx.h"
#include "range.h"
#include "field.h"
#include "map.h"
#include "block.h"
#include "dir.h"

namespace geom {



#define BLOCKFACE_DZYX
//#define BLOCKFACE_ZYXD

#ifdef BLOCKFACE_DZYX
// Specialization for IdxFace
// Enumeration order (to more frequent): dir, z, y, x
template <size_t dim_>
class GBlock<IdxFace, dim_> {
 public:
  static constexpr size_t dim = dim_;
  using Idx = IdxFace;
  using MIdx = GMIdx<dim>;
  using Dir = GDir<dim>;

  class iterator {
    const GBlock* o_; // owner
    MIdx x_;
    Dir d_;

   public:
    explicit iterator(const GBlock* o, MIdx x, Dir d)
        : o_(o), x_(x), d_(d)
    {}
    iterator& operator++() {
      MIdx xd(d_);
      for (size_t n = 0; n < dim; ++n) {
        ++x_[n]; // increment current
        if (x_[n] == o_->b_[n] + o_->cs_[n] + xd[n]) { // if end reached
          if (n < dim - 1) {  
            x_[n] = o_->b_[n];  // reset to begin
          } else {  // if end reached for last dim
            if (size_t(d_) < dim - 1) {
              d_ = Dir(size_t(d_) + 1);
              x_ = o_->b_;
            }
          }
        } else {
          break;
        }
      }  
      return *this;
    }
    bool operator==(const iterator& o) const {
      return x_ == o.x_ && d_ == o.d_;
    }
    bool operator!=(const iterator& o) const {
      return !(*this == o);
    }
    std::pair<MIdx, Dir> operator*() const {
      return std::make_pair(x_, d_);
    }
  };

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
    for (size_t i = 0; i < size_t(dir); ++i) {
      raw += GetNumFaces(i);
    }
    raw += GetFlat(midx, dir);
    return Idx(raw);
  }
  Idx GetIdx(std::pair<MIdx, Dir> p) const {
    return GetIdx(p.first, p.second);
  }
  MIdx GetMIdx(Idx idx) const {
    return GetMIdxDir(idx).first;
  }
  Dir GetDir(Idx idx) const {
    return GetMIdxDir(idx).second;
  }
  MIdx GetSize() const {
    return cs_;
  }
  iterator begin() const {
    return iterator(this, b_, Dir(0));
  }
  iterator end() const {
    MIdx x = b_;
    auto ld = dim - 1; // last direction
    x[ld] = (b_ + cs_)[ld] + 1;
    return iterator(this, x, Dir(ld));
  }

 private:
  MIdx b_;
  MIdx cs_; // cells size
  size_t GetNumFaces(Dir dir) const {
    MIdx bcs = cs_;
    ++bcs[size_t(dir)];
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
    ++bcs[size_t(dir)];
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
    ++bcs[size_t(dir)];
    MIdx midx;
    for (size_t i = 0; i < dim; ++i) {
      midx[i] = raw % bcs[i];
      raw /= bcs[i];
    }
    return b_ + midx;
  }

 public:
  std::pair<MIdx, Dir> GetMIdxDir(Idx idxface) const {
    size_t raw = idxface.GetRaw();
    size_t diridx = 0;
    while (raw >= GetNumFaces(diridx) && diridx < dim) {
      raw -= GetNumFaces(diridx);
      ++diridx;
    }
    Dir dir(diridx); return {GetMIdxFromOffset(raw, dir), dir};
  }
};
#endif

#ifdef BLOCKFACE_ZYXD
// Specialization for IdxFace
// Enumeration order (to more frequent): z, y, x, dir
template <size_t dim_>
class GBlock<IdxFace, dim_> { // [s]andbox
 public:
  static constexpr size_t dim = dim_;
  using Idx = IdxFace;
  using MIdx = GMIdx<dim>;
  using Dir = GDir<dim>;
  
  static_assert(dim == 3, "GBlock<IdxFace,...> implemented only for dim=3");

  class iterator {
    const GBlock* o_; // owner
    MIdx x_;
    Dir d_;

   public:
    explicit iterator(const GBlock* o, MIdx x, Dir d)
        : o_(o), x_(x), d_(d)
    {}
    iterator& operator++() {
      auto& s = o_->cs_;
      auto& x = x_;
      auto& b = o_->b_;
      auto& e = o_->e_;
      auto& d = d_;

      if (x[2] == e[2]) { // z-boundary
        ++x[0];
        if (x[0] == e[0]) {
          x[0] = b[0];
          ++x[1];  // end would be: (b[0], e[1], e[2])
        }
      } else if (x[1] == e[1]) { // y-edge
        ++x[0];
        if (x[0] == e[0]) { // next x-end
          ++x[2];
          x[1] = b[1];
          x[0] = b[0];
          d = Dir(0);
          if (x[2] == e[2]) { // next z-edge
            d = Dir(2);
          } 
        } 
      } else if (x[0] == e[0]) { // x-edge
          x[0] = b[0];
          ++x[1];
          if (x[1] == e[1]) { // next y-edge
            d = Dir(1);
          }
      } else if (d == Dir(2)) { // dirz
        ++x[0];
        d = Dir(0);
      } else {
        d = Dir(size_t(d) + 1);
      }

      return *this;
    }
    bool operator==(const iterator& o) const {
      return x_ == o.x_ && d_ == o.d_;
    }
    bool operator!=(const iterator& o) const {
      return !(*this == o);
    }
    std::pair<MIdx, Dir> operator*() const {
      return std::make_pair(x_, d_);
    }
  };

  GBlock(MIdx begin, MIdx cells)
      : b_(begin), cs_(cells), e_(b_ + cs_)
  {}
  GBlock(MIdx cells)
      : GBlock(MIdx(0), cells)
  {}
  GBlock()
      : GBlock(MIdx(0), MIdx(0))
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
  Idx GetIdx(MIdx x, Dir d) const {
    x -= b_;
    size_t r = 0;
    auto& s = cs_;
    r += (s[1] * s[0] + (s[1] + 1) * s[0] + (s[0] + 1) * s[1]) * x[2]; 
    if (x[2] == s[2]) {
      r += s[0] * x[1] + x[0];   // z-boundary
    } else {
      r += (3 * s[0] + 1) * x[1];
      if (x[1] == s[1]) {
        r += x[0]; // y-boundary
      } else {
        r += 3 * x[0] + size_t(d);
      }
    }
    return Idx(r);
  }
  Idx GetIdx(std::pair<MIdx, Dir> p) const {
    return GetIdx(p.first, p.second);
  }
  MIdx GetMIdx(Idx idx) const {
    return GetMIdxDir(idx).first;
  }
  Dir GetDir(Idx idx) const {
    return GetMIdxDir(idx).second;
  }
  MIdx GetSize() const {
    return cs_;
  }
  iterator begin() const {
    return iterator(this, b_, Dir(0));
  }
  iterator end() const {
    return iterator(this, MIdx(b_[0], e_[1], e_[2]), Dir(2));
  }
  std::pair<MIdx, Dir> GetMIdxDir(Idx i) const {
    size_t r = i.GetRaw();
    MIdx x;
    Dir d;
    auto& s = cs_;
    const size_t n = size();
    if (r >= n - s[0] * s[1]) {  // z-boundary
      r -= n - s[0] * s[1];
      x[2] = s[2];
      x[1] = r / s[0];
      x[0] = r % s[0];
      d = Dir(2);
    } else {
      auto w = (s[1] * s[0] + (s[1] + 1) * s[0] + (s[0] + 1) * s[1]);
      x[2] = r / w;
      r = r % w;
      if (r >= w - s[0]) { // y-boundary
        r -= w - s[0];
        x[1] = s[1];
        x[0] = r;
        d = Dir(1);
      } else {
        auto q = 3 * s[0] + 1;
        x[1] = r / q;
        r = r % q;
        x[0] = r / 3;
        d = Dir(r % 3);
      }
    }
    return std::make_pair(b_ + x, d);
  }

 private:
  MIdx b_;  // begin
  MIdx cs_; // cells size
  MIdx e_;  // end
  size_t GetNumFaces(Dir dir) const {
    MIdx bcs = cs_;
    ++bcs[size_t(dir)];
    size_t res = 1;
    for (size_t i = 0; i < dim; ++i) {
      res *= bcs[i];
    }
    return res;
  }
  size_t GetNumFaces(size_t i) const {
    return GetNumFaces(Dir(i));
  }
};
#endif

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


#if 0
// Specialization for IdxFace
// Implementation via convertion Idx->MIdx->Idx
// (slow)
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
#endif


#if 0
// Specialization for IdxFace
// Implementation with optimization: 
// don't recompute until reaching the end in x direction
// (causes segfault for Hydro)
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
    IdxFace r_;   // current result
    I in_;        // next iterator to trigger recompute of r_ via MIdx
   public:
    IdxFace IdxFromInner(I i) {
      auto x = bi_.GetMIdx(*i); 
      auto d = bi_.GetDir(*i); 
      return ba_.GetIdx(x, d); 
    }
    I GetLast(I i) {
      using Dir = typename B::Dir;
      auto x = bi_.GetMIdx(*i); 
      auto d = bi_.GetDir(*i); 
      x[0] += bi_.GetSize()[0];
      if (d != Dir::i) {
        x[0] -= 1;
      }
      return I(bi_.GetIdx(x, d).GetRaw());
    }
    explicit iterator(const B& ba, const B& bi, const I& i)
        : ba_(ba), bi_(bi), i_(i), r_(IdxFromInner(i_)), in_(i_)
    {
        in_ = GetLast(in_);
    }
    iterator& operator++() {
      if (i_ == in_) {
        ++i_;
        r_ = IdxFromInner(i_);
        ++in_;
        in_ = GetLast(in_);
      } else {
        ++i_;
        r_.AddRaw(1);
      }
      return *this;
    }
    /*
    iterator& operator--() {
      --i_;
      return *this;
    }
    */
    bool operator==(const iterator& o /*other*/) const {
      return i_ == o.i_;
    }
    bool operator!=(const iterator& o /*other*/) const {
      return !(*this == o);
    }
    IdxFace operator*() const {
      return r_;
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
#endif

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
  FieldNode<Vect> fn_node_;
  BlockCells b_cells_;
  BlockFaces b_faces_;
  // inner 
  BlockCells b_incells_;
  BlockFaces b_infaces_;
  BlockNodes b_innodes_;
  // support
  BlockCells b_sucells_;
  BlockFaces b_sufaces_;
  BlockNodes b_sunodes_;


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
  // TODO: rename to GetBlockAllCells()
  const BlockCells& GetBlockCells() const {
    return b_cells_;
  }
  const BlockFaces& GetBlockFaces() const {
    return b_faces_;
  }
  const BlockNodes& GetBlockNodes() const {
    return b_nodes_;
  }
  // TODO: rename to GetBlockCells()
  const BlockCells& GetInBlockCells() const {
    return b_incells_;
  }
  const BlockFaces& GetInBlockFaces() const {
    return b_infaces_;
  }
  const BlockNodes& GetInBlockNodes() const {
    return b_innodes_;
  }
  // TODO: rename to GetBlockCells()
  const BlockCells& GetSuBlockCells() const {
    return b_sucells_;
  }
  const BlockFaces& GetSuBlockFaces() const {
    return b_sufaces_;
  }
  const BlockNodes& GetSuBlockNodes() const {
    return b_sunodes_;
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
  // Type-cast to GRange required for GField initialization
  template <class Idx>
  operator GRange<Idx>() const {
    return GetAll<Idx>(); 
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

  // Returns range of support indices
  template <class Idx>
  GRangeIn<Idx, dim> GetSu() const { 
    return GetSu((Idx*)0);
  }
  GRangeIn<IdxCell, dim> GetSu(IdxCell*) const {
    return GRangeIn<IdxCell, dim>(GetBlockCells(), GetSuBlockCells());
  }
  GRangeIn<IdxFace, dim> GetSu(IdxFace*) const {
    return GRangeIn<IdxFace, dim>(GetBlockFaces(), GetSuBlockFaces());
  }
  GRangeIn<IdxNode, dim> GetSu(IdxNode*) const {
    return GRangeIn<IdxNode, dim>(GetBlockNodes(), GetSuBlockNodes());
  }
  GRangeIn<IdxCell, dim> SuCells() const {
    return GetSu<IdxCell>();
  }
  GRangeIn<IdxFace, dim> SuFaces() const {
    return GetSu<IdxFace>();
  }
  GRangeIn<IdxNode, dim> SuNodes() const {
    return GetSu<IdxNode>();
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
    for (auto idxcell : Cells()) {
      if (point.dist(GetCenter(idxcell)) <
          point.dist(GetCenter(idxcell_nearest))) {
        idxcell_nearest = idxcell;
      }
    }
    return idxcell_nearest;
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
  struct LS { // linear system ax=r
    std::vector<MIdx> st; // stencil
    std::vector<Scal>* a; // matrix coeffs of size n * st.size()
    std::vector<Scal>* b; // rhs of size n
    std::vector<Scal>* x; // solution and initial guess of size n
    enum class T { gen, symm };
    T t = T::gen;
  };
  Sem GetSem(std::string name="") {
    return susp_.GetSem(name);
  }
  bool Pending() const {
    return susp_.Pending();
  }
  const Suspender& GetSusp() const {
    return susp_;
  } 
  std::string GetCurName() const {
    return susp_.GetCurName();
  }
  size_t GetDepth() const {
    return susp_.GetDepth();
  }
  void Comm(FieldCell<Scal>* u) {
    vcm_.push_back(u);
  }
  void Dump(FieldCell<Scal>* u, std::string name) {
    vd_.push_back(std::make_pair(u, name));
  }
  void Reduce(Scal* u, std::string op) {
    vrd_.push_back(std::make_pair(u, op));
  }
  void Solve(const LS& ls) {
    vls_.push_back(ls);
  }
  const std::vector<FieldCell<Scal>*>& GetComm() const {
    return vcm_;
  }
  const std::vector<std::pair<FieldCell<Scal>*, std::string>>& GetDump() const {
    return vd_;
  }
  void ClearComm() {
    vcm_.clear();
  }
  void ClearDump() {
    vd_.clear();
  }
  const std::vector<std::pair<Scal*, std::string>>& GetReduce() const {
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
  std::vector<std::pair<FieldCell<Scal>*, std::string>> vd_; // fields for dump
  std::vector<std::pair<Scal*, std::string>> vrd_; // scalars for reduce
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
    , b_sucells_(b_cells_.GetBegin() + MIdx(1), 
                 b_cells_.GetDimensions() - MIdx(2))
    , b_sufaces_(b_cells_.GetBegin() + MIdx(1), 
                 b_cells_.GetDimensions() - MIdx(2))
    , b_sunodes_(b_nodes_.GetBegin() + MIdx(1), 
                 b_nodes_.GetDimensions() - MIdx(2))
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
    for (auto midx : BlockCells(mb, me - mb + MIdx(dir))) {
      IdxFace idxface = b_faces_.GetIdx(midx, dir);
      ff_direction_[idxface] = dir;
      ff_neighbour_cell_[idxface] = {{
          b_cells_.GetIdx(midx - MIdx(dir)),
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

      const size_t sd(dir);
      if (midx[sd] == mb[sd]) {
        ff_neighbour_cell_[idxface][0] = IdxCell::None();
      }
      if (midx[sd] == me[sd]) {
        ff_neighbour_cell_[idxface][1] = IdxCell::None();
      }
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
    for (auto midx : BlockCells(mb, me - mb + MIdx(dir))) {
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


} // namespace geom
