#pragma once

#include <utility>

#include "block.h"
#include "dir.h"

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
    MIdx w_;
    Dir d_;

   public:
    explicit iterator(const GBlock* o, MIdx w, Dir d)
        : o_(o), w_(w), d_(d)
    {}
    iterator& operator++() {
      MIdx wd(d_);
      for (size_t n = 0; n < dim; ++n) {
        ++w_[n]; // increment current
        if (w_[n] == o_->b_[n] + o_->cs_[n] + wd[n]) { // if end reached
          if (n < dim - 1) {  
            w_[n] = o_->b_[n];  // reset to begin
          } else {  // if end reached for last dim
            if (size_t(d_) < dim - 1) {
              d_ = Dir(size_t(d_) + 1);
              w_ = o_->b_;
            }
          }
        } else {
          break;
        }
      }  
      return *this;
    }
    bool operator==(const iterator& o) const {
      return w_ == o.w_ && d_ == o.d_;
    }
    bool operator!=(const iterator& o) const {
      return !(*this == o);
    }
    std::pair<MIdx, Dir> operator*() const {
      return std::make_pair(w_, d_);
    }
    // Returns number of calls ++i_ for which
    // GetIdx(*++i_).GetRaw() == GetIdx(*i_).GetRaw() + 1
    // (i.e. increment of iterator equivalent to increment of index)
    size_t GetLite() const {
      MIdx wd(d_);
      assert(w_[0] + 1 <= o_->b_[0] + o_->cs_[0] + wd[0]);
      return o_->b_[0] + o_->cs_[0] + wd[0] - w_[0] - 1;
    }
    void IncLite(size_t a) {
      w_[0] += a;
    }
  };

  // b: begin, lower corner cell index
  // cs: cells size
  GBlock(MIdx b, MIdx cs) 
      : b_(b), cs_(cs)
  {
    // number of faces
    nfa_ = 0;
    for (size_t d = 0; d < dim; ++d) {
      nf_[d] = CalcNumFaces(d);
      nfa_ += nf_[d];
    }
    // cumulative number of faces
    nfc_[0] = 0;
    for (size_t d = 1; d < dim; ++d) {
      nfc_[d] = nfc_[d - 1] + nf_[d - 1];
    }
  }
  GBlock(MIdx cs) : GBlock(MIdx(0), cs) {}
  GBlock() : GBlock(MIdx(0), MIdx(0)) {}
  size_t size() const {
    return nfa_;
  }
  operator GRange<Idx>() const {
    return GRange<Idx>(0, size());
  }
  Idx GetIdx(MIdx w, Dir d) const {
    size_t r = nfc_[size_t(d)] + GetFlat(w, size_t(d));
    return Idx(r);
  }
  Idx GetIdx(std::pair<MIdx, Dir> p) const {
    return GetIdx(p.first, p.second);
  }
  Dir GetDir(Idx f) const {
    size_t r = f.GetRaw();
    size_t d = 0;
    while (r >= GetNumFaces(d) && d < dim) {
      r -= GetNumFaces(d);
      ++d;
    }
    return Dir(d);
  }
  std::pair<MIdx, Dir> GetMIdxDir(Idx f) const {
    size_t r = f.GetRaw();
    size_t d = 0;
    while (r >= GetNumFaces(d) && d < dim) {
      r -= GetNumFaces(d);
      ++d;
    }
    return {GetMIdxFromOffset(r, d), Dir(d)};
  }
  MIdx GetMIdx(Idx f) const {
    return GetMIdxDir(f).first;
  }
  MIdx GetSize() const {
    return cs_;
  }
  iterator begin() const {
    return iterator(this, b_, Dir(0));
  }
  iterator end() const {
    MIdx w = b_;
    auto ld = dim - 1; // last direction
    w[ld] = (b_ + cs_)[ld] + 1;
    return iterator(this, w, Dir(ld));
  }

 private:
  const MIdx b_;
  const MIdx cs_; // cells size
  MIdx nf_; // number of faces in each direction
  size_t nfa_; // total number of faces
  MIdx nfc_; // cumulative number of faces up to direction: sum(nf_[0:d])
  size_t GetNumFaces(size_t d) const {
    return nf_[d];
  }
  size_t CalcNumFaces(size_t d) const {
    MIdx bcs = cs_;
    ++bcs[d];
    return bcs.prod();
  }
  size_t GetFlat(MIdx w, size_t d) const {
    w -= b_;
    MIdx bcs = cs_;
    ++bcs[d];
    size_t r = 0;
    for (size_t i = dim; i != 0; ) {
      --i;
      r *= bcs[i];
      r += w[i];
    }
    return r;
  }
  // r: offset from first face with direction d
  // d: direction
  MIdx GetMIdxFromOffset(size_t r, size_t d) const {
    MIdx bcs = cs_;
    ++bcs[d];
    MIdx w;
    for (size_t i = 0; i < dim; ++i) {
      w[i] = r % bcs[i];
      r /= bcs[i];
    }
    return b_ + w;
  }
};

// Specialization for IdxFace
// Enumeration order (to more frequent): dir, z, y, x
template <size_t dim_>
class GBlockPad<IdxFace, dim_> {
 public:
  static constexpr size_t dim = dim_;
  using Idx = IdxFace;
  using MIdx = GMIdx<dim>;
  using Dir = GDir<dim>;

  // b: begin, lower corner
  // rs: raw size (underlying block that fits all faces_)
  GBlockPad(MIdx b, MIdx rs) 
      : b_(b), rs_(rs), br_(b_, rs_)
  {
    // number of faces in each direction
    nfa_ = 0;
    for (size_t d = 0; d < dim; ++d) {
      nf_[d] = CalcNumFaces(d);
      nfa_ += nf_[d];
    }
    // cumulative number of faces
    nfc_[0] = 0;
    for (size_t d = 1; d < dim; ++d) {
      nfc_[d] = nfc_[d - 1] + nf_[d - 1];
    }
  }
  GBlockPad(MIdx rs) : GBlockPad(MIdx(0), rs) {}
  GBlockPad() : GBlockPad(MIdx(0), MIdx(0)) {}
  size_t size() const {
    return nfa_;
  }
  operator GRange<Idx>() const {
    return GRange<Idx>(0, size());
  }
  Idx GetIdx(MIdx w, Dir d) const {
    size_t r = GetBase(size_t(d)) + br_.GetIdx(w);
    return Idx(r);
  }
  Idx GetIdx(std::pair<MIdx, Dir> p) const {
    return GetIdx(p.first, p.second);
  }
  Dir GetDir(Idx f) const {
    size_t r = f.GetRaw();
    size_t d = 0;
    while (r >= GetNumFaces(d) && d < dim) {
      r -= GetNumFaces(d);
      ++d;
    }
    return Dir(d);
  }
  std::pair<MIdx, Dir> GetMIdxDir(Idx f) const {
    size_t r = f.GetRaw();
    size_t d = 0;
    while (r >= GetNumFaces(d) && d < dim) {
      r -= GetNumFaces(d);
      ++d;
    }
    return {br_.GetMIdx(r), Dir(d)};
  }
  MIdx GetMIdx(Idx f) const {
    size_t r = f.GetRaw();
    size_t d = 0;
    while (r >= GetNumFaces(d) && d < dim) {
      r -= GetNumFaces(d);
      ++d;
    }
    return br_.GetMIdx(r);
  }

 private:
  const MIdx b_;
  const MIdx rs_; // size of cells with padding
  const GBlock<size_t, dim> br_; // block of cells with padding
  MIdx nf_; // number of faces in each direction with padding
  size_t nfa_; // total number of indices with padding
  MIdx nfc_; // cumulative number of faces up to direction: sum(nf_[0:d])
             // (base index for faces of direction d)
  size_t GetNumFaces(size_t d) const {
    return nf_[d];
  }
  size_t GetBase(size_t d) const {
    return nfc_[d];
  }
  size_t CalcNumFaces(size_t /*d*/) const {
    return br_.size();
  }
};

template <size_t dim>
using GBlockFaces = GBlock<IdxFace, dim>;
