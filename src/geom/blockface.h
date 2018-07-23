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
    // Returns number of calls ++i_ for which
    // GetIdx(*++i_).GetRaw() == GetIdx(*i_).GetRaw() + 1
    // (i.e. increment of iterator equivalent to increment of index)
    size_t GetLite() const {
      MIdx xd(d_);
      assert(x_[0] + 1 <= o_->b_[0] + o_->cs_[0] + xd[0]);
      return o_->b_[0] + o_->cs_[0] + xd[0] - x_[0] - 1;
    }
    void IncLite(size_t a) {
      x_[0] += a;
    }
  };

  // b: begin, lower corner
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

template <size_t dim>
using GBlockFaces = GBlock<IdxFace, dim>;

