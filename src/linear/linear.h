// Created by Petr Karnakov on 25.04.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <map>
#include <memory>
#include <vector>

#include "debug/isnan.h"
#include "geom/mesh.h"

template <class Scal, class Idx>
struct GTerm {
  GTerm() = default;

  GTerm(Scal a, Idx idx) : a(a), idx(idx) {}

  bool operator<(const GTerm& o) const {
    return idx.GetRaw() < o.idx.GetRaw();
  }
  bool operator==(const GTerm& o) const {
    return idx.GetRaw() == o.idx.GetRaw();
  }

  Scal a;
  Idx idx;
};

// Linear expression.
// t[i].a * x[t[i].idx] + b
// Size_: maximum size of expression
template <class Scal_, class Idx_, size_t Size_>
class Expression {
 public:
  using Scal = Scal_;
  using Idx = Idx_;
  using Term = GTerm<Scal, Idx>;

  Expression() : b_(0), srt_(true), s_(0) {}

  Expression(const Expression& o) = default;
  Expression& operator=(const Expression& o) = default;

  explicit Expression(Scal b) : b_(b), srt_(true), s_(0) {}

  explicit Expression(const Term& t) : tt_{t}, b_(0.), srt_(true), s_(1) {}

  void Clear() {
    b_ = 0;
    srt_ = true;
    s_ = 0;
  }
  size_t size() const {
    return s_;
  }
  bool empty() const {
    return size() == 0;
  }
  Term& operator[](size_t k) {
    return tt_[k];
  }
  const Term& operator[](size_t k) const {
    return tt_[k];
  }
  void InsertTerm(const Term& term) {
    tt_[s_++] = term;
    srt_ = false;
  }
  void InsertTerm(Scal a, Idx idx) {
    InsertTerm(Term(a, idx));
  }
  Scal& Constant() {
    return b_;
  }
  Scal GetConstant() const {
    return b_;
  }
  void SetConstant(Scal b) {
    b_ = b;
  }
  template <class Fld>
  typename Fld::Value Evaluate(const Fld& u) const {
    typename Fld::Value r(b_); // result
    for (size_t i = 0; i < s_; ++i) {
      if (tt_[i].a != 0.) {
        r += u[tt_[i].idx] * tt_[i].a;
      }
    }
    return r;
  }
  Scal Coeff(Idx idx) const {
    for (size_t i = 0; i < s_; ++i) {
      if (tt_[i].idx == idx) {
        return tt_[i].a;
      }
    }
    return 0.;
  }
  Scal CoeffSum() const {
    Scal r = 0.;
    for (size_t i = 0; i < s_; ++i) {
      r += tt_[i].a;
    }
    return r;
  }
  void SortTerms(bool force = false) {
    if (!srt_ || force) {
      std::sort(tt_.begin(), tt_.begin() + s_);
      srt_ = true;
    }
  }
  Expression& operator*=(Scal k) {
    for (size_t i = 0; i < s_; ++i) {
      tt_[i].a *= k;
    }
    b_ *= k;
    return *this;
  }
  Expression operator*(Scal k) const {
    Expression tmp(*this);
    tmp *= k;
    return tmp;
  }
  Expression& operator/=(Scal k) {
    for (size_t i = 0; i < s_; ++i) {
      tt_[i].a /= k;
    }
    b_ /= k;
    return *this;
  }
  Expression operator/(Scal k) const {
    Expression tmp(*this);
    tmp /= k;
    return tmp;
  }
  Expression& operator+=(Expression o) {
    Expression& a = *this;
    Expression& b = o;

    if (a.empty()) {
      std::swap(a.tt_, b.tt_);
      std::swap(a.s_, b.s_);
    } else if (b.empty()) {
      // nop
    } else {
      a.SortTerms();
      b.SortTerms();
      Cont ab;
      size_t abs = 0;
      auto ai = a.tt_.cbegin(), bi = b.tt_.cbegin();
      const auto ae = a.tt_.cbegin() + a.s_;
      const auto be = b.tt_.cbegin() + b.s_;

      while (ai != ae || bi != be) {
        // current term
        auto t = (ai != ae && (bi == be || *ai < *bi)) ? *ai : *bi;
        Scal sum = 0;
        while (ai != ae && *ai == t) {
          sum += ai->a;
          ++ai;
        }
        while (bi != be && *bi == t) {
          sum += bi->a;
          ++bi;
        }
        ab[abs++] = Term(sum, t.idx);
      }
      std::swap(a.tt_, ab);
      std::swap(a.s_, abs);
      a.srt_ = true;
    }

    a.b_ += b.b_;

    return *this;
  }
  Expression operator+(Expression o) const {
    o += *this;
    return o;
  }
  Expression& operator-=(Expression o) {
    o *= -1.;
    *this += o;
    return *this;
  }
  Expression operator-(Expression o) const {
    Expression tmp(*this);
    tmp -= o;
    return tmp;
  }
  size_t Find(Idx idx) const {
    for (size_t i = 0; i < s_; ++i) {
      if (tt_[i].idx == idx) {
        return i;
      }
    }
    return -1;
  }
  // Replaces implicit term with explicit value.
  // (avoid changing the stencil, required by Hypre)
  void SetKnownValue(Idx idx, Scal v) {
    for (size_t i = 0; i < s_; ++i) {
      if (tt_[i].idx == idx) {
        b_ += v * tt_[i].a;
        tt_[i].a = 0.;
      }
    }
  }
  // Replaces expression with u[idx]=v
  void SetKnownValueDiag(Idx idx, Scal v) {
    for (size_t i = 0; i < s_; ++i) {
      tt_[i].a = (tt_[i].idx == idx ? 1. : 0.);
    }
    b_ = -v;
  }

 private:
  static constexpr size_t Size = Size_;
  using Cont = std::array<Term, Size>;

  template <class, class, size_t>
  friend class Expression;

  Cont tt_; // terms
  Scal b_; // free term
  bool srt_; // 1: terms sorted by idx
  size_t s_;
};

template <class Scal, class Idx, size_t Size>
std::ostream& operator<<(
    std::ostream& out, const Expression<Scal, Idx, Size>& e) {
  for (size_t i = 0; i < e.size(); ++i) {
    out << e[i].a << "*[" << e[i].idx.GetRaw() << "] + ";
  }
  out << e.GetConstant();
  return out;
}

template <class M, class S>
void PrintSystem(S& s, M& m, std::ostream& o) {
  auto bc = m.GetIndexCells();
  for (auto i : m.Cells()) {
    o << bc.GetMIdx(i) << " ";
    auto& e = s[i];
    for (size_t j = 0; j < e.size(); ++j) {
      o << e[j].a << "*[" << bc.GetMIdx(e[j].idx) << "] + ";
    }
    o << "\n";
  }
}

template <class M, class Expr, class Scal = typename M::Scal>
typename M::LS ConvertLs(
    const FieldCell<Expr>& fce, std::vector<Scal>& la, std::vector<Scal>& lb,
    std::vector<Scal>& lx, const M& m) {
  using LS = typename M::LS;
  using MIdx = typename M::MIdx;
  using IdxCell = IdxCell;
  auto& bc = m.GetIndexCells();
  LS l;
  // Get stencil from first inner cell
  {
    IdxCell c = *m.Cells().begin();
    auto& e = fce[c];
    for (size_t j = 0; j < e.size(); ++j) {
      MIdx dm = bc.GetMIdx(e[j].idx) - bc.GetMIdx(c);
      l.st.emplace_back(dm);
    }
  }

  size_t n = m.GetInBlockCells().size();
  la.resize(n * l.st.size());
  lb.resize(n);
  lx.resize(n);

  // fill matrix coeffs
  {
    size_t i = 0;
    for (auto c : m.Cells()) {
      auto& e = fce[c];
      for (size_t j = 0; j < e.size(); ++j) {
        // Check stencil
        if (e[j].idx != bc.GetIdx(bc.GetMIdx(c) + MIdx(l.st[j]))) {
          std::cerr << "***"
                    << " MIdx(c)=" << bc.GetMIdx(c)
                    << " MIdx(e[j].idx)=" << bc.GetMIdx(e[j].idx)
                    << " l.st[j]=" << MIdx(l.st[j]) << std::endl;

          throw std::runtime_error("ConvertLs: nonmatching stencils");
        }
        la[i] = e[j].a;
        ++i;
      }
    }
    assert(i == n * l.st.size());
  }

  // fill rhs and zero solution
  {
    size_t i = 0;
    for (auto c : m.Cells()) {
      auto& e = fce[c];
      lb[i] = -e.GetConstant();
      lx[i] = 0.;
      ++i;
    }
    assert(i == lb.size());
  }

  l.a = &la;
  l.b = &lb;
  l.x = &lx;
  return l;
}

// V is expression: v[0] * c + v[1] * cxm + ... v[6] * czp + v[7]
// lx: initial guess
template <
    class M, class Scal = typename M::Scal,
    class V = GVect<Scal, M::dim * 2 + 2>>
typename M::LS ConvertLsCompact(
    const FieldCell<V>& fce, std::vector<Scal>& la, std::vector<Scal>& lb,
    std::vector<Scal>& lx, const M& m) {
  using LS = typename M::LS;
  LS l;

  // stencil
  // XXX: adhoc, assumed order in m.GetCell(c)
  /*
  l.st.emplace_back(0, 0, 0);   // 0
  l.st.emplace_back(-1, 0, 0);  // 1
  l.st.emplace_back(1, 0, 0);   // 2
  l.st.emplace_back(0, -1, 0);  // 3
  l.st.emplace_back(0, 1, 0);   // 4
  l.st.emplace_back(0, 0, -1);  // 5
  l.st.emplace_back(0, 0, 1);   // 6
  */

  l.st.emplace_back(0, 0, -1); // 5
  l.st.emplace_back(0, -1, 0); // 3
  l.st.emplace_back(-1, 0, 0); // 1
  l.st.emplace_back(0, 0, 0); // 0
  l.st.emplace_back(1, 0, 0); // 2
  l.st.emplace_back(0, 1, 0); // 4
  l.st.emplace_back(0, 0, 1); // 6
  assert(l.st.size() == V::dim - 1);

  size_t n = m.GetInBlockCells().size(); // number of equations
  la.resize(n * l.st.size());
  lb.resize(n);
  lx.resize(n);

  // fill matrix coeffs
  {
    size_t i = 0;
    for (auto c : m.Cells()) {
      auto& e = fce[c];
      la[i++] = e[5];
      la[i++] = e[3];
      la[i++] = e[1];
      la[i++] = e[0];
      la[i++] = e[2];
      la[i++] = e[4];
      la[i++] = e[6];
    }
    assert(i == n * l.st.size());
  }

  // fill rhs and zero solution
  {
    size_t i = 0;
    for (auto c : m.Cells()) {
      auto& e = fce[c];
      lb[i] = -e[V::dim - 1];
      ++i;
    }
    assert(i == lb.size());
  }

  l.a = &la;
  l.b = &lb;
  l.x = &lx;
  return l;
}

template <class Scal, class Idx, size_t Size>
bool IsNan(const Expression<Scal, Idx, Size>& e) {
  if (IsNan(e.GetConstant())) {
    return true;
  }
  for (size_t i = 0; i < e.size(); ++i) {
    if (IsNan(e[i].a)) {
      return true;
    }
  }
  return false;
}
