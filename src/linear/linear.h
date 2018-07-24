#pragma once

#include <cassert>
#include <algorithm>
#include <vector>
#include <array>
#include <iostream>
#include <cstdint>
#include <memory>
#include <map>

#include "geom/mesh.h"

namespace solver {

template <class Scal, class Idx>
struct Term {
  Scal coeff;
  Idx idx;
  Term() = default;
  Term(Scal coeff, Idx idx)
      : coeff(coeff)
      , idx(idx)
  {}
  bool operator<(const Term& other) const {
    return idx.GetRaw() < other.idx.GetRaw();
  }
  bool operator==(const Term& other) const {
    return idx.GetRaw() == other.idx.GetRaw();
  }
};

template <class _Scal, class _Idx, size_t _ExprSize>
class Expression {
 public:
  using Scal = _Scal;
  using Idx = _Idx;
  using TermType = Term<Scal, Idx>;

 private:
  static constexpr size_t ExprSize = _ExprSize;
  using TermContainer = std::array<TermType, ExprSize>;
      //boost::container::static_vector<TermType, ExprSize>;
  TermContainer terms_;
  Scal constant_;
  bool sorted_;
  size_t size_;
  template <class OtherScal, class OtherIdx, size_t OtherExprSize>
  friend class Expression;

 public:
  Expression()
      : constant_(0)
      , sorted_(true)
      , size_(0)
  {}
  template <class OtherScal, class OtherIdx, size_t OtherExprSize>
  Expression(const Expression<OtherScal, OtherIdx, OtherExprSize>& other)
      : terms_(other.terms_),
        constant_(other.constant_),
        sorted_(other.sorted_),
        size_(other.size_) { }
  explicit Expression(Scal constant)
      : constant_(constant)
      , sorted_(true)
      , size_(0)
  {}
  explicit Expression(const TermType& term)
      : terms_{term},
        constant_(0.),
        sorted_(true),
        size_(1)
  {}
  void Clear() {
    constant_ = 0;
    sorted_ = true;
    size_ = 0;
  }
  size_t size() const {
    return size_;
  }
  bool empty() const {
    return size_ == 0;
  }
  TermType& operator[](size_t k) {
#ifdef __RANGE_CHECK
    assert(k >=0 && k < size_);
#endif
    return terms_[k];
  }
  const TermType& operator[](size_t k) const {
#ifdef __RANGE_CHECK
    assert(k >=0 && k < size_);
#endif
    return terms_[k];
  }
  void InsertTerm(const TermType& term) {
    terms_[size_++] = term;
    sorted_ = false;
  }
  void InsertTerm(Scal coeff, Idx idx) {
    InsertTerm(TermType(coeff, idx));
  }
  Scal& Constant() {
    return constant_;
  }
  Scal GetConstant() const {
    return constant_;
  }
  void SetConstant(Scal constant) {
    constant_ = constant;
  }
  template <class Field, class Value = typename Field::Value>
  Value Evaluate(const Field& field) const {
    Value res(constant_);
    for (size_t i = 0; i < size_; ++i) {
      res += field[terms_[i].idx] * terms_[i].coeff;
    }
    return res;
  }
  Scal Coeff(Idx idx) const {
    for (size_t i = 0; i < size_; ++i) {
      if (terms_[i].idx == idx) {
        return terms_[i].coeff;
      }
    }
    return 0.;
  }
  Scal CoeffSum() const {
    Scal res = 0.;
    for (size_t i = 0; i < size_; ++i) {
      res += terms_[i].coeff;
    }
    return res;
  }
  void SortTerms(bool force = false) {
    if (!sorted_ || force) {
      std::sort(terms_.begin(), terms_.begin() + size_);
      sorted_ = true;
    }
  }
  Expression& operator*=(Scal k) {
    for (size_t i = 0; i < size_; ++i) {
      terms_[i].coeff *= k;
    }
    constant_ *= k;
    return *this;
  }
  Expression operator*(Scal k) const {
    Expression tmp(*this);
    tmp *= k;
    return tmp;
  }
  Expression& operator/=(Scal k) {
    for (size_t i = 0; i < size_; ++i) {
      terms_[i].coeff /= k;
    }
    constant_ /= k;
    return *this;
  }
  Expression operator/(Scal k) const {
    Expression tmp(*this);
    tmp /= k;
    return tmp;
  }
  Expression& operator+=(Expression other) {
    Expression& a = *this;
    Expression& b = other;

    if (a.empty()) {
      std::swap(a.terms_, b.terms_);
      std::swap(a.size_, b.size_);
    } else if (b.empty()) {
      // nop
    } else {
      a.SortTerms();
      b.SortTerms();
      TermContainer ab;
      size_t absize = 0;
      auto a_it = a.terms_.cbegin(), b_it = b.terms_.cbegin();
      const auto a_end = a.terms_.cbegin() + a.size_;
      const auto b_end = b.terms_.cbegin() + b.size_;

      while (a_it != a_end || b_it != b_end) {
        auto sample = (a_it != a_end && (b_it == b_end || *a_it < *b_it))
            ? *a_it
            : *b_it;
        Scal sum = 0;
        while (a_it != a_end && *a_it == sample) {
          sum += a_it->coeff;
          ++a_it;
        }
        while (b_it != b_end && *b_it == sample) {
          sum += b_it->coeff;
          ++b_it;
        }
        ab[absize++] = TermType(sum, sample.idx);
      }
      std::swap(a.terms_, ab);
      std::swap(a.size_, absize);
      a.sorted_ = true;
    }

    a.constant_ += b.constant_;

    return *this;
  }
  Expression operator+(Expression other) const {
    other += *this;
    return other;
  }
  Expression& operator-=(Expression other) {
    other *= -1.;
    *this += other;
    return *this;
  }
  Expression operator-(Expression other) const {
    Expression tmp(*this);
    tmp -= other;
    return tmp;
  }
  size_t Find(Idx idx) const {
    for (size_t i = 0; i < size_; ++i) {
      if (terms_[i].idx == idx) {
        return i;
      }
    }
    return -1;
  }
  /* // Replaces terms list removing target
  void SetKnownValue(Idx idx, Scal value) {
    TermContainer new_terms;
    size_t newsize = 0;
    for (size_t i = 0; i < size_; ++i) {
      if (terms_[i].idx == idx) {
        constant_ += value * terms_[i].coeff;
      } else {
        new_terms[newsize++] = terms_[i];
      }
    }
    std::swap(terms_, new_terms);
    std::swap(size_, newsize);
  }
  */
  // Keeps terms list setting coeff to zero
  // (avoid changing the stencil, required by Hypre)
  void SetKnownValue(Idx idx, Scal value) {
    for (size_t i = 0; i < size_; ++i) {
      if (terms_[i].idx == idx) {
        constant_ += value * terms_[i].coeff;
        terms_[i].coeff = 0.;
      }
    }
  }
  // Replaces expression with u[idx]=value
  void SetKnownValueDiag(Idx idx, Scal value) {
    for (size_t i = 0; i < size_; ++i) {
      terms_[i].coeff = (terms_[i].idx == idx ? 1. : 0.);
    }
    constant_ = -value;
  }
  // Remove terms to get triangular matrix
  void RestrictTerms(Idx lower, Idx upper) {
    TermContainer new_terms;
    size_t newsize = 0;
    for (size_t i = 0; i < size_; ++i) {
      if (lower.GetRaw() <= terms_[i].idx.GetRaw() &&
          terms_[i].idx.GetRaw() <= upper.GetRaw()) {
        new_terms[newsize++] = terms_[i];
      }
    }
    std::swap(terms_, new_terms);
    std::swap(size_, newsize);
  }
};

template <class Scal, class Idx, size_t ExprSize>
std::ostream& operator<<(std::ostream& out,
                         const Expression<Scal, Idx, ExprSize>& expr) {
  for (size_t i = 0; i < expr.size(); ++i) {
    out << expr[i].coeff << "*[" << expr[i].idx.GetRaw() << "] + ";
  }
  out << expr.GetConstant();
  return out;
}

template <class M, class S>
void PrintSystem(S& s, M& m, std::ostream& o) {
  auto bc = m.GetIndexCells();
  for (auto i : m.Cells()) {
    o << bc.GetMIdx(i) << " "; 
    auto& e = s[i];
    for (size_t j = 0; j < e.size(); ++j) { 
      o << e[j].coeff << "*[" << bc.GetMIdx(e[j].idx) << "] + ";
    }
    o << "\n";
  }
}


template <class System, class Result>
void Transpose(const System& system, Result& result) {
  using ResExpr = typename Result::Value;
  using Term = typename ResExpr::TermType;
  result.Reinit(system.GetRange(), ResExpr());
  for (auto idx : system.GetRange()) {
    auto& eqn = system[idx];
    for (size_t i = 0; i < eqn.size(); ++i) {
      result[eqn[i].idx] += ResExpr(Term(eqn[i].coeff, idx));
    }
  }
}

template <class System>
void Normalize(System& system) {
  for (auto idx : system.GetRange()) {
    auto& eqn = system[idx];
    decltype(eqn[0].coeff) sum = 0.;
    for (size_t i = 0; i < eqn.size(); ++i) {
      sum += eqn[i].coeff;
    }
    for (size_t i = 0; i < eqn.size(); ++i) {
      eqn[i].coeff /= sum;
    }
    eqn.Constant() /= sum;
  }
}

template <class Scal, class Idx, class Expr>
class LinearSolver {
  template <class T>
  using Field = GField<T, Idx>;

 public:
  virtual Field<Scal> Solve(const Field<Expr>&) = 0;
  virtual ~LinearSolver() {}
};

class LinearSolverFactoryGeneric {
 public:
  virtual ~LinearSolverFactoryGeneric() {}
};
/*
// TODO Check that the matrix structure retained
template <class Scal, class Idx, class Expr>
class Pardiso : public LinearSolver<Scal, Idx, Expr> {
 public:
  using IparmType = std::array<MKL_INT, 64>;

 private:
  template <class T>
  using Field = FieldGeneric<T, Idx>;
  IparmType iparm;

  MKL_INT maxfct, mnum, phase, error, msglvl;
  MKL_INT mtype;               // Matrix type
  MKL_INT nrhs;                // Number of RHS'
  MKL_INT num_eqn;             // Number of equations
  Range<Idx> idx_range;
  std::vector<MKL_INT> mkl_ia; // ia[i] is index of the first nonzero
                               // element in equation i, i=0..num_eqn
  std::vector<MKL_INT> mkl_ja; // ja[j] is the solution vector position
                               // corresponding to a[j], j=0..num_nonzero-1
  std::vector<Scal> mkl_a;     // Matrix content (nonzero elements)
  std::vector<Scal> mkl_rhs;
  double ddum;                 // Dummy double
  MKL_INT idum;                // Dummy int
  std::array<std::int64_t, 64> pt; // Parameters storage
  int num_nonzero;             // Number of nonzero elements
  bool initialized;
  Scal last_hash_;

  void SetIparm(size_t one_based_idx, MKL_INT value) {
    iparm[one_based_idx - 1] = value;
  }

  void SetParameters() {
    iparm.fill(0);
    SetIparm(1, 1);        // No solver default
    SetIparm(4, 0);        // Iterative Krylov Subspace iterations
    SetIparm(10, 13);      // Perturb the pivot elements with 1E-13
    SetIparm(14, 0);       // Output: Number of perturbed pivots
    SetIparm(18, 0);       // Output: Number of nonzeros in the factor LU
    SetIparm(19, 0);       // Output: Mflops for LU factorization

    MKL_INT single_precision;
    if (sizeof(Scal) == 4) {
      single_precision = 1;
    } else if (sizeof(Scal) == 8) {
      single_precision = 0;
    } else {
      assert(false); // sizeof(Scal) should be either 4 or 8")
    }
    SetIparm(28, single_precision);

    SetIparm(35, 1);      // Zero-based indexing
    maxfct = 1;           // Maximum number of numerical factorizations.
    mnum = 1;             // Which factorization to use.

#if __DEBUG > 1
    msglvl = 1;           // Write statistics
    SetIparm(27, 1);      // Check the matrix
#else
    msglvl = 0;           // No statistical information
    SetIparm(27, 0);      // No matrix check
#endif

    error = 0;            // Initialize error flag

    mtype = 2;           // Real symmetric matrix
    nrhs = 1;             // Use one RHS
  }
  void FillIndexArrays(const Field<Expr>& system)
  {
    idx_range = system.GetRange();
    num_eqn = system.size();
    num_nonzero = 0;

    for (size_t i = 0; i < system.size(); ++i) {
      num_nonzero += system[Idx(i)].size();
    }

    mkl_ia.resize(num_eqn + 1);
    mkl_ja.resize(num_nonzero);
    mkl_a.resize(num_nonzero);
    mkl_rhs.resize(num_eqn);

    // matrix structure
    size_t m=0; // index of current nonzero element
    for(size_t i = 0; i < system.size(); ++i) {
      const auto &expr = system[Idx(i)];
      mkl_ia[i] = m;
      for(size_t k = 0; k < expr.size(); ++k) {
        mkl_ja[m] = expr[k].idx.GetRaw();
        ++m;
      }
    }
    mkl_ia[num_eqn] = num_nonzero;
  }
  void FillContent(const Field<Expr>& system)
  {
    size_t m = 0;
    for(size_t i = 0; i < system.size(); ++i) {
      const auto &expr = system[Idx(i)];
      mkl_rhs[i] = -expr.GetConstant();
      for(size_t k = 0; k < expr.size(); ++k) {
        mkl_a[m] = expr[k].coeff;
        ++m;
      }
    }
    assert(m == mkl_a.size());
  }
  Scal CalcMatrixHash() {
    Scal res = 0.;
    for (auto a : mkl_a) {
      res *= -0.9;
      res += a;
    }
    return res;
  }
  void Init(const Field<Expr>& system) {
    if (initialized) {
      assert(system.size() == static_cast<size_t>(num_eqn));
      FillContent(system);
      return;
    }

    FillIndexArrays(system);
    FillContent(system);

    // Reordering and Symbolic Factorization.
    // This step also allocates all memory that is necessary
    // for the factorization.
    phase = 11;
    mkl::PARDISO(pt.data(), &maxfct, &mnum, &mtype, &phase,
            &num_eqn, mkl_a.data(), mkl_ia.data(), mkl_ja.data(),
            &idum, &nrhs, iparm.data(), &msglvl, &ddum, &ddum, &error);
    initialized=true;
  }

 public:
  Pardiso(const std::map<size_t, MKL_INT>& _iparm_one_based, MKL_INT _mtype)
      : initialized(false),
        last_hash_(3.14) {
    SetParameters();
    pt.fill(0);

    if (_mtype != 0) {
      mtype = _mtype;
    }

    for (auto it = _iparm_one_based.begin();
        it != _iparm_one_based.end(); ++it) {
      SetIparm(it->first, it->second);
    }
  }
  Field<Scal> Solve(const Field<Expr>& system) override {
    Init(system);

    Field<Scal> res(idx_range);

    //Numerical factorization
    auto hash = CalcMatrixHash();
    if (hash != last_hash_) {
      last_hash_ = hash;
      phase = 22;
      mkl::PARDISO(
          pt.data(), &maxfct, &mnum, &mtype, &phase,
          &num_eqn, mkl_a.data(), mkl_ia.data(), mkl_ja.data(),
          &idum, &nrhs, iparm.data(), &msglvl, &ddum, &ddum, &error);
    }

    // Back substitution and iterative refinement
    phase = 33;
    mkl::PARDISO(
        pt.data(), &maxfct, &mnum, &mtype, &phase,
        &num_eqn, mkl_a.data(), mkl_ia.data(), mkl_ja.data(),
        &idum, &nrhs, iparm.data(), &msglvl, mkl_rhs.data(),
        res.data(), &error);

    return res;
  }
  ~Pardiso() {
    if(!initialized) {
      return;
    }

    phase = -1;           // Release internal memory.
    mkl::PARDISO(
        pt.data(), &maxfct, &mnum, &mtype, &phase,
        &num_eqn, mkl_a.data(), mkl_ia.data(), mkl_ja.data(),
        &idum, &nrhs, iparm.data(), &msglvl, &ddum, &ddum, &error);
    initialized = false;
  }
};

class PardisoFactory : public LinearSolverFactoryGeneric {
 private:
  std::map<size_t, MKL_INT> iparm_one_based_;
  MKL_INT mtype_;

 public:
  PardisoFactory(const std::map<size_t, MKL_INT>& iparm_one_based,
                 MKL_INT mtype)
      : iparm_one_based_(iparm_one_based),
        mtype_(mtype) {}
  template <class Scal, class Idx, class Expr>
  std::shared_ptr<LinearSolver<Scal, Idx, Expr>> Create() const {
    return std::make_shared<
        Pardiso<Scal, Idx, Expr>>(iparm_one_based_, mtype_);
  }
};*/

template <class Scal, class Idx, class Expr>
class LuDecomposition : public LinearSolver<Scal, Idx, Expr> {
  template <class T>
  using Field = GField<T, Idx>;

 public:
  Field<Scal> Solve(const Field<Expr>& system) override {
    Field<Scal> res(system.GetRange(), 0);

    // forward step
    for(size_t i = 0; i < system.size(); ++i) {
      Scal sum = 0;
      const auto &expr = system[Idx(i)];
      size_t k = 0;
      while (k < expr.size() && expr[k].idx.GetRaw() < i) {
        sum += expr[k].coeff * res[expr[k].idx];
        ++k;
      }
      assert(k < expr.size() && expr[k].idx.GetRaw() == i);
      Scal coeff_diag = expr[k].coeff;
      res[Idx(i)] = (-expr.GetConstant() - sum) / coeff_diag;
    }

    // backward step
    for (size_t i = system.size(); i > 0; ) {
      --i;
      Scal sum = 0;
      const auto &expr = system[Idx(i)];
      size_t k = expr.size() - 1;
      while (k >= 0 && expr[k].idx.GetRaw() > i) {
        sum += expr[k].coeff * res[expr[k].idx];
        --k;
      }
      assert(k >= 0 && expr[k].idx.GetRaw() == i);
      Scal coeff_diag = expr[k].coeff;
      res[Idx(i)] -= sum / coeff_diag;
    }

    return res;
  }
};

class LuDecompositionFactory : public LinearSolverFactoryGeneric {
 private:
 public:
  template <class Scal, class Idx, class Expr>
  std::shared_ptr<LinearSolver<Scal, Idx, Expr>> Create() const {
    return std::make_shared<LuDecomposition<Scal, Idx, Expr>>();
  }
};

template <class Scal, class Idx, class Expr>
class LuDecompositionRelaxed : public LinearSolver<Scal, Idx, Expr> {
  template <class T>
  using Field = GField<T, Idx>;
  Scal tolerance_;
  size_t num_iters_limit_;
  Scal relaxation_factor_;

 public:
  LuDecompositionRelaxed(Scal tolerance, size_t num_iters_limit,
                         Scal relaxation_factor)
      : tolerance_(tolerance),
        num_iters_limit_(num_iters_limit),
        relaxation_factor_(relaxation_factor) {}
  Field<Scal> Solve(const Field<Expr>& system) override {
    auto range = system.GetRange();
    Field<Scal> res(range, 0);

    Field<Scal> corr(range, 0);
    Field<Scal> f(range, 0);
    for (auto idx : range) {
      f[idx] = system[idx].GetConstant();
    }

    size_t iter = 0;
    Scal diff = 0.;
    do {
      diff = 0.;
      // forward step
      for(size_t i = 0; i < system.size(); ++i) {
        Idx idx(i);
        const auto &expr = system[idx];
        Scal sum = 0;
        size_t k = 0;
        while (k < expr.size() && expr[k].idx.GetRaw() < i) {
          sum += expr[k].coeff * corr[expr[k].idx];
          ++k;
        }
        // TODO: measure assertion overhead
        assert(k < expr.size() && expr[k].idx.GetRaw() == i);
        Scal coeff_diag = expr[k].coeff + relaxation_factor_;
        corr[idx] = (-f[idx] - sum) / coeff_diag;
      }

      // backward step
      for (size_t i = system.size(); i > 0; ) {
        --i;
        Idx idx(i);
        const auto &expr = system[idx];
        Scal sum = 0;
        size_t k = expr.size() - 1;
        while (k >= 0 && expr[k].idx.GetRaw() > i) {
          sum += expr[k].coeff * res[expr[k].idx];
          --k;
        }
        assert(k >= 0 && expr[k].idx.GetRaw() == i);
        Scal coeff_diag = expr[k].coeff + relaxation_factor_;
        corr[idx] -= sum / coeff_diag;
      }

      for (auto idx : range) {
        res[idx] += corr[idx];
        diff = std::max(diff, std::abs(corr[idx]));
      }
      for (auto idx : range) {
        f[idx] = system[idx].Evaluate(res);
      }
    } while (diff > tolerance_ && iter++ < num_iters_limit_);

    std::cout << "iter = " << iter << ", diff = " << diff << std::endl;

    return res;
  }
};

class LuDecompositionRelaxedFactory : public LinearSolverFactoryGeneric {
 private:
  double tolerance_;
  size_t num_iters_limit_;
  double relaxation_factor_;
 public:
  LuDecompositionRelaxedFactory(double tolerance, size_t num_iters_limit,
                                double relaxation_factor)
      : tolerance_(tolerance),
        num_iters_limit_(num_iters_limit),
        relaxation_factor_(relaxation_factor) {}
  template <class Scal, class Idx, class Expr>
  std::shared_ptr<LinearSolver<Scal, Idx, Expr>> Create() const {
    return std::make_shared<LuDecompositionRelaxed<Scal, Idx, Expr>>(
        tolerance_, num_iters_limit_, relaxation_factor_);
  }
};

template <class Scal, class Idx, class Expr>
class GaussSeidel : public LinearSolver<Scal, Idx, Expr> {
  template <class T>
  using Field = GField<T, Idx>;
  Scal tolerance_;
  size_t num_iters_limit_;
  Scal relaxation_factor_;

 public:
  GaussSeidel(Scal tolerance, size_t num_iters_limit,
                         Scal relaxation_factor)
      : tolerance_(tolerance),
        num_iters_limit_(num_iters_limit),
        relaxation_factor_(relaxation_factor) {}
  Field<Scal> Solve(const Field<Expr>& system) override {
    Field<Scal> res(system.GetRange(), 0);

    size_t iter = 0;
    Scal diff = 0.;
    do {
      diff = 0.;
      for (size_t raw = 0; raw < system.size(); ++raw) {
        Idx idx(raw);
        const Expr& eqn = system[idx];
        Scal sum = 0.;
        Scal diag_coeff = 0.;
        for (size_t i = 0; i < eqn.size(); ++i) {
          const auto& term = eqn[i];
          if (term.idx != idx) {
            sum += term.coeff * res[term.idx];
          } else {
            diag_coeff = term.coeff;
          }
        }
        Scal value = -(eqn.GetConstant() + sum) / diag_coeff;
        Scal corr = value - res[idx];
        diff = std::max(diff, std::abs(corr));
        res[idx] += corr * relaxation_factor_;
      }
    } while (diff > tolerance_ && iter++ < num_iters_limit_);

    std::cout << "iter = " << iter << ", diff = " << diff << std::endl;

    return res;
  }
};

class GaussSeidelFactory : public LinearSolverFactoryGeneric {
 private:
  double tolerance_;
  size_t num_iters_limit_;
  double relaxation_factor_;
 public:
  GaussSeidelFactory(double tolerance, size_t num_iters_limit,
                     double relaxation_factor)
      : tolerance_(tolerance),
        num_iters_limit_(num_iters_limit),
        relaxation_factor_(relaxation_factor) {}
  template <class Scal, class Idx, class Expr>
  std::shared_ptr<LinearSolver<Scal, Idx, Expr>> Create() const {
    return std::make_shared<GaussSeidel<Scal, Idx, Expr>>(
        tolerance_, num_iters_limit_, relaxation_factor_);
  }
};

template <class Scal, class Idx, class Expr>
class Jacobi : public LinearSolver<Scal, Idx, Expr> {
  template <class T>
  using Field = GField<T, Idx>;
  Scal tolerance_;
  size_t num_iters_limit_;
  Scal relaxation_factor_;

 public:
  Jacobi(Scal tolerance, size_t num_iters_limit,
                         Scal relaxation_factor)
      : tolerance_(tolerance),
        num_iters_limit_(num_iters_limit),
        relaxation_factor_(relaxation_factor) {}
  Field<Scal> Solve(const Field<Expr>& system) override {
    Field<Scal> res(system.GetRange(), 0);
    auto next = res;

    size_t iter = 0;
    Scal diff = 0.;
    do {
      diff = 0.;
      for (size_t raw = 0; raw < system.size(); ++raw) {
        Idx idx(raw);
        const Expr& eqn = system[idx];
        Scal sum = 0.;
        Scal diag_coeff = 0.;
        for (size_t i = 0; i < eqn.size(); ++i) {
          const auto& term = eqn[i];
          if (term.idx != idx) {
            sum += term.coeff * res[term.idx];
          } else {
            diag_coeff = term.coeff;
          }
        }
        Scal value = -(eqn.GetConstant() + sum) / diag_coeff;
        Scal corr = value - res[idx];
        diff = std::max(diff, std::abs(corr));
        next[idx] = res[idx] + corr * relaxation_factor_;
      }
      std::swap(res, next);
    } while (diff > tolerance_ && iter++ < num_iters_limit_);

    std::cout << "iter = " << iter << ", diff = " << diff << std::endl;

    return res;
  }
};

class JacobiFactory : public LinearSolverFactoryGeneric {
 private:
  double tolerance_;
  size_t num_iters_limit_;
  double relaxation_factor_;
 public:
  JacobiFactory(double tolerance, size_t num_iters_limit,
                     double relaxation_factor)
      : tolerance_(tolerance),
        num_iters_limit_(num_iters_limit),
        relaxation_factor_(relaxation_factor) {}
  template <class Scal, class Idx, class Expr>
  std::shared_ptr<LinearSolver<Scal, Idx, Expr>> Create() const {
    return std::make_shared<Jacobi<Scal, Idx, Expr>>(
        tolerance_, num_iters_limit_, relaxation_factor_);
  }
};

class LinearSolverFactory {
  std::shared_ptr<const LinearSolverFactoryGeneric> p_generic_factory_;
  template <class Factory, class Scal, class Idx, class Expr>
  bool TryCreate(std::shared_ptr<LinearSolver<Scal, Idx, Expr>>& res) const {
    if (auto p_factory_ =
        dynamic_cast<const Factory*>(p_generic_factory_.get())) {
      res = p_factory_->template Create<Scal, Idx, Expr>();
      return true;
    }
    return false;
  }
 public:
  template <class Factory>
  explicit LinearSolverFactory(
      std::shared_ptr<const Factory> p_generic_factory)
      : p_generic_factory_(p_generic_factory) {}
  template <class Scal, class Idx, class Expr>
  std::shared_ptr<LinearSolver<Scal, Idx, Expr>> Create() const {
    std::shared_ptr<LinearSolver<Scal, Idx, Expr>> res;
    bool found = false;

    /*found = found ||
        TryCreate<PardisoFactory, Scal, Idx, Expr>(res);*/
    found = found ||
        TryCreate<LuDecompositionFactory, Scal, Idx, Expr>(res);
    found = found ||
        TryCreate<LuDecompositionRelaxedFactory, Scal, Idx, Expr>(res);
    found = found ||
        TryCreate<GaussSeidelFactory, Scal, Idx, Expr>(res);
    found = found ||
        TryCreate<JacobiFactory, Scal, Idx, Expr>(res);

    if (!found) {
      throw std::runtime_error(
        "LinearSolverFactory: Create() is undefined");
    }
    return res;
  }
};

template <class M, class Expr, class Scal=typename M::Scal>
typename M::LS ConvertLs(
    const FieldCell<Expr>& fce,
    std::vector<Scal>& la, 
    std::vector<Scal>& lb, 
    std::vector<Scal>& lx, 
    const M& m) {
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
  lb.resize(n, 1.);
  lx.resize(n, 0.);

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
              << " l.st[j]=" << MIdx(l.st[j]) 
              << std::endl;
          assert(false);
        }
        la[i] = e[j].coeff;
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

} // namespace solver
