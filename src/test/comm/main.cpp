// Created by Petr Karnakov on 30.05.2018
// Copyright 2018 ETH Zurich

#undef NDEBUG
#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

#include "distr/distrsolver.h"
#include "geom/mesh.h"
#include "kernel/kernelmeshpar.h"
#include "solver/pois.h"
#include "solver/solver.h"
#include "util/format.h"
#include "util/mpi.h"
#include "util/suspender.h"

template <class T>
std::ostream& operator<<(std::ostream& o, const std::vector<T>& v) {
  std::string p = "";
  for (auto a : v) {
    o << p << a;
    p = " ";
  }
  return o;
}

struct GPar {};

template <class M_>
class Test : public KernelMeshPar<M_, GPar> {
 public:
  using P = KernelMeshPar<M_, GPar>;
  using M = M_;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using MIdx = typename M::MIdx;
  using Par = GPar;
  static constexpr size_t dim = M::dim;

  using P::P;
  void Run() override;

 protected:
  using P::m;
  using P::var;

 private:
  void TestComm();
  void TestReduce();
  void TestScatter();
  void TestShared();

  bool Cmp(Scal a, Scal b) {
    return std::abs(a - b) < 1e-10;
  }
  bool Cmp(Vect a, Vect b) {
    return a.dist(b) < 1e-10;
  }
  bool Cmp(MIdx a, MIdx b) {
    return a == b;
  }
  template <class T>
  bool Cmp(const std::vector<T>& a, const std::vector<T>& b) {
    return a == b;
  }

  FieldCell<Scal> fc_;
  FieldCell<Vect> fcv_;
  Scal r_; // test Reduce
  std::pair<Scal, int> rsi_; // test Reduce minloc
  std::vector<Scal> rvs_; // reduction vector<Scal> (concatenation)
  std::vector<int> rvi_; // reduction vector<int> (concatenation)
  std::vector<std::vector<int>> rvvi_; // reduction vector<vector<int>>
  std::vector<std::vector<Scal>> rvvs_; // reduction vector<vector<int>>
  static std::ofstream fout;
};

template <class M>
std::ofstream Test<M>::fout;

template <class Idx, class M>
typename M::Scal DiffMax(
    const GField<typename M::Scal, Idx>& u,
    const GField<typename M::Scal, Idx>& v, const M& m) {
  using Scal = typename M::Scal;
  Scal r = 0;
  for (auto i : m.template GetRangeIn<Idx>()) {
    r = std::max(r, std::abs(u[i] - v[i]));
  }
  return r;
}

template <class Idx, class M>
typename M::Scal DiffMean(
    const GField<typename M::Scal, Idx>& u,
    const GField<typename M::Scal, Idx>& v, const M& m) {
  using Scal = typename M::Scal;
  Scal r = 0;
  Scal w = 0;
  for (auto i : m.template GetRangeIn<Idx>()) {
    r += std::abs(u[i] - v[i]);
    w += 1.;
  }
  return r / w;
}

template <class Idx, class M>
typename M::Scal Max(const GField<typename M::Scal, Idx>& u, const M& m) {
  using Scal = typename M::Scal;
  Scal r = 0;
  for (auto i : m.template GetRangeIn<Idx>()) {
    r = std::max(r, u[i]);
  }
  return r;
}

template <class Idx, class M>
typename M::Scal Mean(const GField<typename M::Scal, Idx>& u, const M& m) {
  using Scal = typename M::Scal;
  Scal r = 0;
  Scal w = 0.;
  for (auto i : m.template GetRangeIn<Idx>()) {
    r += u[i];
    w += 1.;
  }
  return r / w;
}

#define CMP(a, b)   \
  if (!Cmp(a, b)) { \
    fassert(false); \
  }

// Print CMP
#define PCMP(a, b)                                                          \
  {                                                                         \
    auto _fmt = fout.flags();                                               \
    fout << std::scientific << std::setprecision(16) << "CMP " << #a << "=" \
         << a << "\n    " << #b << "=" << b << std::endl;                   \
    CMP(a, b);                                                              \
    fout.flags(_fmt);                                                       \
  }

// Print CMP standard flags
#define PCMPF(a, b)                                             \
  fout << "CMP " << std::setw(15) << #a << "=" << a << "\n    " \
       << std::setw(15) << #b << "=" << b << std::endl;         \
  CMP(a, b);

template <class M>
void Test<M>::TestComm() {
  auto sem = m.GetSem("Comm");
  auto func = [&](Vect x) {
    auto gl = m.GetGlobalLength();
    x = (x + 16 * gl) / gl;
    x = (x - Vect(MIdx(x))) * gl;
    return std::sin(x[0]) * std::cos(x[1]) * std::exp(x[2]);
  };
  auto funcv = [=](Vect v) {
    return Vect(func(v), func(v * 2.), func(v * 3.));
  };
  auto& bc = m.GetIndexCells();
  if (sem("init")) {
    fc_.Reinit(m);
    fcv_.Reinit(m);
    for (auto c : m.Cells()) {
      auto x = m.GetCenter(c);
      fc_[c] = func(x);
      fcv_[c] = funcv(x);
    }
    m.Comm(&fc_);
    m.Comm(&fcv_);
  }
  if (sem("check")) {
    for (auto c : m.AllCells()) {
      auto x = m.GetCenter(c);
      if (!Cmp(fc_[c], func(x))) {
        fout << bc.GetMIdx(c) << " " << fc_[c] << " != " << func(x) << " "
             << std::endl;
        fassert(false);
      }
      if (!Cmp(fcv_[c], funcv(x))) {
        fout << bc.GetMIdx(c) << " " << fcv_[c] << " != " << funcv(x) << " "
             << std::endl;
        fassert(false);
      }
    }
  }
}

template <class M>
void Test<M>::TestReduce() {
  auto sem = m.GetSem("Reduce");
  MIdx nprocs(var.Int["px"], var.Int["py"], var.Int["pz"]);
  MIdx nblocks(var.Int["bx"], var.Int["by"], var.Int["bz"]);
  GBlock<IdxCell, dim> blocks_all(nprocs * nblocks);
  GIndex<size_t, dim> indexc(nprocs * nblocks);
  GBlock<IdxCell, dim> procs(nprocs);
  GBlock<IdxCell, dim> blocks(nblocks);
  const MIdx block_index =
      m.GetInBlockCells().GetBegin() / m.GetInBlockCells().GetSize();
  auto func = [](MIdx w) {
    return std::sin(w[0] + 1.7) * std::cos(w[1]) * std::exp(w[2] * 0.1);
  };
  if (sem("sum")) {
    r_ = func(block_index);
    m.Reduce(&r_, Reduction::sum);
  }
  if (sem("sum-check")) {
    Scal exact = 0.;
    for (auto b : blocks_all) {
      exact += func(b);
    }
    PCMP(r_, exact);
  }
  if (sem("prod")) {
    r_ = func(block_index);
    m.Reduce(&r_, Reduction::prod);
  }
  if (sem("prod-check")) {
    Scal exact = 1.;
    for (auto b : blocks_all) {
      exact *= func(b);
    }
    PCMP(r_, exact);
  }
  if (sem("max")) {
    r_ = func(block_index);
    m.Reduce(&r_, Reduction::max);
  }
  if (sem("max-check")) {
    Scal s = -std::numeric_limits<Scal>::max();
    for (auto b : blocks_all) {
      s = std::max(s, func(b));
    }
    PCMP(r_, s);
  }
  if (sem("min")) {
    r_ = func(block_index);
    m.Reduce(&r_, Reduction::min);
  }
  if (sem("min-check")) {
    Scal s = std::numeric_limits<Scal>::max();
    for (auto b : blocks_all) {
      s = std::min(s, func(b));
    }
    PCMP(r_, s);
  }
  if (sem("minloc")) {
    MIdx w(block_index);
    rsi_ = std::make_pair(func(w), indexc.GetIdx(w));
    m.Reduce(&rsi_, Reduction::minloc);
  }
  if (sem("minloc-check")) {
    Scal s = std::numeric_limits<Scal>::max();
    MIdx ws;
    for (auto w : blocks_all) {
      if (func(w) < s) {
        ws = w;
        s = func(w);
      }
    }
    PCMP(rsi_.first, s);
    PCMP(indexc.GetMIdx(rsi_.second), ws);
  }
  if (sem("maxloc")) {
    MIdx w(block_index);
    rsi_ = std::make_pair(func(w), indexc.GetIdx(w));
    m.Reduce(&rsi_, Reduction::maxloc);
  }
  if (sem("maxloc-check")) {
    Scal s = -std::numeric_limits<Scal>::max();
    MIdx ws;
    for (auto w : blocks_all) {
      if (func(w) > s) {
        ws = w;
        s = func(w);
      }
    }
    PCMP(rsi_.first, s);
    PCMP(indexc.GetMIdx(rsi_.second), ws);
  }
  if (sem("cat")) {
    MIdx w(block_index);
    rvs_.resize(0);
    rvs_.push_back(indexc.GetIdx(w));
    m.Reduce(&rvs_, Reduction::concat);
  }
  if (sem("cat-check")) {
    std::vector<Scal> exact;
    for (auto wp : procs) {
      for (auto wb : blocks) {
        auto w = nblocks * wp + wb;
        exact.push_back(indexc.GetIdx(w));
      }
    }
    if (m.IsRoot()) {
      std::sort(rvs_.begin(), rvs_.end());
      std::sort(exact.begin(), exact.end());
      PCMPF(rvs_, exact);
    }
  }
  const size_t q = 10;
  if (sem("cati")) {
    MIdx w(block_index);
    rvi_.resize(0);
    size_t i = indexc.GetIdx(w);
    for (size_t j = 0; j < i % q; ++j) {
      rvi_.push_back(q * i + j);
    }
    m.Reduce(&rvi_, Reduction::concat);
  }
  if (sem("cati-check")) {
    std::vector<int> rvi;
    for (auto wp : procs) {
      for (auto wb : blocks) {
        auto w = nblocks * wp + wb;
        size_t i = indexc.GetIdx(w);
        for (size_t j = 0; j < i % q; ++j) {
          rvi.push_back(q * i + j);
        }
      }
    }
    if (m.IsRoot()) {
      std::sort(rvi_.begin(), rvi_.end());
      std::sort(rvi.begin(), rvi.end());
      PCMP(rvi_, rvi);
    }
  }
  if (sem("catvi")) {
    MIdx w(block_index);
    rvvi_.resize(0);
    size_t i = indexc.GetIdx(w);
    for (size_t j = 0; j < i % q; ++j) {
      rvvi_.push_back(std::vector<int>({int(q * i + j)}));
    }
    m.Reduce(&rvvi_, Reduction::concat);
  }
  if (sem("catvi-check")) {
    std::vector<std::vector<int>> rvvi;
    for (auto wp : procs) {
      for (auto wb : blocks) {
        auto w = nblocks * wp + wb;
        size_t i = indexc.GetIdx(w);
        for (size_t j = 0; j < i % q; ++j) {
          rvvi.push_back(std::vector<int>({int(q * i + j)}));
        }
      }
    }
    if (m.IsRoot()) {
      std::sort(rvvi_.begin(), rvvi_.end());
      std::sort(rvvi.begin(), rvvi.end());
      PCMP(rvvi_, rvvi);
    }
  }
  if (sem("bcast-catvi")) {
    MIdx w(block_index);
    rvvi_.resize(0);
    size_t i = indexc.GetIdx(block_index);
    for (size_t j = 0; j < (i + 5) % q; ++j) {
      rvvi_.push_back(std::vector<int>({int(100 * i + j)}));
    }
    m.Bcast(&rvvi_);
  }
  if (sem("bcast-catvi-check")) {
    // init for root nblocks
    std::vector<std::vector<int>> rvvi;
    size_t i = indexc.GetIdx(MIdx(0));
    for (size_t j = 0; j < (i + 5) % q; ++j) {
      rvvi.push_back(std::vector<int>({int(100 * i + j)}));
    }
    PCMP(rvvi_, rvvi);
  }
}

template <class M>
void Test<M>::TestScatter() {
  auto sem = m.GetSem("scatter");
  MIdx nprocs(var.Int["px"], var.Int["py"], var.Int["pz"]);
  MIdx nblocks(var.Int["bx"], var.Int["by"], var.Int["bz"]);
  GBlock<IdxCell, dim> blocks_all(nprocs * nblocks);
  GIndex<size_t, dim> indexc(nprocs * nblocks);
  GBlock<IdxCell, dim> procs(nprocs);
  GBlock<IdxCell, dim> blocks(nblocks);
  const MIdx block_index =
      m.GetInBlockCells().GetBegin() / m.GetInBlockCells().GetSize();
  auto GetBlockData = [](size_t i) {
    std::vector<Scal> r;
    r.push_back(Scal(i));
    for (size_t j = 0; j < i; ++j) {
      r.push_back(Scal(j));
    }
    return r;
  };
  if (sem("vector<Scal>")) {
    rvvs_.clear();
    for (auto wp : procs) {
      for (auto wb : blocks) {
        auto w = nblocks * wp + wb;
        size_t i = indexc.GetIdx(w);
        rvvs_.push_back(GetBlockData(i));
      }
    }
    if (m.IsRoot()) {
      fout << "rvvs=" << rvvs_ << std::endl;
      m.Scatter({&rvvs_, &rvs_});
    } else {
      m.Scatter({nullptr, &rvs_});
    }
  }
  if (sem("check")) {
    assert(rvs_.size() > 0);
    PCMPF(rvs_, GetBlockData(std::lround(rvs_[0])));
  }
  if (sem("gather")) {
    MIdx w(block_index);
    size_t i = indexc.GetIdx(w);
    rvvs_ = {GetBlockData(i)};
    m.Reduce(&rvvs_, Reduction::concat);
  }
  if (sem("scatter")) {
    if (m.IsRoot()) {
      fout << "rvvs=" << rvvs_ << std::endl;
      m.Scatter({&rvvs_, &rvs_});
    } else {
      m.Scatter({nullptr, &rvs_});
    }
  }
  if (sem("check")) {
    MIdx w(block_index);
    size_t i = indexc.GetIdx(block_index);
    PCMPF(rvs_, GetBlockData(i));
  }
}

template <class M>
void Test<M>::TestShared() {
  auto sem = m.GetSem("shared");
  struct {
    FieldCell<Scal> fc;
  } * ctx(sem);
  auto& t = *ctx;
  auto& ms = m.GetShared();
  if (sem()) {
    if (m.IsRoot()) {
      fout << "\nTestShared\n";
      fout << util::Format(
          "m.size={} ms.size={}\n", m.GetInBlockCells().GetSize(),
          ms.GetInBlockCells().GetSize());
    }
  }
  // Prints slice of field
  auto print = [&](const FieldCell<Scal>& fc) {
    MIdx last_w(ms.GetAllBlockCells().GetBegin());
    for (auto c : ms.AllCellsM()) {
      const MIdx w(c);
      if (w[2] != 0) {
        continue;
      }
      if (w[1] != last_w[1]) {
        fout << '\n';
      }
      fout << fc[c] << ' ';
      last_w = w;
    }
    fout << '\n';
  };
  if (sem() && m.IsLead()) {
    const auto rank = MpiWrapper(m.GetMpiComm()).GetCommRank();
    t.fc.Reinit(ms, 0);
    for (auto c : ms.CellsM()) {
      t.fc[c] = (rank + 1) * 10 + IdxCell(c).raw() % 10;
    }
    fout << "before Comm\n";
    print(t.fc);
    ms.Comm(&t.fc);
  }
  if (sem() && m.IsLead()) {
    fout << "after Comm\n";
    print(t.fc);
  }
}

template <class M>
void Test<M>::Run() {
  auto sem = m.GetSem("Run");
  if (sem()) {
    if (m.IsLead()) {
      fout.open(
          util::Format("out_{}", MpiWrapper(m.GetMpiComm()).GetCommRank()));
    }
  }
  if (sem.Nested()) {
    TestComm();
  }
  if (sem.Nested()) {
    TestReduce();
  }
  if (sem.Nested()) {
    TestScatter();
  }
  if (sem.Nested()) {
    TestShared();
  }
}

void Main(MPI_Comm comm, Vars& var) {
  using M = MeshCartesian<double, 3>;
  using K = Test<M>;
  using Par = typename K::Par;
  Par par;

  DistrSolver<M, K> ds(comm, var, par);
  ds.Run();
}

int main(int argc, const char** argv) {
  return RunMpi(argc, argv, Main);
}
