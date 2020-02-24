// Created by Petr Karnakov on 25.04.2018
// Copyright 2018 ETH Zurich

// vim: expandtab:smarttab:sw=2:ts=2
#pragma once

#include <array>
#include <cassert>
#include <limits>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

#include "blockface.h"
#include "distr/reduce.h"
#include "field.h"
#include "idx.h"
#include "map.h"
#include "range.h"
#include "rangein.h"
#include "transform.h"
#include "util/histogram.h"
#include "util/suspender.h"
#include "vect.h"

// Returns column of cells cmm,cm,cp,cpp.
// nci: 0 or 1 such that m.GetCell(f, nci) == cp
template <class M>
void GetCellColumn(
    const M& m, IdxFace f, size_t nci, IdxCell& cmm, IdxCell& cm, IdxCell& cp,
    IdxCell& cpp) {
  const size_t d{m.GetIndexFaces().GetDir(f)};
  cp = m.GetCell(f, nci);
  cm = m.GetCell(f, 1 - nci);
  cpp = m.GetCell(cp, d * 2 + nci);
  cmm = m.GetCell(cm, d * 2 + 1 - nci);
}

// TODO: Neighbour faces iterator introducing (cell, face) pairs
// TODO: consider computing some on-the-fly to reduce memory access
template <class Scal_, size_t dim_>
class MeshStructured {
 public:
  using Scal = Scal_;
  static constexpr size_t dim = dim_;
  using Vect = generic::Vect<Scal, dim>;
  using Dir = GDir<dim>;
  using MIdx = GMIdx<dim>;
  using BlockCells = GBlockCells<dim>;
  using BlockFaces = GBlockFaces<dim>;
  using BlockNodes = GBlockNodes<dim>;
  using IndexCells = GIndex<IdxCell, dim>;
  using IndexFaces = GIndex<IdxFace, dim>;
  using IndexNodes = GIndex<IdxNode, dim>;
  using M = MeshStructured;
  static constexpr size_t kCellNumNeighbourFaces = 6;
  static constexpr size_t kCellNumNeighbourNodes = 8;
  static constexpr size_t kFaceNumNeighbourNodes = 4;
  static constexpr size_t kFaceNumNeighbourCells = 2;

  // b: begin, lower corner cell index
  // cs: inner cells size
  // dom: domain, rectangle covering inner cells
  // hl: halo cells from each side
  // isroot: root block
  // gs: global mesh size
  // id: unique id
  MeshStructured(
      MIdx b, MIdx cs, Rect<Vect> dom, int hl, bool isroot, bool islead,
      MIdx gs, int id);
  MeshStructured(const MeshStructured&) = delete;
  MeshStructured(MeshStructured&&) = default;
  MeshStructured& operator=(const MeshStructured&) = delete;
  MeshStructured& operator=(MeshStructured&&) = delete;
  MIdx GetGlobalSize() const {
    return gs_;
  }
  int GetId() const {
    return id_;
  }
  Vect GetGlobalLength() const {
    return gl_;
  }
  Vect GetCellSize() const {
    return h_;
  }
  IntIdx GetHash(MIdx w) {
    // XXX: adhoc, hash for cell index, assume mesh size <= mn
    const size_t mn = 1000;
    return (w[2] * mn + w[1]) * mn + w[0];
  }
  IntIdx GetHash(IdxCell c) {
    return GetHash(GetIndexCells().GetMIdx(c));
  }
  // Indexer
  const IndexCells& GetIndexCells() const {
    return bcr_;
  }
  const IndexFaces& GetIndexFaces() const {
    return bfr_;
  }
  const IndexNodes& GetIndexNodes() const {
    return bnr_;
  }

  // In blocks
  const BlockCells& GetInBlockCells() const {
    return bci_;
  }
  const BlockFaces& GetInBlockFaces() const {
    return bfi_;
  }
  const BlockNodes& GetInBlockNodes() const {
    return bni_;
  }

  // Su blocks
  const BlockCells& GetSuBlockCells() const {
    return bcs_;
  }
  const BlockFaces& GetSuBlockFaces() const {
    return bfs_;
  }
  const BlockNodes& GetSuBlockNodes() const {
    return bns_;
  }

  // All blocks
  const BlockCells& GetAllBlockCells() const {
    return bca_;
  }
  const BlockFaces& GetAllBlockFaces() const {
    return bfa_;
  }
  const BlockNodes& GetAllBlockNodes() const {
    return bna_;
  }
  Vect GetCenter(IdxCell c) const {
    return (dom_.lb + hh_) + Vect(bcr_.GetMIdx(c) - mb_) * h_;
  }
  Vect GetCenter(IdxFace f) const {
    auto p = bfr_.GetMIdxDir(f);
    const MIdx& w = p.first;
    size_t d(p.second);
    Vect r = (dom_.lb + hh_) + Vect(w - mb_) * h_;
    r[d] -= hh_[d];
    return r;
  }
  const Vect& GetSurface(IdxFace f) const {
    return vs_[size_t(bfr_.GetDir(f))];
  }
  Vect GetNode(IdxNode n) const {
    return dom_.lb + Vect(bnr_.GetMIdx(n) - mb_) * h_;
  }
  Scal GetVolume(IdxCell) const {
    return vol_;
  }
  Scal GetArea(IdxFace f) const {
    return va_[size_t(bfr_.GetDir(f))];
  }
  Scal GetAreaFraction(IdxFace) const {
    return 1;
  }
  Scal GetVolumeFraction(IdxCell) const {
    return 1;
  }
  IdxCell GetCell(IdxCell c, size_t q) const {
    assert(q < kCellNumNeighbourFaces);
    return IdxCell(size_t(c) + cnc_[q]);
  }
  IdxFace GetFace(IdxCell c, size_t q) const {
    assert(q < kCellNumNeighbourFaces);
    return IdxFace(size_t(c) + cnf_[q]);
  }
  Scal GetOutwardFactor(IdxCell, size_t q) const {
    assert(q < kCellNumNeighbourFaces);
    return co_[q];
  }
  Vect GetOutwardSurface(IdxCell c, size_t n) const {
    assert(n < kCellNumNeighbourFaces);
    return GetSurface(GetFace(c, n)) * GetOutwardFactor(c, n);
  }
  IdxNode GetNode(IdxCell c, size_t q) const {
    assert(q < kCellNumNeighbourNodes);
    return IdxNode(size_t(c) + cnn_[q]);
  }
  Dir GetDir(IdxFace f) const {
    return bfr_.GetDir(f);
  }
  IdxCell GetCell(IdxFace f, size_t q) const {
    auto qm = kFaceNumNeighbourCells;
    assert(q < qm);
    size_t d(bfr_.GetDir(f));
    return IdxCell(size_t(f) + fnc_[d * qm + q]);
  }
  Vect GetVectToCell(IdxFace f, size_t n) const {
    assert(n < kFaceNumNeighbourCells);
    return GetCenter(GetCell(f, n)) - GetCenter(f);
  }
  IdxNode GetNode(IdxFace f, size_t q) const {
    auto qm = kFaceNumNeighbourNodes;
    assert(q < qm);
    size_t d(bfr_.GetDir(f));
    return IdxNode(size_t(f) + fnn_[d * qm + q]);
  }
  size_t GetNumFaces(IdxCell) const {
    return kCellNumNeighbourFaces;
  }
  // Neighbour cell indices
  GRange<size_t> Nci(IdxCell c) const {
    return GRange<size_t>(0, GetNumFaces(c));
  }
  // sw: stencil half-width, results in stencil [-sw,sw]
  template <size_t sw>
  TransformIterator<IdxCell, GBlock<size_t, dim>> StencilGeneral(
      IdxCell c) const {
    constexpr size_t sn = sw * 2 + 1;
    const GBlock<size_t, dim> bo(MIdx(-sw), MIdx(sn));
    return MakeTransformIterator<IdxCell>(bo, [this, c](MIdx wo) {
      auto& bc = GetIndexCells();
      return bc.GetIdx(bc.GetMIdx(c) + wo);
    });
  }
  // Cell indices over 3x3x3 stencil
  auto Stencil(IdxCell c) const {
    return StencilGeneral<1>(c);
  }
  // Cell indices over 5x5x5 stencil
  auto Stencil5(IdxCell c) const {
    return StencilGeneral<2>(c);
  }
  // Returns id of cell adjacent to c by face f.
  // -1 if f and c are not neighbours
  size_t GetNci(IdxCell c, IdxFace f) const {
    for (size_t q : Nci(c)) {
      if (GetFace(c, q) == f) {
        return q;
      }
    }
    return size_t(-1);
  }
  // Returns id of face opposite to q.
  // XXX: assumes indices of neighbours are -x,+x,-y,+y,-z,+z
  size_t GetOpposite(size_t q) const {
    if (q == size_t(-1)) {
      return q;
    }

    if (q % 2 == 0) {
      return ++q;
    }
    return --q;
  }
  size_t GetNumNodes(IdxCell) const {
    return kCellNumNeighbourNodes;
  }
  size_t GetNumCells(IdxFace) const {
    return kFaceNumNeighbourCells;
  }
  size_t GetNumNodes(IdxFace) const {
    return kFaceNumNeighbourNodes;
  }
  // inner means not halo
  bool IsInner(IdxCell c) const {
    MIdx w = bcr_.GetMIdx(c);
    return bci_.GetBegin() <= w && w < bci_.GetEnd();
  }

  GIndex<IdxCell, dim> GetIndexer(IdxCell*) const {
    return bcr_;
  }
  GIndex<IdxFace, dim> GetIndexer(IdxFace*) const {
    return bfr_;
  }
  GIndex<IdxNode, dim> GetIndexer(IdxNode*) const {
    return bnr_;
  }
  template <class Idx>
  GIndex<Idx, dim> GetIndexer() const {
    return GetIndexer((Idx*)nullptr);
  }

  // Type-cast to GRange required for GField initialization
  template <class Idx>
  operator GRange<Idx>() const {
    return GRange<Idx>(Idx(GetIndexer<Idx>().size()));
  }

  // Returns range of inner indices
  template <class Idx>
  GRangeIn<Idx, dim> GetIn() const {
    return GetIn((Idx*)0);
  }
  GRangeIn<IdxCell, dim> GetIn(IdxCell*) const {
    return GRangeIn<IdxCell, dim>(GetIndexCells(), GetInBlockCells());
  }
  GRangeIn<IdxFace, dim> GetIn(IdxFace*) const {
    return GRangeIn<IdxFace, dim>(GetIndexFaces(), GetInBlockFaces());
  }
  GRangeIn<IdxNode, dim> GetIn(IdxNode*) const {
    return GRangeIn<IdxNode, dim>(GetIndexNodes(), GetInBlockNodes());
  }
  GRangeIn<IdxCell, dim> Cells() const {
    return GetIn<IdxCell>();
  }
  GRangeIn<IdxFace, dim> Faces() const {
    return GetIn<IdxFace>();
  }
  GRangeIn<IdxNode, dim> Nodes() const {
    return GetIn<IdxNode>();
  }

  // Returns range of support indices
  template <class Idx>
  GRangeIn<Idx, dim> GetSu() const {
    return GetSu((Idx*)0);
  }
  GRangeIn<IdxCell, dim> GetSu(IdxCell*) const {
    return GRangeIn<IdxCell, dim>(GetIndexCells(), GetSuBlockCells());
  }
  GRangeIn<IdxFace, dim> GetSu(IdxFace*) const {
    return GRangeIn<IdxFace, dim>(GetIndexFaces(), GetSuBlockFaces());
  }
  GRangeIn<IdxNode, dim> GetSu(IdxNode*) const {
    return GRangeIn<IdxNode, dim>(GetIndexNodes(), GetSuBlockNodes());
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

  // Returns range of all indices
  template <class Idx>
  GRangeIn<Idx, dim> GetAll() const {
    return GetAll((Idx*)0);
  }
  GRangeIn<IdxCell, dim> GetAll(IdxCell*) const {
    return GRangeIn<IdxCell, dim>(GetIndexCells(), GetAllBlockCells());
  }
  GRangeIn<IdxFace, dim> GetAll(IdxFace*) const {
    return GRangeIn<IdxFace, dim>(GetIndexFaces(), GetAllBlockFaces());
  }
  GRangeIn<IdxNode, dim> GetAll(IdxNode*) const {
    return GRangeIn<IdxNode, dim>(GetIndexNodes(), GetAllBlockNodes());
  }
  GRangeIn<IdxCell, dim> AllCells() const {
    return GetAll<IdxCell>();
  }
  GRangeIn<IdxFace, dim> AllFaces() const {
    return GetAll<IdxFace>();
  }
  GRangeIn<IdxNode, dim> AllNodes() const {
    return GetAll<IdxNode>();
  }
  // Wrappers to satisfy the interface of Embed<M>
  Vect GetFaceCenter(IdxFace f) const {
    return GetCenter(f);
  }
  Vect GetFaceCenter(IdxCell) const {
    return GetNan<Vect>();
  }
  auto CFaces() const {
    return GRange<IdxCell>();
  }
  auto SuCFaces() const {
    return GRange<IdxCell>();
  }

  bool IsInside(IdxCell c, Vect vect) const {
    for (auto q : Nci()) {
      IdxFace f = GetFace(c, q);
      if (GetOutwardSurface(c, q).dot(vect - GetCenter(f)) > 0) {
        return false;
      }
    }
    return true;
  }
  IdxCell FindNearestCell(Vect x) const {
    Scal dn = std::numeric_limits<Scal>::max(); // sqrdist to nearest
    IdxCell cn(0); // cell nearest
    for (auto c : Cells()) {
      Scal d = x.sqrdist(GetCenter(c));
      if (d < dn) {
        dn = d;
        cn = c;
      }
    }
    return cn;
  }
  Vect GetNormal(IdxFace f) const {
    return GetSurface(f) / GetArea(f);
  }
  bool IsRoot() const {
    return isroot_;
  }
  bool IsLead() const {
    return islead_;
  }

  // CheckNan flag
  bool CN() const {
    return checknan_;
  }
  void SetCN(bool c) {
    checknan_ = c;
  }
  // Pairs face,nci for which the halos cells
  // are set to nan after each communication
  const std::vector<std::pair<IdxFace, size_t>>& GetNanFaces() const {
    return vfnan_;
  }
  void SetNanFaces(const std::vector<std::pair<IdxFace, size_t>>& vfnan) {
    vfnan_ = vfnan;
  }
  // Maximum number of communication requests (limited in legacy Cubism)
  size_t GetMaxComm() const {
    return maxcomm_;
  }
  void SetMaxComm(size_t maxcomm) {
    maxcomm_ = maxcomm;
  }
  // time sampler
  Sampler& GetSampler() {
    return samp_;
  }
  const Sampler& GetSampler() const {
    return samp_;
  }
  void SeedSample() {
    samp_.SeedSample();
  }
  void CollectSample(const std::string& name) {
    samp_.CollectSample(name);
  }
  // Fills halo cell with garbage.
  // Using actual NaNs not allowed since some code relies on u*0 == 0
  template <class T>
  void ApplyNanFaces(FieldCell<T>& fc) {
    for (auto p : vfnan_) {
      IdxFace f = p.first;
      size_t nci = p.second;
      IdxCell cmm, cm, cp, cpp;
      GetCellColumn(*this, f, nci, cmm, cm, cp, cpp);
      fc[cm] = T(1e10);
      fc[cmm] = T(1e10);
    }
  }

  // Effective dimension
  size_t GetEdim() const {
    return edim_;
  }
  void SetEdim(size_t edim) {
    edim_ = edim;
  }
  MeshStructured& GetMesh() {
    return *this;
  }
  const MeshStructured& GetMesh() const {
    return *this;
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
    std::string prefix = ""; // custom prefix for config (tol, solver, maxiter)
  };
  Sem GetSem(std::string name = "") {
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

  // Comm request
  struct Co {
    virtual ~Co() {}
    // Number of scalar cell fields  (used in Dump and ghost communication)
    virtual size_t GetSize() const = 0;
    // Currently active index and element stride required for dumping individual
    // AoS components in vector fields
    virtual int GetIndex() const = 0;
    virtual int GetStride() const = 0;
    // Pointer to the start of data.  Retain type-safety of data type, requires
    // that Vect::value_type is Scal and must not be changed
    virtual Scal* GetBasePtr() = 0;
  };
  // FieldCell<Scal>
  struct CoFcs : public Co {
    CoFcs(FieldCell<Scal>* f) : f(f) {}
    size_t GetSize() const override {
      return 1;
    }
    int GetIndex() const override {
      return 0;
    }
    int GetStride() const override {
      return 1;
    }
    Scal* GetBasePtr() override {
      return &(*f)[IdxCell(0)];
    }
    FieldCell<Scal>* f;
  };
  // FieldCell<Vect>
  struct CoFcv : public Co {
    // f: vector field
    // i: component (0,1,2), or -1 for all
    CoFcv(FieldCell<Vect>* f, int d) : f(f), d(d) {}
    size_t GetSize() const override {
      return d == -1 ? Vect::dim : 1;
    }
    int GetIndex() const override {
      return d;
    }
    int GetStride() const override {
      return Vect::dim;
    }
    Scal* GetBasePtr() override {
      return &(*f)[IdxCell(0)][0];
    }
    FieldCell<Vect>* f;
    int d;
  };
  void Comm(const std::shared_ptr<Co>& r) {
    vcm_.push_back(r);
  }
  // Comm request for field f
  void Comm(FieldCell<Scal>* f) {
    Comm(std::make_shared<CoFcs>(f));
  }
  // Comm request for component d of f
  void Comm(FieldCell<Vect>* f, int d) {
    Comm(std::make_shared<CoFcv>(f, d));
  }
  // Comm request for all components of f
  void Comm(FieldCell<Vect>* f) {
    Comm(f, -1);
  }
  // f: scalar field
  // n: name
  void Dump(const FieldCell<Scal>* f, std::string n) {
    auto ff = const_cast<FieldCell<Scal>*>(f);
    vd_.push_back(std::make_pair(std::make_shared<CoFcs>(ff), n));
  }
  // f: vector field
  // d: component (0,1,2)
  // n: name
  void Dump(const FieldCell<Vect>* f, int d, std::string n) {
    auto ff = const_cast<FieldCell<Vect>*>(f);
    vd_.push_back(std::make_pair(std::make_shared<CoFcv>(ff, d), n));
  }
  void Dump(const std::shared_ptr<Co>& o, std::string name) {
    vd_.push_back(std::make_pair(o, name));
  }
  // Returns buffers for linear system
  void GetSolveTmp(
      std::vector<Scal>*& a, std::vector<Scal>*& b, std::vector<Scal>*& x) {
    a = &lsa_;
    b = &lsb_;
    x = &lsx_;
  }
  void Solve(const LS& ls) {
    vls_.push_back(ls);
  }
  const std::vector<std::shared_ptr<Co>>& GetComm() const {
    return vcm_;
  }
  const std::vector<std::pair<std::shared_ptr<Co>, std::string>>& GetDump()
      const {
    return vd_;
  }
  void ClearComm() {
    vcm_.clear();
  }
  void ClearDump() {
    vd_.clear();
  }

  using Rd = UReduce<Scal>;
  using Op = typename Rd::Op;
  template <class T>
  using OpT = typename Rd::template OpT<T>;
  using OpS = typename Rd::OpS;
  using OpSum = typename Rd::OpSum;
  using OpProd = typename Rd::OpProd;
  using OpMax = typename Rd::OpMax;
  using OpMin = typename Rd::OpMin;
  using OpSI = typename Rd::OpSI;
  using OpMinloc = typename Rd::OpMinloc;
  using OpMaxloc = typename Rd::OpMaxloc;
  template <class T>
  using OpVT = typename Rd::template OpT<T>;
  using OpCat = typename Rd::OpCat;
  template <class T>
  using OpCatT = typename Rd::template OpCatT<T>;
  template <class T>
  using OpCatVT = typename Rd::template OpCatVT<T>;

  void Reduce(const std::shared_ptr<Op>& o) {
    rd_.Add(o);
  }
  // u: src and dst buffer
  // o: operation
  void Reduce(Scal* u, std::string o) {
    rd_.Add(u, o);
  }
  const std::vector<std::shared_ptr<Op>>& GetReduce() const {
    return rd_.Get();
  }
  void ClearReduce() {
    rd_.Clear();
  }
  void Bcast(const std::shared_ptr<Op>& o) {
    bcast_.push_back(o);
  }
  const std::vector<std::shared_ptr<Op>>& GetBcast() const {
    return bcast_;
  }
  void ClearBcast() {
    bcast_.clear();
  }
  // Scatter request:
  // first: on root, vector for each block; ignored on others
  // second: receive buffer
  using ScatterRequest =
      std::pair<const std::vector<std::vector<Scal>>*, std::vector<Scal>*>;
  // scatter vo, send from root vo[block] to block vo[0]
  // vo: array of size comm_size on root
  //     array of size 1 on others
  void Scatter(const ScatterRequest& req) {
    scatter_.push_back(req);
  }
  const std::vector<ScatterRequest>& GetScatter() const {
    return scatter_;
  }
  void ClearScatter() {
    scatter_.clear();
  }
  const std::vector<LS>& GetSolve() const {
    return vls_;
  }
  Scal GetResidual() const {
    return lin_res_;
  }
  int GetIter() const {
    return lin_iter_;
  }
  void SetResidual(Scal res) {
    lin_res_ = res;
  }
  void SetIter(int iter) {
    lin_iter_ = iter;
  }
  void ClearSolve() {
    vls_.clear();
  }
  // Request timer report to file s
  void TimerReport(const std::string& s) {
    trep_ = s;
  }
  std::string GetTimerReport() const {
    return trep_;
  }
  void ClearTimerReport() {
    trep_ = "";
  }

 private:
  // b:Block, fc:FieldCell, ff:FieldFace, fn:FieldNode
  // inner
  const BlockCells bci_;
  const BlockFaces bfi_;
  const BlockNodes bni_;
  // all
  const BlockCells bca_;
  const BlockFaces bfa_;
  const BlockNodes bna_;
  // support
  const BlockCells bcs_;
  const BlockFaces bfs_;
  const BlockNodes bns_;
  // raw
  const IndexCells bcr_;
  const IndexFaces bfr_;
  const IndexNodes bnr_;

  const bool isroot_;
  const bool islead_;
  const MIdx mb_, me_; // begin,end of bci_
  const Rect<Vect> dom_; // domain covering bci_
  const Vect h_;
  const Vect hh_; // h_/2
  const Scal vol_;
  const MIdx gs_; // global mesh size
  const int id_; // unique id
  const Vect gl_; // global domain length
  bool checknan_; // CheckNan flag
  // pairs face,nci for which the halos cells
  // are set to nan after each communication
  std::vector<std::pair<IdxFace, size_t>> vfnan_;
  size_t maxcomm_ = 0; // maximum number of communication requests
  size_t edim_; // effective dimension
  std::array<Vect, dim> vs_; // surface vectors
  Vect va_; // surface area
  // cell neighbour cell
  std::array<size_t, kCellNumNeighbourFaces> cnc_;
  // cell neighbour face
  std::array<size_t, kCellNumNeighbourFaces> cnf_;
  // cell outward factor
  std::array<Scal, kCellNumNeighbourFaces> co_;
  // cell neighbour node
  std::array<size_t, kCellNumNeighbourNodes> cnn_;
  // face neighbour cell
  std::array<size_t, kFaceNumNeighbourCells * dim> fnc_;
  // face neighbour node
  std::array<size_t, kFaceNumNeighbourNodes * dim> fnn_;

  Suspender susp_;
  std::vector<std::shared_ptr<Co>> vcm_; // comm
  std::vector<std::pair<std::shared_ptr<Co>, std::string>> vd_; // dump
  std::string trep_; // timer report filename
  Rd rd_;
  std::vector<LS> vls_; // solve
  std::vector<std::shared_ptr<Op>> bcast_; // list of broadcast requests
  std::vector<ScatterRequest> scatter_; // list of scatter requests
  // tmp for GetSolveTmp()
  std::vector<Scal> lsa_; // matrix coeffs of size n * st.size()
  std::vector<Scal> lsb_; // rhs of size n
  std::vector<Scal> lsx_; // solution and initial guess of size n
  Scal lin_res_;
  int lin_iter_;

  Sampler samp_; // sample collector for histogram usage, always active
};

template <class _Scal, size_t _dim>
MeshStructured<_Scal, _dim>::MeshStructured(
    MIdx b, MIdx cs, Rect<Vect> dom, int hl, bool isroot, bool islead, MIdx gs,
    int id)
    // inner
    : bci_(b, cs)
    , bfi_(bci_.GetBegin(), bci_.GetSize())
    , bni_(bci_.GetBegin(), bci_.GetSize() + MIdx(1))
    // all
    , bca_(b - MIdx(hl), cs + MIdx(2 * hl))
    , bfa_(bca_.GetBegin(), bca_.GetSize())
    , bna_(bca_.GetBegin(), bca_.GetSize() + MIdx(1))
    // support
    , bcs_(bca_.GetBegin() + MIdx(1), bca_.GetSize() - MIdx(2))
    , bfs_(bcs_.GetBegin(), bcs_.GetSize())
    , bns_(bcs_.GetBegin(), bcs_.GetSize() + MIdx(1))
    // raw
    , bcr_(bca_.GetBegin(), bca_.GetSize() + MIdx(1))
    , bfr_(bcr_.GetBegin(), bcr_.GetSize())
    , bnr_(bcr_.GetBegin(), bcr_.GetSize())
    , isroot_(isroot)
    , islead_(islead)
    , mb_(bci_.GetBegin())
    , me_(bci_.GetEnd())
    , dom_(dom)
    , h_(dom.GetDimensions() / Vect(bci_.GetSize()))
    , hh_(h_ * 0.5)
    , vol_(h_.prod())
    , gs_(gs)
    , id_(id)
    , gl_(Vect(gs) * h_)
    , checknan_(false)
    , edim_(dim) {
  static_assert(dim == 3, "Not implemented for dim != 3");

  // mesh step

  // surface area
  va_[0] = h_[1] * h_[2];
  va_[1] = h_[2] * h_[0];
  va_[2] = h_[0] * h_[1];

  // surface vectors
  vs_[0] = Vect(va_[0], 0., 0.);
  vs_[1] = Vect(0., va_[1], 0.);
  vs_[2] = Vect(0., 0., va_[2]);

  // cell neighbour cell offset
  {
    MIdx w = bcr_.GetBegin(); // any cell
    IdxCell c = bcr_.GetIdx(w);
    for (auto q : Nci(c)) {
      MIdx wo;
      switch (q) {
        case 0:
          wo = MIdx(-1, 0, 0);
          break;
        case 1:
          wo = MIdx(1, 0, 0);
          break;
        case 2:
          wo = MIdx(0, -1, 0);
          break;
        case 3:
          wo = MIdx(0, 1, 0);
          break;
        case 4:
          wo = MIdx(0, 0, -1);
          break;
        default:
        case 5:
          wo = MIdx(0, 0, 1);
          break;
      };
      IdxCell cn = bcr_.GetIdx(w + wo);
      cnc_[q] = cn.GetRaw() - c.GetRaw();
    }
  }

  // cell neighbour face offset
  {
    MIdx w = bcr_.GetBegin(); // any cell
    IdxCell c = bcr_.GetIdx(w);
    for (auto q : Nci(c)) {
      MIdx wo;
      Dir d;
      switch (q) {
        case 0:
          wo = MIdx(0, 0, 0);
          d = Dir::i;
          break;
        case 1:
          wo = MIdx(1, 0, 0);
          d = Dir::i;
          break;
        case 2:
          wo = MIdx(0, 0, 0);
          d = Dir::j;
          break;
        case 3:
          wo = MIdx(0, 1, 0);
          d = Dir::j;
          break;
        case 4:
          wo = MIdx(0, 0, 0);
          d = Dir::k;
          break;
        default:
        case 5:
          wo = MIdx(0, 0, 1);
          d = Dir::k;
          break;
      };
      IdxFace fn = bfr_.GetIdx(std::make_pair(w + wo, d));
      cnf_[q] = fn.GetRaw() - c.GetRaw();
    }
  }

  // cell outward factor
  {
    for (size_t q = 0; q < kCellNumNeighbourFaces; ++q) {
      co_[q] = (q % 2 == 0 ? -1. : 1.);
    }
  }

  // cell neighbour node offset
  {
    MIdx w = bcr_.GetBegin(); // any cell
    IdxCell c = bcr_.GetIdx(w);
    auto qm = kCellNumNeighbourNodes;
    for (size_t q = 0; q < qm; ++q) {
      MIdx wo;
      switch (q) {
        case 0:
          wo = MIdx(0, 0, 0);
          break;
        case 1:
          wo = MIdx(1, 0, 0);
          break;
        case 2:
          wo = MIdx(0, 1, 0);
          break;
        case 3:
          wo = MIdx(1, 1, 0);
          break;
        case 4:
          wo = MIdx(0, 0, 1);
          break;
        case 5:
          wo = MIdx(1, 0, 1);
          break;
        case 6:
          wo = MIdx(0, 1, 1);
          break;
        default:
        case 7:
          wo = MIdx(1, 1, 1);
          break;
      };
      IdxNode nn = bnr_.GetIdx(w + wo);
      cnn_[q] = nn.GetRaw() - c.GetRaw();
    }
  }

  // face neighbour cell offset
  {
    MIdx w = bcr_.GetBegin(); // any cell
    for (size_t d = 0; d < dim; ++d) {
      IdxFace f = bfr_.GetIdx(w, Dir(d));
      auto qm = kFaceNumNeighbourCells;
      for (size_t q = 0; q < qm; ++q) {
        MIdx wo(0);
        switch (q) {
          case 0:
            wo[d] = -1;
            break;
          default:
          case 1:
            wo[d] = 0;
            break;
        };
        IdxCell cn = bcr_.GetIdx(w + wo);
        fnc_[d * qm + q] = cn.GetRaw() - f.GetRaw();
      }
    }
  }

  // face neighbour cell offset
  {
    MIdx w = bcr_.GetBegin(); // any cell
    for (size_t d = 0; d < dim; ++d) {
      IdxFace f = bfr_.GetIdx(w, Dir(d));
      auto qm = kFaceNumNeighbourNodes;
      for (size_t q = 0; q < qm; ++q) {
        MIdx wo;
        if (d == 0) {
          switch (q) {
            case 0:
              wo = MIdx(0, 0, 0);
              break;
            case 1:
              wo = MIdx(0, 1, 0);
              break;
            case 2:
              wo = MIdx(0, 1, 1);
              break;
            default:
            case 3:
              wo = MIdx(0, 0, 1);
              break;
          }
        } else if (d == 1) {
          switch (q) {
            case 0:
              wo = MIdx(0, 0, 0);
              break;
            case 1:
              wo = MIdx(0, 0, 1);
              break;
            case 2:
              wo = MIdx(1, 0, 1);
              break;
            default:
            case 3:
              wo = MIdx(1, 0, 0);
              break;
          }
        } else {
          switch (q) {
            case 0:
              wo = MIdx(0, 0, 0);
              break;
            case 1:
              wo = MIdx(1, 0, 0);
              break;
            case 2:
              wo = MIdx(1, 1, 0);
              break;
            default:
            case 3:
              wo = MIdx(0, 1, 0);
              break;
          }
        }
        IdxNode nn = bnr_.GetIdx(w + wo);
        fnn_[d * qm + q] = nn.GetRaw() - f.GetRaw();
      }
    }
  }
}

// Create uniform mesh
// domain: rectangle covering inner cells
// begin: index of lower inner cells
// s: number of inner cells in each direction
// hl: number of halo layers
// isroot: root block
// islead: lead block
// gs: global mesh size
// id: unique id
template <class M>
M InitUniformMesh(
    Rect<typename M::Vect> domain, typename M::MIdx begin, typename M::MIdx s,
    int hl, bool isroot, bool islead, typename M::MIdx gs, int id) {
  return M(begin, s, domain, hl, isroot, islead, gs, id);
}

template <class M>
M InitUniformMesh(
    const Rect<typename M::Vect>& domain, typename M::MIdx s,
    typename M::MIdx gs, int id) {
  using MIdx = typename M::MIdx;
  return InitUniformMesh<M>(domain, MIdx(0), s, true, true, gs, id);
}
