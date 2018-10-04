#pragma once

#include <vector>
#include <array>
#include <cassert>
#include <utility>
#include <limits>
#include <memory>
#include <stdexcept>

#include "util/suspender.h"
#include "idx.h"
#include "vect.h"
#include "range.h"
#include "field.h"
#include "map.h"
#include "blockface.h"
#include "rangein.h"
#include "distr/reduce.h"


// TODO: Neighbour faces iterator introducing (cell, face) pairs
// TODO: consider computing some on-the-fly to reduce memory access
template <class Scal_, size_t dim_>
class MeshStructured {
 public:
  using Scal = Scal_;
  static constexpr size_t dim = dim_;
  using Vect = GVect<Scal, dim>;
  using Dir = GDir<dim>;
  using MIdx = GMIdx<dim>;
  using BlockCells = GBlockCells<dim>;
  using BlockFaces = GBlockFaces<dim>;
  using BlockNodes = GBlockNodes<dim>;
  using IndexCells = GIndex<IdxCell,dim>;
  using IndexFaces = GIndex<IdxFace,dim>;
  using IndexNodes = GIndex<IdxNode,dim>;
  static constexpr size_t kCellNumNeighbourFaces = 2 * dim;
  static constexpr size_t kCellNumNeighbourNodes = std::pow(2, dim);
  static constexpr size_t kFaceNumNeighbourNodes = std::pow(2, dim - 1);
  static constexpr size_t kFaceNumNeighbourCells = 2;

  // b: begin, lower corner cell index
  // cs: inner cells size
  // dom: domain, rectangle covering inner cells
  // hl: halo cells from each side
  // isroot: root block
  // gs: global mesh size
  MeshStructured(MIdx b, MIdx cs, Rect<Vect> dom, int hl, bool isroot, MIdx gs);
  MIdx GetGlobalSize() const {
    return gs_;
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
  IdxCell GetNeighbourCell(IdxCell c, size_t q) const {
    assert(q < kCellNumNeighbourFaces);
    return IdxCell(size_t(c) + cnc_[q]);
  }
  IdxFace GetNeighbourFace(IdxCell c, size_t q) const {
    assert(q < kCellNumNeighbourFaces);
    return IdxFace(size_t(c) + cnf_[q]);
  }
  Scal GetOutwardFactor(IdxCell, size_t q) const {
    assert(q < kCellNumNeighbourFaces);
    return co_[q];
  }
  Vect GetOutwardSurface(IdxCell c, size_t n) const {
    assert(n < kCellNumNeighbourFaces);
    return GetSurface(GetNeighbourFace(c, n)) * GetOutwardFactor(c, n);
  }
  IdxNode GetNeighbourNode(IdxCell c, size_t q) const {
    assert(q < kCellNumNeighbourNodes);
    return IdxNode(size_t(c) + cnn_[q]);
  }
  Dir GetDir(IdxFace f) const {
    return bfr_.GetDir(f);
  }
  IdxCell GetNeighbourCell(IdxFace f, size_t q) const {
    auto qm = kFaceNumNeighbourCells;
    assert(q < qm);
    size_t d(bfr_.GetDir(f));
    return IdxCell(size_t(f) + fnc_[d * qm + q]);
  }
  Vect GetVectToCell(IdxFace f, size_t n) const {
    assert(n < kFaceNumNeighbourCells);
    return GetCenter(GetNeighbourCell(f, n)) - GetCenter(f);
  }
  IdxNode GetNeighbourNode(IdxFace f, size_t q) const {
    auto qm = kFaceNumNeighbourNodes;
    assert(q < qm);
    size_t d(bfr_.GetDir(f));
    return IdxNode(size_t(f) + fnn_[d * qm + q]);
  }
  size_t GetNumNeighbourFaces(IdxCell) const {
    return kCellNumNeighbourFaces;
  }
  // Neighbour cell indices
  GRange<size_t> Nci(IdxCell c) const {
    return GRange<size_t>(0, GetNumNeighbourFaces(c));
  }
  // Returns id of cell adjacent to c by face f.
  // -1 if f and c are not neighbours
  size_t GetNci(IdxCell c, IdxFace f) const {
    for (size_t q : Nci(c)) {
      if (GetNeighbourFace(c, q) == f) {
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
  bool IsInner(IdxCell c) const {
    MIdx w = bcr_.GetMIdx(c);
    return bci_.GetBegin() <= w && w < bci_.GetEnd();
  }

  // Returns range of raw indices
  template <class Idx>
  GRange<Idx> GetRaw() const { 
    return GetRaw((Idx*)0);
  }
  GRange<IdxCell> GetRaw(IdxCell*) const { 
    return GRange<IdxCell>(0, bcr_.size());
  }
  GRange<IdxFace> GetRaw(IdxFace*) const { 
    return GRange<IdxFace>(0, bfr_.size());
  }
  GRange<IdxNode> GetRaw(IdxNode*) const { 
    return GRange<IdxNode>(0, bnr_.size());
  }
  GRange<IdxCell> RawCells() const {
    return GetRaw<IdxCell>();
  }
  GRange<IdxFace> RawFaces() const {
    return GetRaw<IdxFace>();
  }
  GRange<IdxNode> RawNodes() const {
    return GetRaw<IdxNode>();
  }
  // Type-cast to GRange required for GField initialization
  template <class Idx>
  operator GRange<Idx>() const {
    return GetRaw<Idx>(); 
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

  bool IsInside(IdxCell c, Vect vect) const {
    for (auto q : Nci()) {
      IdxFace f = GetNeighbourFace(c, q);
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
  bool IsRoot() const { return isroot_; }

  // CheckNan flag
  bool CN() const { return checknan_; }
  void SetCN(bool c) { checknan_ = c; }

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

  // Comm request
  struct Co {
    virtual ~Co() {}
    // Number of scalar cell fields  (used in Dump)
    virtual size_t GetSize() const = 0;
  };
  // FieldCell<Scal>
  struct CoFcs : public Co {
    CoFcs(FieldCell<Scal>* f) : f(f) {}
    size_t GetSize() const { return 1; }
    FieldCell<Scal>* f;
  };
  // FieldCell<Vect>
  struct CoFcv : public Co {
    // f: vector field
    // i: component (0,1,2), or -1 for all
    CoFcv(FieldCell<Vect>* f, int d) : f(f), d(d) {}
    size_t GetSize() const { return d == -1 ? Vect::dim : 1; }
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
  void GetSolveTmp(std::vector<Scal>*& a, std::vector<Scal>*& b, 
                   std::vector<Scal>*& x) {
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
  const std::vector<std::pair<
      std::shared_ptr<Co>, std::string>>& GetDump() const {
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
  const std::vector<LS>& GetSolve() const {
    return vls_;
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
  const MIdx gs_;
  const MIdx mb_, me_; // begin,end of bci_
  const Rect<Vect> dom_; // domain covering bci_
  const Vect h_;
  const Vect hh_; // h_/2
  const Scal vol_;
  bool checknan_; // CheckNan flag
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
  std::vector<std::pair<std::shared_ptr<Co>, std::string>> vd_;  // dump
  std::string trep_; // timer report filename
  Rd rd_;
  std::vector<LS> vls_; // solve
  std::vector<std::shared_ptr<Op>> bcast_; // list of broadcast requests
  // tmp for GetSolveTmp()
  std::vector<Scal> lsa_; // matrix coeffs of size n * st.size()
  std::vector<Scal> lsb_; // rhs of size n
  std::vector<Scal> lsx_; // solution and initial guess of size n
};


template <class _Scal, size_t _dim>
MeshStructured<_Scal, _dim>::MeshStructured(
    MIdx b, MIdx cs, Rect<Vect> dom, int hl, bool isroot, MIdx gs)
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
    , gs_(gs)
    , mb_(bci_.GetBegin())
    , me_(bci_.GetEnd())
    , dom_(dom)
    , h_(dom.GetDimensions() / Vect(bci_.GetSize()))
    , hh_(h_ * 0.5)
    , vol_(h_.prod())
    , checknan_(false)
{
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
        case 0: wo = MIdx(-1, 0, 0); break;
        case 1: wo = MIdx(1, 0, 0); break;
        case 2: wo = MIdx(0, -1, 0); break;
        case 3: wo = MIdx(0, 1, 0); break;
        case 4: wo = MIdx(0, 0, -1); break;
        default: 
        case 5: wo = MIdx(0, 0, 1); break;
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
        case 0: wo = MIdx(0, 0, 0); d = Dir::i; break;
        case 1: wo = MIdx(1, 0, 0); d = Dir::i; break;
        case 2: wo = MIdx(0, 0, 0); d = Dir::j; break;
        case 3: wo = MIdx(0, 1, 0); d = Dir::j; break;
        case 4: wo = MIdx(0, 0, 0); d = Dir::k; break;
        default: 
        case 5: wo = MIdx(0, 0, 1); d = Dir::k; break;
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
        case 0: wo = MIdx(0, 0, 0); break;
        case 1: wo = MIdx(1, 0, 0); break;
        case 2: wo = MIdx(0, 1, 0); break;
        case 3: wo = MIdx(1, 1, 0); break;
        case 4: wo = MIdx(0, 0, 1); break;
        case 5: wo = MIdx(1, 0, 1); break;
        case 6: wo = MIdx(0, 1, 1); break;
        default: 
        case 7: wo = MIdx(1, 1, 1); break;
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
          case 0: wo[d] = -1;  break;
          default: 
          case 1: wo[d] = 0;  break;
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
            case 0: wo = MIdx(0,0,0); break;
            case 1: wo = MIdx(0,1,0); break;
            case 2: wo = MIdx(0,1,1); break;
            default:
            case 3: wo = MIdx(0,0,1); break;
          }
        } else if (d == 1) {
          switch (q) {
            case 0: wo = MIdx(0,0,0); break;
            case 1: wo = MIdx(0,0,1); break;
            case 2: wo = MIdx(1,0,1); break;
            default:
            case 3: wo = MIdx(1,0,0); break;
          }
        } else {
          switch (q) {
            case 0: wo = MIdx(0,0,0); break;
            case 1: wo = MIdx(1,0,0); break;
            case 2: wo = MIdx(1,1,0); break;
            default:
            case 3: wo = MIdx(0,1,0); break;
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
// gs: global mesh size
template <class M>
M InitUniformMesh(Rect<typename M::Vect> domain,
                     typename M::MIdx begin,
                     typename M::MIdx s, int hl, bool isroot,
                     typename M::MIdx gs) {
  return M(begin, s, domain, hl, isroot, gs);
}

template <class M>
M InitUniformMesh(const Rect<typename M::Vect>& domain,
                  typename M::MIdx s, typename M::MIdx gs) {
  using MIdx = typename M::MIdx;
  return InitUniformMesh<M>(domain, MIdx(0), s, true, gs);
}



