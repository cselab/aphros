// Created by Petr Karnakov on 30.05.2018
// Copyright 2018 ETH Zurich

#undef NDEBUG
#include <mpi.h>
#include <array>
#include <cassert>
#include <iostream>
#include <limits>
#include <sstream>

#include "distr/reduce.h"
#include "geom/blockface.h"
#include "geom/field.h"
#include "geom/idx.h"
#include "geom/map.h"
#include "geom/notation.h"
#include "geom/range.h"
#include "geom/rangein.h"
#include "geom/rangemulti.h"
#include "geom/transform.h"
#include "geom/vect.h"
#include "util/suspender.h"

namespace generic {

template <class T>
class Range {
 public:
  using Value = T;

  class iterator {
   public:
    explicit iterator(Value pos) : pos_(pos) {}
    iterator() = default;
    iterator(const iterator&) = default;
    iterator(iterator&&) = default;
    iterator& operator=(const iterator&) = default;
    iterator& operator=(iterator&&) = default;
    iterator& operator++() {
      ++pos_;
      return *this;
    }
    bool operator==(const iterator& other) const {
      return pos_ == other.pos_;
    }
    bool operator!=(const iterator& other) const {
      return pos_ != other.pos_;
    }
    Value operator*() const {
      return pos_;
    }

   private:
    Value pos_;
  };

  constexpr Range() : begin_(0), end_(0) {}
  constexpr explicit Range(Value end) : begin_(0), end_(end) {}
  constexpr Range(Value begin, Value end) : begin_(begin), end_(end) {}

  iterator begin() const {
    return iterator(begin_);
  }
  iterator end() const {
    return iterator(end_);
  }
  constexpr size_t size() const {
    return static_cast<size_t>(end_ - begin_);
  }
  void clear() {
    (*this) = Range();
  }
  constexpr bool operator==(const Range& other) const {
    return begin_ == other.begin_ && end_ == other.end_;
  }
  constexpr bool operator!=(const Range& other) const {
    return !(*this == other);
  }

 private:
  Value begin_;
  Value end_;
};

template <class Idx, size_t dim>
using Block = GBlock<Idx, dim>;
} // namespace generic

// Returns true if a < b (lex starting from end)
template <class T, size_t d>
bool Less(const generic::Vect<T, d>& a, const generic::Vect<T, d>& b) {
  int i = a.size();
  while (i--) {
    if (a[i] != b[i]) {
      return a[i] < b[i];
    }
  }
  return false;
}

namespace generic {

template <size_t n>
constexpr std::array<size_t, n> GetDirs();

template <>
constexpr std::array<size_t, 1> GetDirs<1>() {
  return {0};
}

template <>
constexpr std::array<size_t, 2> GetDirs<2>() {
  return {0, 1};
}

template <>
constexpr std::array<size_t, 3> GetDirs<3>() {
  return {0, 1, 2};
}

template <class Scal_, size_t dim_>
class MeshCartesian {
 public:
  using Scal = Scal_;
  static constexpr size_t dim = dim_;
  static constexpr generic::Range<size_t> dirs{dim};
  using Vect = generic::Vect<Scal, dim>;
  using Dir = GDir<dim>;
  using Direction = generic::Direction<Scal, dim>;
  using MIdx = GMIdx<dim>;
  using BlockCells = GBlockCells<dim>;
  using BlockFaces = GBlockFaces<dim>;
  using BlockNodes = GBlockNodes<dim>;
  using IndexCells = GIndex<IdxCell, dim>;
  using IndexFaces = GIndex<IdxFace, dim>;
  using IndexNodes = GIndex<IdxNode, dim>;
  using M = MeshCartesian;
  enum class Type { regular, cut, excluded };

  static constexpr size_t kCellNumNeighborFaces = dim * 2;
  static constexpr size_t kCellNumNeighborNodes = std::pow<size_t>(2, dim);
  static constexpr size_t kFaceNumNeighborNodes = std::pow<size_t>(2, dim - 1);
  static constexpr size_t kFaceNumNeighborCells = 2;
  static constexpr size_t kNumStencil = std::pow<size_t>(3, dim);
  static constexpr size_t kNumStencil5 = std::pow<size_t>(5, dim);
  static constexpr bool kIsEmbed = false;

  struct Flags {
    size_t edim = dim; // effective dimension
    bool linreport = false;
    bool check_nan = false; // check_nan field read by CHECKNAN macro
    bool check_symmetry = false;
    Scal check_symmetry_dump_threshold = 1e-5;
    Scal nan_faces_value = 1e100;
    std::array<bool, dim> is_periodic = {};
    Vect global_origin = Vect(0); // origin of global mesh
    MIdx global_blocks = MIdx(0); // number of blocks in global mesh
    Vect block_length = Vect(0); // length of one block (local mesh)
    static int GetIdFromBlock(MIdx block, MIdx global_blocks) {
      const auto& w = block;
      const auto& wmax = global_blocks;
      const int id = w[0] + w[1] * wmax[0] + w[2] * wmax[0] * wmax[1];
      return id;
    };
    int GetIdFromBlock(MIdx block) const {
      return GetIdFromBlock(block, global_blocks);
    }
    MIdx GetBlockFromPoint(Vect x) const {
      MIdx w((x - global_origin) / block_length);
      w = w.max(MIdx(0));
      w = w.min(global_blocks - MIdx(1));
      return w;
    }
    MPI_Comm comm;
  };

  // b: begin, lower corner cell index
  // cs: inner cells size
  // domain: domain, rectangle covering inner cells
  // halos: halo cells from each side
  // isroot: root block
  // gs: global mesh size
  // id: unique id
  MeshCartesian(
      MIdx b, MIdx cs, Rect<Vect> domain, int halos, bool isroot, bool islead,
      MIdx gs, int id);
  MeshCartesian(const MeshCartesian&) = delete;
  MeshCartesian(MeshCartesian&&);
  MeshCartesian& operator=(const MeshCartesian&) = delete;
  MeshCartesian& operator=(MeshCartesian&&) = delete;
  ~MeshCartesian();
  MIdx GetGlobalSize() const {
    return global_size_;
  }
  // Returns linear index of block between 0 and number of blocks - 1.
  int GetId() const {
    return id_;
  }
  int GetIdFromPoint(Vect x) const {
    return flags.GetIdFromBlock(flags.GetBlockFromPoint(x));
  }
  Vect GetGlobalLength() const {
    return global_length_;
  }
  Vect GetCellSize() const {
    return cell_size_;
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
    return indexc_;
  }
  const IndexFaces& GetIndexFaces() const {
    return indexf_;
  }
  const IndexNodes& GetIndexNodes() const {
    return indexn_;
  }

  // In blocks
  const BlockCells& GetInBlockCells() const {
    return blockci_;
  }
  const BlockFaces& GetInBlockFaces() const {
    return blockfi_;
  }
  const BlockNodes& GetInBlockNodes() const {
    return blockni_;
  }

  // Su blocks
  const BlockCells& GetSuBlockCells() const {
    return blockcs_;
  }
  const BlockFaces& GetSuBlockFaces() const {
    return blockfs_;
  }
  const BlockNodes& GetSuBlockNodes() const {
    return blockns_;
  }

  // All blocks
  const BlockCells& GetAllBlockCells() const {
    return blockca_;
  }
  const BlockFaces& GetAllBlockFaces() const {
    return blockfa_;
  }
  const BlockNodes& GetAllBlockNodes() const {
    return blockna_;
  }
  Vect GetCenter(IdxCell c) const {
    return fc_center_[c];
  }
  Rect<Vect> GetBoundingBox() const {
    return domain_;
  }
  Rect<Vect> GetGlobalBoundingBox() const {
    return Rect<Vect>(Vect(0), GetGlobalLength());
  }
  Vect GetCenter(IdxFace f) const {
    auto p = indexf_.GetMIdxDir(f);
    const MIdx& w = p.first;
    size_t d(p.second);
    Vect r =
        (domain_.low + half_cell_size_) + Vect(w - incells_begin_) * cell_size_;
    r[d] -= half_cell_size_[d];
    return r;
  }
  const Vect& GetSurface(IdxFace f) const {
    return face_surface_[size_t(indexf_.GetDir(f))];
  }
  Vect GetNode(IdxNode n) const {
    return domain_.low + Vect(indexn_.GetMIdx(n) - incells_begin_) * cell_size_;
  }
  Scal GetVolume(IdxCell) const {
    return cell_volume_;
  }
  Scal GetArea(IdxFace f) const {
    return face_area_[size_t(indexf_.GetDir(f))];
  }
  Scal GetAreaFraction(IdxFace) const {
    return 1;
  }
  Scal GetVolumeFraction(IdxCell) const {
    return 1;
  }
  IdxCell GetCell(IdxCell c, size_t q) const {
    assert(q < kCellNumNeighborFaces);
    return IdxCell(size_t(c) + cell_cell_[q]);
  }
  IdxFace GetFace(IdxCell c, size_t q) const {
    assert(q < kCellNumNeighborFaces);
    return IdxFace(size_t(c) + cell_face_[q]);
  }
  IdxFace GetFace(IdxFace f, size_t q) const {
    const auto qm = kCellNumNeighborFaces;
    assert(q < qm);
    return IdxFace(size_t(f) + face_face_[size_t(indexf_.GetDir(f)) * qm + q]);
  }
  Scal GetOutwardFactor(IdxCell, size_t q) const {
    assert(q < kCellNumNeighborFaces);
    return cell_outward_[q];
  }
  Vect GetOutwardSurface(IdxCell c, size_t n) const {
    assert(n < kCellNumNeighborFaces);
    return GetSurface(GetFace(c, n)) * GetOutwardFactor(c, n);
  }
  IdxNode GetNode(IdxCell c, size_t q) const {
    assert(q < kCellNumNeighborNodes);
    return IdxNode(size_t(c) + cell_node_[q]);
  }
  Dir GetDir(IdxFace f) const {
    return indexf_.GetDir(f);
  }
  IdxCell GetCell(IdxFace f, size_t q) const {
    auto qm = kFaceNumNeighborCells;
    assert(q < qm);
    return IdxCell(size_t(f) + face_cell_[size_t(indexf_.GetDir(f)) * qm + q]);
  }
  Vect GetVectToCell(IdxFace f, size_t n) const {
    assert(n < kFaceNumNeighborCells);
    return GetCenter(GetCell(f, n)) - GetCenter(f);
  }
  IdxNode GetNode(IdxFace f, size_t q) const {
    auto qm = kFaceNumNeighborNodes;
    assert(q < qm);
    size_t d(indexf_.GetDir(f));
    return IdxNode(size_t(f) + face_node_[d * qm + q]);
  }
  size_t GetNumFaces(IdxCell) const {
    return kCellNumNeighborFaces;
  }
  // Neighbor cell indices
  GRange<size_t> Nci(IdxCell c) const {
    return GRange<size_t>(0, GetNumFaces(c));
  }
  IdxCell GetCellFromPoint(Vect x) const {
    auto& ic = GetIndexCells();
    auto& bc = GetInBlockCells();
    const Vect h = GetCellSize();
    MIdx w((x - domain_.low) / h);
    w = w.max(MIdx(0));
    w = w.min(bc.GetSize() - MIdx(1));
    return ic.GetIdx(bc.GetBegin() + w);
  }
  // Cell indices over 3x3x3 stencil
  auto Stencil(IdxCell c) const {
    return MakeTransformIterator<IdxCell>(
        GRange<size_t>(0, kNumStencil),
        [this, c](size_t i) { return IdxCell(size_t(c) + stencil_[i]); });
  }
  const auto& GetStencilOffsets() const {
    return stencil_;
  }
  // Cell indices over 5x5x5 stencil
  auto Stencil5(IdxCell c) const {
    return MakeTransformIterator<IdxCell>(
        GRange<size_t>(0, kNumStencil5),
        [this, c](size_t i) { return IdxCell(size_t(c) + stencil5_[i]); });
  }
  const auto& GetStencil5Offsets() const {
    return stencil5_;
  }
  // Returns id of cell adjacent to c by face f.
  // -1 if f and c are not neighbors
  size_t GetNci(IdxCell c, IdxFace f) const {
    for (size_t q : Nci(c)) {
      if (GetFace(c, q) == f) {
        return q;
      }
    }
    return -1;
  }
  // Returns id of face opposite to q.
  size_t GetOpposite(size_t q) const {
    if (q == size_t(-1)) {
      return -1;
    }
    return q ^ 1;
  }
  size_t GetNumNodes(IdxCell) const {
    return kCellNumNeighborNodes;
  }
  size_t GetNumCells(IdxFace) const {
    return kFaceNumNeighborCells;
  }
  size_t GetNumNodes(IdxFace) const {
    return kFaceNumNeighborNodes;
  }
  // Returns column of cells cmm,cm,cp,cpp.
  // nci: 0 or 1 such that m.GetCell(f, nci) == cp
  std::array<IdxCell, 4> GetCellColumn(IdxFace f, size_t nci) const {
    const size_t d{GetIndexFaces().GetDir(f)};
    const IdxCell cm = GetCell(f, 1 - nci);
    const IdxCell cp = GetCell(f, nci);
    const IdxCell cmm = GetCell(cm, d * 2 + 1 - nci);
    const IdxCell cpp = GetCell(cp, d * 2 + nci);
    return {cmm, cm, cp, cpp};
  }
  // Returns true if cell is not a halo cell.
  bool IsInner(IdxCell c) const {
    MIdx w = indexc_.GetMIdx(c);
    return blockci_.GetBegin() <= w && w < blockci_.GetEnd();
  }
  bool IsBoundary(IdxFace f, size_t& nci) const {
    auto p = GetIndexFaces().GetMIdxDir(f);
    const size_t d(p.second);
    if (flags.is_periodic[d] || d > GetEdim()) {
      return false;
    }
    const auto w = p.first;
    if (w[d] == 0) {
      nci = 1;
      return true;
    } else if (w[d] == GetGlobalSize()[d]) {
      nci = 0;
      return true;
    }
    return false;
  }
  GIndex<IdxCell, dim> GetIndexer(IdxCell*) const {
    return indexc_;
  }
  GIndex<IdxFace, dim> GetIndexer(IdxFace*) const {
    return indexf_;
  }
  GIndex<IdxNode, dim> GetIndexer(IdxNode*) const {
    return indexn_;
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
  auto GetRangeIn() const {
    return GetRangeIn((Idx*)0);
  }
  auto GetRangeIn(IdxCell*) const {
    return GRangeIn<IdxCell, dim>(GetIndexCells(), GetInBlockCells());
  }
  auto GetRangeIn(IdxFace*) const {
    return GRangeIn<IdxFace, dim>(GetIndexFaces(), GetInBlockFaces());
  }
  auto GetRangeIn(IdxNode*) const {
    return GRangeIn<IdxNode, dim>(GetIndexNodes(), GetInBlockNodes());
  }
  auto Cells() const {
    return GetRangeIn<IdxCell>();
  }
  auto Cells4() const {
    return generic::RangeMulti<IdxCell, dim, 4>(
        indexc_.GetIdx(blockci_.GetBegin()), blockci_.GetSize(),
        indexc_.GetSize());
  }
  auto Faces() const {
    return GetRangeIn<IdxFace>();
  }
  auto Nodes() const {
    return GetRangeIn<IdxNode>();
  }

  // Returns range of support indices
  template <class Idx>
  auto GetRangeSu() const {
    return GetRangeSu((Idx*)0);
  }
  auto GetRangeSu(IdxCell*) const {
    return GRangeIn<IdxCell, dim>(GetIndexCells(), GetSuBlockCells());
  }
  auto GetRangeSu(IdxFace*) const {
    return GRangeIn<IdxFace, dim>(GetIndexFaces(), GetSuBlockFaces());
  }
  auto GetRangeSu(IdxNode*) const {
    return GRangeIn<IdxNode, dim>(GetIndexNodes(), GetSuBlockNodes());
  }
  auto SuCells() const {
    return GetRangeSu<IdxCell>();
  }
  auto SuCells4() const {
    return generic::RangeMulti<IdxCell, dim, 4>(
        indexc_.GetIdx(blockcs_.GetBegin()), blockcs_.GetSize(),
        indexc_.GetSize());
  }
  auto SuFaces() const {
    return GetRangeSu<IdxFace>();
  }
  auto SuNodes() const {
    return GetRangeSu<IdxNode>();
  }

  // Returns range of all indices
  template <class Idx>
  auto GetRangeAll() const {
    return GetRangeAll((Idx*)0);
  }
  auto GetRangeAll(IdxCell*) const {
    return GRangeIn<IdxCell, dim>(GetIndexCells(), GetAllBlockCells());
  }
  auto GetRangeAll(IdxFace*) const {
    return GRangeIn<IdxFace, dim>(GetIndexFaces(), GetAllBlockFaces());
  }
  auto GetRangeAll(IdxNode*) const {
    return GRangeIn<IdxNode, dim>(GetIndexNodes(), GetAllBlockNodes());
  }
  auto AllCells() const {
    return GetRangeAll<IdxCell>();
  }
  auto AllFaces() const {
    return GetRangeAll<IdxFace>();
  }
  auto AllNodes() const {
    return GetRangeAll<IdxNode>();
  }

  // Expression on face: v[0] * cm + v[1] * cp + v[2]
  using ExprFace = generic::Vect<Scal, 3>;
  // Expression on cell: v[0] * c + v[1] * cxm + ... + v[6] * czp + v[7]
  using Expr = generic::Vect<Scal, M::dim * 2 + 2>;

  // Wrappers to satisfy the interface of Embed<M>
  bool IsCell(IdxCell) const {
    return true;
  }
  bool IsCell(IdxFace) const {
    return false;
  }
  std::vector<Vect> GetFacePoly(IdxFace f) const {
    std::vector<Vect> xx;
    for (size_t e = 0; e < GetNumNodes(f); ++e) {
      auto n = GetNode(f, e);
      xx.push_back(GetNode(n));
    }
    return xx;
  }
  std::vector<Vect> GetFacePoly(IdxCell) const {
    return {};
  }
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
  Vect GetNormal(IdxCell) const {
    return GetNan<Vect>();
  }
  Type GetType(IdxCell) const {
    return Type::regular;
  }
  Type GetType(IdxFace) const {
    return Type::regular;
  }
  bool IsRegular(IdxCell) const {
    return true;
  }
  bool IsCut(IdxCell) const {
    return false;
  }
  bool IsExcluded(IdxCell) const {
    return false;
  }
  IdxCell GetRegularNeighbor(IdxCell c) const {
    return c;
  }
  Vect GetSurface(IdxCell) const {
    return GetNan<Vect>();
  }
  Scal GetArea(IdxCell) const {
    return GetNan<Scal>();
  }
  template <class F>
  void LoopFaces(F lambda) const {
    for (auto f : Faces()) {
      lambda(f);
    }
  }
  template <class F>
  void LoopSuFaces(F lambda) const {
    for (auto f : SuFaces()) {
      lambda(f);
    }
  }
  template <class F>
  void LoopNci(IdxCell c, F lambda) const {
    for (auto q : Nci(c)) {
      lambda(q);
    }
  }
  template <class F>
  void LoopNciFaces(IdxCell c, F lambda) const {
    for (auto q : Nci(c)) {
      lambda(q);
    }
  }
  template <class F>
  void LoopNciEmbed(IdxCell, F) const {
    return;
  }
  void AppendExpr(Expr& sum, const ExprFace& v, size_t q) const {
    sum[0] += v[1 - q % 2];
    sum[1 + q] += v[q % 2];
    sum.back() += v.back();
  }
  void AppendExpr(Expr& sum, const ExprFace& v, size_t q, IdxCell) const {
    AppendExpr(sum, v, q);
  }
  bool IsInside(IdxCell c, Vect vect) const {
    for (auto q : Nci(c)) {
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
  // Pairs face,nci for which the halos cells
  // are set to nan after each communication
  const std::vector<std::pair<IdxFace, size_t>>& GetNanFaces() const;
  void SetNanFaces(const std::vector<std::pair<IdxFace, size_t>>& vfnan);
  // Fills halo cell with garbage.
  // Using actual NaNs not allowed since some code relies on u*0 == 0
  template <class T>
  void ApplyNanFaces(FieldCell<T>& fc);

  size_t GetEdim() const {
    return flags.edim;
  }
  MeshCartesian& GetMesh() {
    return *this;
  }
  const MeshCartesian& GetMesh() const {
    return *this;
  }

  // Notation
  static Direction direction(size_t i) {
    return Direction(i);
  }
  static Direction direction(size_t i, size_t forward) {
    return Direction(i, forward);
  }
  IdxCellMesh<M> operator()(IdxCell c) const {
    return IdxCellMesh<M>(c, *this);
  }
  IdxFaceMesh<M> operator()(IdxFace f) const {
    return IdxFaceMesh<M>(f, *this);
  }
  IdxCellMesh<M> cell() const {
    return IdxCellMesh<M>(*this);
  }
  IdxFaceMesh<M> face() const {
    return IdxFaceMesh<M>(*this);
  }
  auto CellsM() const {
    return MakeTransformIterator<IdxCellMesh<M>>(
        Cells(), [this](IdxCell c) { return IdxCellMesh<M>(c, *this); });
  }
  auto SuCellsM() const {
    return MakeTransformIterator<IdxCellMesh<M>>(
        SuCells(), [this](IdxCell c) { return IdxCellMesh<M>(c, *this); });
  }
  auto FacesM() const {
    return MakeTransformIterator<IdxFaceMesh<M>>(
        Faces(), [this](IdxFace f) { return IdxFaceMesh<M>(f, *this); });
  }
  auto SuFacesM() const {
    return MakeTransformIterator<IdxFaceMesh<M>>(
        SuFaces(), [this](IdxFace f) { return IdxFaceMesh<M>(f, *this); });
  }

 public:
  using Sem = Suspender::Sem;
  Sem GetSem(std::string name = "") {
    return susp_.GetSem(name);
  }
  bool Pending() const {
    return susp_.Pending();
  }
  const Suspender& GetSuspender() const {
    return susp_;
  }

  struct CommRequest {
    virtual ~CommRequest() {}
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
  struct CommRequestScal : public CommRequest {
    CommRequestScal(FieldCell<Scal>* field_) : field(field_) {}
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
      return &(*field)[IdxCell(0)];
    }
    FieldCell<Scal>* field;
  };
  // FieldCell<Vect>
  struct CommRequestVect : public CommRequest {
    // f: vector field
    // i: component (0,1,2), or -1 for all
    CommRequestVect(FieldCell<Vect>* field_, int d_) : field(field_), d(d_) {}
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
      return &(*field)[IdxCell(0)][0];
    }
    FieldCell<Vect>* field;
    int d;
  };
  void Comm(std::unique_ptr<CommRequest>&& r);
  void Comm(FieldCell<Scal>* field);
  void Comm(FieldCell<Vect>* field, int component);
  void Comm(FieldCell<Vect>* field);
  const std::vector<std::unique_ptr<CommRequest>>& GetComm() const;
  void ClearComm();

  void Dump(const FieldCell<Scal>* field, std::string name);
  void Dump(const FieldCell<Vect>* field, int component, std::string name);
  void Dump(std::unique_ptr<CommRequest>&& request, std::string name);
  const std::vector<std::pair<std::unique_ptr<CommRequest>, std::string>>&
  GetDump() const;
  void ClearDump();

  using Op = typename UReduce<Scal>::Op;

  MPI_Comm GetMpiComm() const {
    return flags.comm;
  }
  void Reduce(std::unique_ptr<Op>&& o);
  void Reduce(Scal* u, std::string o);
  const std::vector<std::unique_ptr<Op>>& GetReduce() const;
  void Reduce(Scal* buf, ReductionType::Sum);
  void Reduce(Scal* buf, ReductionType::Prod);
  void Reduce(Scal* buf, ReductionType::Max);
  void Reduce(Scal* buf, ReductionType::Min);
  void Reduce(std::pair<Scal, int>* buf, ReductionType::MaxLoc);
  void Reduce(std::pair<Scal, int>* buf, ReductionType::MinLoc);
  template <class T>
  void Reduce(std::vector<T>* buf, ReductionType::Concat) {
    Reduce(std::make_unique<typename UReduce<Scal>::template OpCatT<T>>(buf));
  }
  template <class T>
  void Reduce(std::vector<std::vector<T>>* buf, ReductionType::Concat) {
    Reduce(std::make_unique<typename UReduce<Scal>::template OpCatVT<T>>(buf));
  }
  void ClearReduce();
  void ReduceToLead(std::unique_ptr<Op>&& o);
  template <class T>
  void GatherToLead(std::vector<T>* buf) {
    ReduceToLead(
        std::make_unique<typename UReduce<Scal>::template OpCatT<T>>(buf));
  }
  template <class T>
  void GatherToLead(std::vector<std::vector<T>>* buf) {
    ReduceToLead(
        std::make_unique<typename UReduce<Scal>::template OpCatVT<T>>(buf));
  }
  const std::vector<std::unique_ptr<Op>>& GetReduceToLead() const;
  void ClearReduceToLead();
  void Bcast(std::unique_ptr<Op>&& o);
  template <class T>
  void Bcast(std::vector<T>* buf) {
    Bcast(std::make_unique<typename UReduce<Scal>::template OpCatT<T>>(buf));
  }
  template <class T>
  void Bcast(std::vector<std::vector<T>>* buf) {
    Bcast(std::make_unique<typename UReduce<Scal>::template OpCatVT<T>>(buf));
  }
  const std::vector<std::unique_ptr<Op>>& GetBcast() const;
  void ClearBcast();
  // Scatter request:
  // first: sendbuffer on root, vector for each block; ignored on others
  // second: receive buffer
  using ScatterRequest =
      std::pair<const std::vector<std::vector<Scal>>*, std::vector<Scal>*>;
  void Scatter(const ScatterRequest& req);
  const std::vector<ScatterRequest>& GetScatter() const;
  void ClearScatter();
  // Request timer report to file s
  void TimerReport(const std::string& s) {
    timer_report_path_ = s;
  }
  std::string GetTimerReport() const {
    return timer_report_path_;
  }
  void ClearTimerReport() {
    timer_report_path_ = "";
  }

 public:
  Flags flags;

 private:
  // b:Block, fc:FieldCell, ff:FieldFace, fn:FieldNode
  // inner
  const BlockCells blockci_;
  const BlockFaces blockfi_;
  const BlockNodes blockni_;
  // all
  const BlockCells blockca_;
  const BlockFaces blockfa_;
  const BlockNodes blockna_;
  // support
  const BlockCells blockcs_;
  const BlockFaces blockfs_;
  const BlockNodes blockns_;
  // index
  const IndexCells indexc_;
  const IndexFaces indexf_;
  const IndexNodes indexn_;

  const bool isroot_;
  const bool islead_;
  const MIdx incells_begin_, incells_end_;
  const Rect<Vect> domain_; // domain covering blockci_
  const Vect cell_size_;
  const Vect half_cell_size_;
  const Scal cell_volume_;
  const MIdx global_size_;
  const int id_; // unique id
  const Vect global_length_; // global domain length
  bool checknan_; // CheckNan flag
  // pairs face,nci for which the halos cells
  // are set to nan after each communication
  std::array<Vect, dim> face_surface_; // surface vectors
  Vect face_area_;
  // offsets
  std::array<size_t, kCellNumNeighborFaces> cell_cell_;
  std::array<size_t, kCellNumNeighborFaces> cell_face_;
  std::array<Scal, kCellNumNeighborFaces> cell_outward_;
  std::array<size_t, kCellNumNeighborNodes> cell_node_;
  std::array<size_t, kFaceNumNeighborCells * dim> face_cell_;
  std::array<size_t, kCellNumNeighborFaces * dim> face_face_;
  std::array<size_t, kFaceNumNeighborNodes * dim> face_node_;
  std::array<size_t, kNumStencil> stencil_; // 3x3x3 stencil
  std::array<size_t, kNumStencil5> stencil5_; // 5x5x5 stencil

  FieldCell<Vect> fc_center_;

  Suspender susp_;
  std::string timer_report_path_;
  struct Imp;
  std::unique_ptr<Imp> imp;
};

template <class Scal, size_t dim>
constexpr generic::Range<size_t> MeshCartesian<Scal, dim>::dirs;

template <class _Scal, size_t _dim>
struct MeshCartesian<_Scal, _dim>::Imp {
  using Owner = MeshCartesian<_Scal, _dim>;
  Imp(Owner* owner_) : owner(owner_) {}

  // sw: stencil half-width, results in stencil [-sw,sw]
  template <size_t sw>
  TransformIterator<IdxCell, GBlock<size_t, dim>> StencilGeneral(
      IdxCell c) const {
    constexpr size_t sn = sw * 2 + 1;
    const GBlock<size_t, dim> bo(MIdx(-sw), MIdx(sn));
    return MakeTransformIterator<IdxCell>(bo, [this, c](MIdx wo) {
      auto& indexc = owner->GetIndexCells();
      return indexc.GetIdx(indexc.GetMIdx(c) + wo);
    });
  }

  Owner* owner;

  // Requests
  std::vector<std::unique_ptr<CommRequest>> commreq;
  std::vector<std::pair<std::unique_ptr<CommRequest>, std::string>> dump;
  UReduce<Scal> reduce;
  UReduce<Scal> reduce_lead;
  std::vector<std::unique_ptr<typename UReduce<Scal>::Op>> bcast;
  std::vector<ScatterRequest> scatter;
  std::vector<std::pair<IdxFace, size_t>> vfnan;
};

template <class _Scal, size_t _dim>
MeshCartesian<_Scal, _dim>::MeshCartesian(
    MIdx b, MIdx cs, Rect<Vect> domain, int halos, bool isroot, bool islead,
    MIdx gs, int id)
    // inner
    : blockci_(b, cs)
    , blockfi_(blockci_.GetBegin(), blockci_.GetSize())
    , blockni_(blockci_.GetBegin(), blockci_.GetSize() + MIdx(1))
    // all
    , blockca_(b - MIdx(halos), cs + MIdx(2 * halos))
    , blockfa_(blockca_.GetBegin(), blockca_.GetSize())
    , blockna_(blockca_.GetBegin(), blockca_.GetSize() + MIdx(1))
    // support
    , blockcs_(blockca_.GetBegin() + MIdx(1), blockca_.GetSize() - MIdx(2))
    , blockfs_(blockcs_.GetBegin(), blockcs_.GetSize())
    , blockns_(blockcs_.GetBegin(), blockcs_.GetSize() + MIdx(1))
    // index
    , indexc_(blockca_.GetBegin(), blockca_.GetSize() + MIdx(1))
    , indexf_(indexc_.GetBegin(), indexc_.GetSize())
    , indexn_(indexc_.GetBegin(), indexc_.GetSize())
    , isroot_(isroot)
    , islead_(islead)
    , incells_begin_(blockci_.GetBegin())
    , incells_end_(blockci_.GetEnd())
    , domain_(domain)
    , cell_size_(domain.GetDimensions() / Vect(blockci_.GetSize()))
    , half_cell_size_(cell_size_ * 0.5)
    , cell_volume_(cell_size_.prod())
    , global_size_(gs)
    , id_(id)
    , global_length_(Vect(gs) * cell_size_)
    , face_area_(Vect(cell_volume_) / cell_size_)
    , imp(new Imp(this)) {
  // surface vectors
  for (auto d : dirs) {
    face_surface_[d] = Vect::GetUnit(d) * face_area_[d];
  }

  {
    const MIdx w = indexc_.GetBegin();
    const IdxCell c = indexc_.GetIdx(w);
    for (auto d : dirs) {
      const auto wd = MIdx::GetUnit(d);

      // cell neighbor cell offset
      cell_cell_[2 * d] = size_t(indexc_.GetIdx(w - wd)) - size_t(c);
      cell_cell_[2 * d + 1] = size_t(indexc_.GetIdx(w + wd)) - size_t(c);

      // cell neighbor face offset
      cell_face_[2 * d] = size_t(indexf_.GetIdx(w, Dir(d))) - size_t(c);
      cell_face_[2 * d + 1] =
          size_t(indexf_.GetIdx(w + wd, Dir(d))) - size_t(c);

      // cell outward factor
      cell_outward_[2 * d] = -1;
      cell_outward_[2 * d + 1] = 1;
    }
  }

  { // cell neighbor node offset
    const MIdx w = indexc_.GetBegin();
    const IdxCell c = indexc_.GetIdx(w);
    const generic::Block<size_t, dim> block(MIdx(2));
    size_t i = 0;
    for (MIdx offset : block) {
      cell_node_[i++] = size_t(indexn_.GetIdx(w + offset)) - size_t(c);
    }
  }

  { // face neighbor cell offset
    const MIdx w = indexc_.GetBegin();
    for (auto d : dirs) {
      const auto wd = MIdx::GetUnit(d);
      const IdxFace f = indexf_.GetIdx(w, Dir(d));
      face_cell_[2 * d] = size_t(indexc_.GetIdx(w - wd)) - size_t(f);
      face_cell_[2 * d + 1] = size_t(indexc_.GetIdx(w)) - size_t(f);
    }
  }

  { // face neighbor face offset
    const MIdx w = indexc_.GetBegin();
    for (auto df : dirs) {
      const IdxFace f = indexf_.GetIdx(w, Dir(df));
      for (auto d : dirs) {
        const auto wd = MIdx::GetUnit(d);
        face_face_[df * kCellNumNeighborFaces + 2 * d] =
            size_t(indexf_.GetIdx(w - wd, Dir(df))) - size_t(f);

        face_face_[df * kCellNumNeighborFaces + 2 * d + 1] =
            size_t(indexf_.GetIdx(w + wd, Dir(df))) - size_t(f);
      }
    }
  }

  { // face neighbor node offset
    using HyperMIdx = GMIdx<dim - 1>;
    std::array<HyperMIdx, kFaceNumNeighborNodes> hyper_offsets{};
    // hyper_block:
    // 2D: (0) (1)
    // 3D: (0 0) (1 0) (0 1) (1 1)
    // 4D: (0 0 0) (1 0 0) (0 1 0) (1 1 0) (0 0 1) (1 0 1) (0 1 1) (1 1 1)
    // hyper_offsets:
    // 2D: (0) (1)
    // 3D: (0 0) (1 0) (1 1) (0 1)
    // 4D: (0 0 0) (1 0 0) (1 1 0) (0 1 0) (0 0 1) (1 0 1) (1 1 1) (0 1 1)
    {
      const generic::Block<size_t, dim - 1> hyper_block(HyperMIdx(2));
      fassert_equal(hyper_block.size(), kFaceNumNeighborNodes);
      size_t i = 0;
      for (auto wm : hyper_block) {
        // swap every second pair
        hyper_offsets[i ^ ((i & 2) >> 1)] = wm;
        ++i;
      }
    }

    const MIdx base = indexc_.GetBegin();
    for (auto d : dirs) {
      const IdxFace f = indexf_.GetIdx(base, Dir(d));
      for (size_t q = 0; q < kFaceNumNeighborNodes; ++q) {
        MIdx offset(0);
        for (size_t dm = 0; dm + 1 < dim; ++dm) {
          offset[(d + dm + 1) % dim] = hyper_offsets[q][dm];
        }
        face_node_[d * kFaceNumNeighborNodes + q] =
            size_t(indexn_.GetIdx(base + offset)) - size_t(f);
      }
    }
  }
  { // stencil 3x3x3 offsets
    const IdxCell c = indexc_.GetIdx(blockci_.GetBegin());
    size_t i = 0;
    for (auto cn : imp->template StencilGeneral<1>(c)) {
      stencil_[i++] = size_t(cn) - size_t(c);
    }
    fassert_equal(i, kNumStencil);
  }
  { // stencil 5x5x5 offsets
    const IdxCell c = indexc_.GetIdx(blockci_.GetBegin());
    size_t i = 0;
    for (auto cn : imp->template StencilGeneral<2>(c)) {
      stencil5_[i++] = size_t(cn) - size_t(c);
    }
    fassert_equal(i, kNumStencil5);
  }

  { // cell centers
    fc_center_.Reinit(*this);
    for (auto c : AllCells()) {
      fc_center_[c] = (domain_.low + half_cell_size_) +
                      Vect(indexc_.GetMIdx(c) - incells_begin_) * cell_size_;
    }
  }
}

template <class Scal, size_t dim>
MeshCartesian<Scal, dim>::~MeshCartesian() = default;

template <class Scal, size_t dim>
MeshCartesian<Scal, dim>::MeshCartesian(MeshCartesian&&) = default;

template <class M>
M InitMesh(
    Rect<typename M::Vect> domain, typename M::MIdx begin, typename M::MIdx s,
    int halos, bool isroot, bool islead, typename M::MIdx gs, int id) {
  return {begin, s, domain, halos, isroot, islead, gs, id};
}

template <class Scal, size_t dim>
void MeshCartesian<Scal, dim>::Comm(std::unique_ptr<CommRequest>&& r) {
  imp->commreq.emplace_back(std::move(r));
}
template <class Scal, size_t dim>
void MeshCartesian<Scal, dim>::Comm(FieldCell<Scal>* f) {
  Comm(std::make_unique<CommRequestScal>(f));
}
template <class Scal, size_t dim>
void MeshCartesian<Scal, dim>::Comm(FieldCell<Vect>* f, int d) {
  Comm(std::make_unique<CommRequestVect>(f, d));
}
template <class Scal, size_t dim>
void MeshCartesian<Scal, dim>::Comm(FieldCell<Vect>* f) {
  Comm(f, -1);
}
template <class Scal, size_t dim>
void MeshCartesian<Scal, dim>::Dump(const FieldCell<Scal>* f, std::string n) {
  auto ff = const_cast<FieldCell<Scal>*>(f);
  imp->dump.emplace_back(std::make_unique<CommRequestScal>(ff), n);
}
template <class Scal, size_t dim>
void MeshCartesian<Scal, dim>::Dump(
    const FieldCell<Vect>* f, int d, std::string n) {
  auto ff = const_cast<FieldCell<Vect>*>(f);
  imp->dump.emplace_back(std::make_unique<CommRequestVect>(ff, d), n);
}
template <class Scal, size_t dim>
void MeshCartesian<Scal, dim>::Dump(
    std::unique_ptr<CommRequest>&& o, std::string name) {
  imp->dump.emplace_back(std::move(o), name);
}
template <class Scal, size_t dim>
auto MeshCartesian<Scal, dim>::GetComm() const
    -> const std::vector<std::unique_ptr<CommRequest>>& {
  return imp->commreq;
}
template <class Scal, size_t dim>
auto MeshCartesian<Scal, dim>::GetDump() const -> const
    std::vector<std::pair<std::unique_ptr<CommRequest>, std::string>>& {
  return imp->dump;
}
template <class Scal, size_t dim>
void MeshCartesian<Scal, dim>::ClearComm() {
  imp->commreq.clear();
}
template <class Scal, size_t dim>
void MeshCartesian<Scal, dim>::ClearDump() {
  imp->dump.clear();
}
template <class Scal, size_t dim>
void MeshCartesian<Scal, dim>::Reduce(
    std::unique_ptr<typename UReduce<Scal>::Op>&& o) {
  imp->reduce.Add(std::move(o));
}
template <class Scal, size_t dim>
void MeshCartesian<Scal, dim>::Reduce(Scal* u, std::string o) {
  imp->reduce.Add(u, o);
}
template <class Scal, size_t dim>
const std::vector<std::unique_ptr<typename UReduce<Scal>::Op>>&
MeshCartesian<Scal, dim>::GetReduce() const {
  return imp->reduce.Get();
}
template <class Scal, size_t dim>
void MeshCartesian<Scal, dim>::Reduce(Scal* buf, ReductionType::Sum) {
  Reduce(buf, "sum");
}
template <class Scal, size_t dim>
void MeshCartesian<Scal, dim>::Reduce(Scal* buf, ReductionType::Prod) {
  Reduce(buf, "prod");
}
template <class Scal, size_t dim>
void MeshCartesian<Scal, dim>::Reduce(Scal* buf, ReductionType::Max) {
  Reduce(buf, "max");
}
template <class Scal, size_t dim>
void MeshCartesian<Scal, dim>::Reduce(Scal* buf, ReductionType::Min) {
  Reduce(buf, "min");
}
template <class Scal, size_t dim>
void MeshCartesian<Scal, dim>::Reduce(
    std::pair<Scal, int>* buf, ReductionType::MaxLoc) {
  Reduce(std::make_unique<typename UReduce<Scal>::OpMaxloc>(buf));
}
template <class Scal, size_t dim>
void MeshCartesian<Scal, dim>::Reduce(
    std::pair<Scal, int>* buf, ReductionType::MinLoc) {
  Reduce(std::make_unique<typename UReduce<Scal>::OpMinloc>(buf));
}
template <class Scal, size_t dim>
void MeshCartesian<Scal, dim>::ClearReduce() {
  imp->reduce.Clear();
}
template <class Scal, size_t dim>
void MeshCartesian<Scal, dim>::ReduceToLead(
    std::unique_ptr<typename UReduce<Scal>::Op>&& o) {
  imp->reduce_lead.Add(std::move(o));
}
template <class Scal, size_t dim>
const std::vector<std::unique_ptr<typename UReduce<Scal>::Op>>&
MeshCartesian<Scal, dim>::GetReduceToLead() const {
  return imp->reduce_lead.Get();
}
template <class Scal, size_t dim>
void MeshCartesian<Scal, dim>::ClearReduceToLead() {
  imp->reduce_lead.Clear();
}

template <class Scal, size_t dim>
void MeshCartesian<Scal, dim>::Scatter(const ScatterRequest& req) {
  imp->scatter.push_back(req);
}
template <class Scal, size_t dim>
auto MeshCartesian<Scal, dim>::GetScatter() const
    -> const std::vector<ScatterRequest>& {
  return imp->scatter;
}
template <class Scal, size_t dim>
void MeshCartesian<Scal, dim>::ClearScatter() {
  imp->scatter.clear();
}
template <class Scal, size_t dim>
void MeshCartesian<Scal, dim>::Bcast(
    std::unique_ptr<typename UReduce<Scal>::Op>&& o) {
  imp->bcast.emplace_back(std::move(o));
}
template <class Scal, size_t dim>
auto MeshCartesian<Scal, dim>::GetBcast() const
    -> const std::vector<std::unique_ptr<Op>>& {
  return imp->bcast;
}
template <class Scal, size_t dim>
void MeshCartesian<Scal, dim>::ClearBcast() {
  imp->bcast.clear();
}
template <class Scal, size_t dim>
auto MeshCartesian<Scal, dim>::GetNanFaces() const
    -> const std::vector<std::pair<IdxFace, size_t>>& {
  return imp->vfnan;
}
template <class Scal, size_t dim>
void MeshCartesian<Scal, dim>::SetNanFaces(
    const std::vector<std::pair<IdxFace, size_t>>& vfnan) {
  imp->vfnan = vfnan;
}
template <class Scal, size_t dim>
template <class T>
void MeshCartesian<Scal, dim>::ApplyNanFaces(FieldCell<T>& fc) {
  for (auto p : imp->vfnan) {
    IdxFace f = p.first;
    size_t nci = p.second;
    auto cc = GetCellColumn(f, nci);
    fc[cc[0]] = T(flags.nan_faces_value);
    fc[cc[1]] = T(flags.nan_faces_value);
  }
}

} // namespace generic

const int dim = 4;
using MIdx = GMIdx<dim>;
using Dir = GDir<dim>;
using Scal = double;
using Vect = generic::Vect<Scal, dim>;

void TestBlock() {
  std::cout << __func__ << std::endl;
  const size_t halos = 1;
  MIdx oi(0); // origin inner
  MIdx si(2); // size inner
  MIdx oa = oi - MIdx(halos); // origin all
  MIdx sa = si + MIdx(2 * halos); // size all

  GBlockFaces<dim> bi(oi, si);
  GIndex<IdxFace, dim> ba(oa, sa);

  GRangeIn<IdxFace, dim> ri(ba, bi);

  const MIdx xp0 = oa - MIdx(Dir(0));
  MIdx xp = xp0;
  Dir dp(0); // direction

  // Check that whole inner block covered with ascending indices
  for (auto i : ri) {
    auto x = ba.GetMIdx(i);
    auto d = ba.GetDir(i);

    // Next direction, reset xp
    if (dp < d) {
      xp = xp0;
    }

    // std::cerr << x << " " << d.GetLetter() << std::endl;
    assert(Less(xp, x));
    assert(oi <= x && x < oi + si + MIdx(d));

    xp = x;
    dp = d;
  }
}

bool Cmp(Scal a, Scal b) {
  return std::abs(a - b) < 1e-12;
}

template <class T>
bool Cmp(T* a, T* b) {
  return a == b;
}

bool Cmp(size_t a, size_t b) {
  return a == b;
}

#define CMP(a, b) assert(Cmp(a, b));

// Print CMP
#define PCMP(a, b)                                                    \
  std::cerr << #a << "=" << a << ", " << #b << "=" << b << std::endl; \
  CMP(a, b);

void TestMesh() {
  std::cout << __func__ << std::endl;
  const Rect<Vect> dom(Vect(0), Vect(1));
  using M = generic::MeshCartesian<Scal, dim>;
  const MIdx begin(0);
  const MIdx size(2);
  const int halos = 2;
  Vect doms = dom.GetDimensions();
  const Vect h = dom.GetDimensions() / Vect(size);
  M m = generic::InitMesh<M>(dom, begin, size, halos, true, true, size, 0);

  // Total volume
  Scal v = 0.;
  for (auto i : m.Cells()) {
    v += m.GetVolume(i);
  }
  PCMP(v, doms.prod());

  // Cell volume
  {
    IdxCell c(0);
    PCMP(m.GetVolume(c), h.prod());
    for (auto i : m.AllCells()) {
      CMP(m.GetVolume(i), m.GetVolume(c));
    }
  }

  // Face area
  for (int q = 0; q < dim; ++q) {
    Dir d(q);
    IdxFace f(0);
    // Find any face with direction d
    for (IdxFace i : m.AllFaces()) {
      if (m.GetDir(i) == d) {
        f = i;
        break;
      }
    }
    assert(m.GetDir(f) == d);

    PCMP(m.GetArea(f), h.prod() / h[q]);
    for (auto i : m.AllFaces()) {
      if (m.GetDir(i) == d) {
        CMP(m.GetArea(i), m.GetArea(f));
      }
    }
  }

  // Number of elements
  auto sh = size + MIdx(halos * 2); // size with halos
  PCMP(m.GetAllBlockCells().size(), sh.prod());
  size_t nf = 0;
  for (int q = 0; q < dim; ++q) {
    auto w = sh;
    ++w[q];
    nf += w.prod();
  }
  PCMP(m.GetAllBlockFaces().size(), nf);
  PCMP(m.GetAllBlockNodes().size(), (sh + MIdx(1)).prod());

  // Distance between centers
  for (auto i : m.Cells()) {
    Vect xi = m.GetCenter(i);
    for (auto n : m.Nci(i)) {
      Dir d(n / 2);
      Scal k = (n % 2 == 0 ? -1. : 1.);
      auto j = m.GetCell(i, n);
      Vect xj = m.GetCenter(j);
      PCMP((xj - xi)[size_t(d)], h[size_t(d)] * k);
    }
  }

  // Index of opposite face
  for (auto d : m.dirs) {
    fassert_equal(m.GetOpposite(2 * d), 2 * d + 1);
    fassert_equal(m.GetOpposite(2 * d + 1), 2 * d);
  }
  fassert_equal(m.GetOpposite(-1), size_t(-1));

  {
    std::cout << "Cell neighbor cells\n";
    const auto& ic = m.GetIndexCells();
    for (auto c : m.Cells()) {
      std::cout << "c=" << ic.GetMIdx(c) << ":";
      for (auto q : m.Nci(c)) {
        std::cout << " " << ic.GetMIdx(m.GetCell(c, q));
      }
      std::cout << std::endl;
      break;
    }
  }

  {
    std::cout << "\nFace neighbor faces\n";
    const auto& index = m.GetIndexFaces();
    auto str = [&](IdxFace f) {
      std::stringstream s;
      s << index.GetMIdx(f) << index.GetDir(f).GetLetter();
      return s.str();
    };
    for (auto d : m.dirs) {
      const MIdx w = m.GetInBlockCells().GetBegin();
      const IdxFace f = index.GetIdx(w, Dir(d));
      std::cout << "f=" << str(f) << ":";
      for (auto q : m.Nci(IdxCell(0))) {
        std::cout << " " << str(m.GetFace(f, q));
      }
      std::cout << std::endl;
    }
  }

  {
    std::cout << "\nCell neighbor cells\n";
    IdxCell c(0);
    for (auto q : m.Nci(c)) {
      IdxFace f = m.GetFace(c, q);
      PCMP(m.GetNci(c, f), q);
    }
  }

  // Comm
  FieldCell<Scal> fc;
  m.Comm(&fc);
  auto p = dynamic_cast<typename M::CommRequestScal*>(m.GetComm()[0].get());
  CMP(p->field, &fc);

  // Stencil indices
  {
    std::cout << "\nStencilIndices" << std::endl;
    auto& indexc = m.GetIndexCells();
    IdxCell c(indexc.GetIdx(m.GetInBlockCells().GetBegin() + MIdx(1)));
    std::cout << "c: " << c.GetRaw() << " " << indexc.GetMIdx(c) << std::endl;
    for (IdxCell cn : m.Stencil(c)) {
      std::cout << cn.GetRaw() << " " << indexc.GetMIdx(cn) << std::endl;
    }
  }
}

void TestMeshIndices() {
  std::cout << __func__ << std::endl;
  const Rect<Vect> dom(Vect(0), Vect(1));
  using M = generic::MeshCartesian<Scal, dim>;
  const MIdx begin(0);
  const MIdx size(2);
  const int halos = 0;
  M m = generic::InitMesh<M>(dom, begin, size, halos, true, true, size, 0);

  {
    std::cout << "\nm.GetAllBlockCells()" << std::endl;
    auto blockca = m.GetAllBlockCells();
    auto indexc = m.GetIndexCells();
    for (auto w : blockca) {
      std::cout << w << " " << indexc.GetIdx(w).GetRaw() << std::endl;
    }
    std::cout << std::endl;
  }

  {
    std::cout << "\nm.GetAllBlockFaces()" << std::endl;
    auto blockfa = m.GetAllBlockFaces();
    auto indexf = m.GetIndexFaces();
    for (auto p : blockfa) {
      std::cout << p.first << " " << p.second.GetLetter() << ","
                << indexf.GetIdx(p).GetRaw() << std::endl;
    }
    std::cout << std::endl;
  }

  {
    std::cout << "\n(MIdx,Dir) <-> IdxFace" << std::endl;
    GBlock<IdxFace, dim> blockf(begin, size);
    GIndex<IdxFace, dim> indexf(begin, size + MIdx(1));
    for (auto p : blockf) {
      auto f = indexf.GetIdx(p);
      auto pp = indexf.GetMIdxDir(f);
      auto ff = indexf.GetIdx(pp);
      std::cout << p.first << "," << p.second.GetLetter() << " " << f.GetRaw()
                << " | " << pp.first << "," << pp.second.GetLetter() << " "
                << ff.GetRaw() << std::endl;
      assert(f == ff);
      assert(p == pp);
    }
    std::cout << std::endl;

    using It = typename GBlock<IdxFace, dim>::iterator;
    It it(&blockf, MIdx(0), Dir(0));
    int count = 0;
    for (it = blockf.begin(); it != blockf.end(); ++it) {
      std::cout << (*it).first << " " << (*it).second.GetLetter() << std::endl;
      ++count;
    }
    assert(count == dim * size.prod() + (MIdx(size.prod()) / size).sum());
  }
}

int main() {
  TestBlock();
  TestMesh();
  TestMeshIndices();
}
