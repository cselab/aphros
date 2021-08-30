// Created by Petr Karnakov on 25.04.2018
// Copyright 2018 ETH Zurich

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
#include "loop.h"
#include "map.h"
#include "notation.h"
#include "range.h"
#include "rangein.h"
#include "rangemulti.h"
#include "transform.h"
#include "util/mpi.h"
#include "util/suspender.h"
#include "vect.h"

// Returns column of cells cmm,cm,cp,cpp.
// nfi: 0 or 1 such that m.GetCell(f, nfi) == cp
template <class M>
inline void GetCellColumn(
    const M& m, IdxFace f, Side s, IdxCell& cmm, IdxCell& cm, IdxCell& cp,
    IdxCell& cpp) {
  const size_t d = m.GetIndexFaces().GetDir(f).raw();
  cp = m.GetCell(f, s);
  cm = m.GetCell(f, s.opposite());
  cpp = m.GetCell(cp, IdxNci(2 * d + s.raw()));
  cmm = m.GetCell(cm, IdxNci(2 * d + 1 - s.raw()));
}

template <class Scal_, size_t dim_>
class MeshCartesian {
 public:
  using Scal = Scal_;
  static constexpr size_t dim = dim_;
  static constexpr generic::Range<size_t> dirs{dim};
  using Vect = generic::Vect<Scal, dim>;
  using Dir = GDir<dim>;
  using Direction = generic::Direction<dim>;
  using MIdx = generic::MIdx<dim>;
  using BlockCells = GBlockCells<dim>;
  using BlockFaces = GBlockFaces<dim>;
  using BlockNodes = GBlockNodes<dim>;
  using IndexCells = GIndex<IdxCell, dim>;
  using IndexFaces = GIndex<IdxFace, dim>;
  using IndexNodes = GIndex<IdxNode, dim>;
  using M = MeshCartesian;
  enum class Type { regular, cut, excluded };

  template <class T>
  static constexpr T Pow(T base, size_t pow) {
    return pow == 0 ? 1 : base * Pow(base, pow - 1);
  }

  static constexpr size_t kCellNumNeighborFaces = dim * 2;
  static constexpr size_t kCellNumNeighborNodes = Pow(2, dim);
  static constexpr size_t kFaceNumNeighborNodes = Pow(2, dim - 1);
  static constexpr size_t kFaceNumNeighborCells = 2;
  static constexpr size_t kNumStencil = Pow(3, dim);
  static constexpr size_t kNumStencil5 = Pow(5, dim);
  static constexpr bool kIsEmbed = false;

  using Cell = typename geom::Loop<M>::Cell;
  using Face = typename geom::Loop<M>::Face;

  template <class Func>
  void ForEachCell(Func func) const {
    geom::Loop<M>::ForEachCell(*this, func);
  }

  template <class Func>
  void ForEachFace(Func func) const {
    geom::Loop<M>::ForEachFace(*this, func);
  }

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
    Scal particles_halo_radius = 2; // surplus to the size of the bounding box
                                    // for transfer of particles
                                    // relative to the cell size
    static int GetIdFromBlock(MIdx block, MIdx global_blocks) {
      return GIndex<int, dim>(global_blocks).GetIdx(block);
    };
    int GetIdFromBlock(MIdx block) const {
      return GetIdFromBlock(block, global_blocks);
    }
    static MIdx GetBlockFromId(int id, MIdx global_blocks) {
      return GIndex<int, dim>(global_blocks).GetMIdx(id);
    };
    MIdx GetBlockFromId(int id) const {
      return GetBlockFromId(id, global_blocks);
    }
    MIdx GetBlockFromPoint(Vect x) const {
      MIdx w((x - global_origin) / block_length);
      w = w.max(MIdx(0));
      w = w.min(global_blocks - MIdx(1));
      return w;
    }
    MPI_Comm comm;
  };

  // begin: begin, lower corner cell index
  // size: inner cells size
  // domain: domain, rectangle covering inner cells
  // halos: halo cells from each side
  // isroot: root block
  // gs: global mesh size
  // id: unique id
  MeshCartesian(
      MIdx begin, MIdx size, Rect<Vect> domain, int halos, bool isroot,
      bool islead, MIdx gs, int id);
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
  void SetShared(M* mshared) {
    mshared_ = mshared;
  }
  M& GetShared() {
    return *mshared_;
  }
  const M& GetShared() const {
    return *mshared_;
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
    return GIndex<size_t, dim>(global_size_.max(MIdx(1000))).GetIdx(w);
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
    return ff_center_[f];
  }
  const Vect& GetSurface(IdxFace f) const {
    return face_surface_[indexf_.GetDir(f).raw()];
  }
  Vect GetNode(IdxNode n) const {
    return domain_.low + Vect(indexn_.GetMIdx(n) - incells_begin_) * cell_size_;
  }
  Scal GetVolume(IdxCell) const {
    return cell_volume_;
  }
  Scal GetArea(IdxFace f) const {
    return face_area_[indexf_.GetDir(f).raw()];
  }
  Scal GetAreaFraction(IdxFace) const {
    return 1;
  }
  Scal GetVolumeFraction(IdxCell) const {
    return 1;
  }
  IdxCell GetCell(IdxCell c, IdxNci q) const {
    return IdxCell(c.raw() + cell_cell_[q.raw()]);
  }
  IdxFace GetFace(IdxCell c, IdxNci q) const {
    return IdxFace(c.raw() + cell_face_[q.raw()]);
  }
  IdxFace GetFace(IdxFace f, IdxNci q) const {
    const auto qm = kCellNumNeighborFaces;
    return IdxFace(
        f.raw() + face_face_[indexf_.GetDir(f).raw() * qm + q.raw()]);
  }
  int GetOutwardFactor(IdxCell, IdxNci q) const {
    return -1 + (q.raw() % 2) * 2;
  }
  Vect GetOutwardSurface(IdxCell c, IdxNci q) const {
    return GetSurface(GetFace(c, q)) * GetOutwardFactor(c, q);
  }
  IdxNode GetNode(IdxCell c, size_t q) const {
    assert(q < kCellNumNeighborNodes);
    return IdxNode(c.raw() + cell_node_[q]);
  }
  Dir GetDir(IdxFace f) const {
    return indexf_.GetDir(f);
  }
  IdxCell GetCell(IdxFace f, Side s) const {
    const auto qm = kFaceNumNeighborCells;
    return IdxCell(
        f.raw() + face_cell_[indexf_.GetDir(f).raw() * qm + s.raw()]);
  }
  Vect GetVectToCell(IdxFace f, Side s) const {
    return GetCenter(GetCell(f, s)) - GetCenter(f);
  }
  IdxNode GetNode(IdxFace f, size_t q) const {
    auto qm = kFaceNumNeighborNodes;
    assert(q < qm);
    return IdxNode(f.raw() + face_node_[indexf_.GetDir(f).raw() * qm + q]);
  }
  size_t GetNumFaces(IdxCell) const {
    return kCellNumNeighborFaces;
  }
  // Neighbor cell indices
  generic::Range<IdxNci> Nci(IdxCell c) const {
    return generic::Range<IdxNci>(IdxNci(GetNumFaces(c)));
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
        generic::Range<size_t>(kNumStencil),
        [this, c](size_t i) { return IdxCell(c.raw() + stencil_[i]); });
  }
  const auto& GetStencilOffsets() const {
    return stencil_;
  }
  // Cell indices over 5x5x5 stencil
  auto Stencil5(IdxCell c) const {
    return MakeTransformIterator<IdxCell>(
        generic::Range<size_t>(kNumStencil5),
        [this, c](size_t i) { return IdxCell(c.raw() + stencil5_[i]); });
  }
  const auto& GetStencil5Offsets() const {
    return stencil5_;
  }
  // Returns index `q` such that `m.GetFace(c, q) == f`.
  // -1 if f and c are not neighbors
  IdxNci GetNci(IdxCell c, IdxFace f) const {
    for (auto q : Nci(c)) {
      if (GetFace(c, q) == f) {
        return q;
      }
    }
    return IdxNci(-1);
  }
  IdxNci GetOpposite(IdxNci q) const {
    return IdxNci(q.raw() ^ 1);
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
  std::array<IdxCell, 4> GetCellColumn(IdxFace f, size_t nfi) const {
    const auto d = GetIndexFaces().GetDir(f);
    const IdxCell cm = GetCell(f, 1 - nfi);
    const IdxCell cp = GetCell(f, nfi);
    const IdxCell cmm = GetCell(cm, IdxNci(2 * d.raw() + 1 - nfi));
    const IdxCell cpp = GetCell(cp, IdxNci(2 * d.raw() + nfi));
    return {cmm, cm, cp, cpp};
  }
  // Returns true if cell is not a halo cell.
  bool IsInner(IdxCell c) const {
    MIdx w = indexc_.GetMIdx(c);
    return blockci_.GetBegin() <= w && w < blockci_.GetEnd();
  }
  bool IsInnerPoint(Vect x) {
    return domain_.low <= x && x < domain_.high;
  }
  bool IsBoundary(IdxFace f, /*out*/ size_t& nci) const {
    auto p = GetIndexFaces().GetMIdxDir(f);
    const size_t d = p.second.raw();
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
  void AppendExpr(Expr& sum, const ExprFace& v, IdxNci q) const {
    sum[0] += v[1 - q.raw() % 2];
    sum[1 + q.raw()] += v[q.raw() % 2];
    sum.back() += v.back();
  }
  void AppendExpr(Expr& sum, const ExprFace& v, IdxNci q, IdxCell) const {
    AppendExpr(sum, v, q);
  }
  bool IsInside(IdxCell c, Vect x) const {
    for (auto q : Nci(c)) {
      const IdxFace f = GetFace(c, q);
      if (GetOutwardSurface(c, q).dot(x - GetCenter(f)) > 0) {
        return false;
      }
    }
    return true;
  }
  IdxCell FindNearestCell(Vect x) const {
    Scal dn = std::numeric_limits<Scal>::max(); // sqrdist to nearest
    IdxCell cn(0); // cell nearest
    for (auto c : Cells()) {
      const Scal d = x.sqrdist(GetCenter(c));
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
  auto CellsM() const {
    return MakeTransformIterator<IdxCellMesh<M>>(
        Cells(), [this](IdxCell c) { return IdxCellMesh<M>(c, *this); });
  }
  auto SuCellsM() const {
    return MakeTransformIterator<IdxCellMesh<M>>(
        SuCells(), [this](IdxCell c) { return IdxCellMesh<M>(c, *this); });
  }
  auto AllCellsM() const {
    return MakeTransformIterator<IdxCellMesh<M>>(
        AllCells(), [this](IdxCell c) { return IdxCellMesh<M>(c, *this); });
  }
  auto FacesM() const {
    return MakeTransformIterator<IdxFaceMesh<M>>(
        Faces(), [this](IdxFace f) { return IdxFaceMesh<M>(f, *this); });
  }
  auto SuFacesM() const {
    return MakeTransformIterator<IdxFaceMesh<M>>(
        SuFaces(), [this](IdxFace f) { return IdxFaceMesh<M>(f, *this); });
  }
  auto AllFacesM() const {
    return MakeTransformIterator<IdxFaceMesh<M>>(
        AllFaces(), [this](IdxFace f) { return IdxFaceMesh<M>(f, *this); });
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

  enum class CommStencil { none, full_two, full_one, direct_two, direct_one };
  struct CommRequest {
    CommRequest() = default;
    CommRequest(CommStencil stencil_) : stencil(stencil_) {}
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
    CommStencil stencil = CommStencil::full_two;
  };
  // FieldCell<Scal>
  struct CommRequestScal : public CommRequest {
    CommRequestScal(FieldCell<Scal>* field_, CommStencil stencil_)
        : CommRequest(stencil_), field(field_) {}
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
  struct CommRequestVect : public CommRequest {
    // f: vector field
    // d_: component (0,1,2), or -1 for all
    CommRequestVect(FieldCell<Vect>* field_, int d_, CommStencil stencil_)
        : CommRequest(stencil_), field(field_), d(d_) {}
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
  void Comm(
      FieldCell<Scal>* field, CommStencil stencil = CommStencil::full_two);
  void Comm(
      FieldCell<Vect>* field, int component,
      CommStencil stencil = CommStencil::full_two);
  void Comm(
      FieldCell<Vect>* field, CommStencil stencil = CommStencil::full_two);
  const std::vector<std::unique_ptr<CommRequest>>& GetComm() const;
  void ClearComm();

  struct CommPartRequest {
    std::vector<Vect>* x = nullptr; // particle position
    std::vector<bool>* is_inner = nullptr; // true if particle is owned by block
    std::vector<std::vector<Scal>*> attr_scal; // scalar attributes
    std::vector<std::vector<Vect>*> attr_vect; // vector attributes
  };
  // Requests exchange of particles between blocks.
  // Particles with `is_inner` set to `true` are transfered to blocks
  // within `flags.particles_halo_radius`.
  // Then `is_inner` is set to `true` for particles with `IsInnerPoint()` true.
  void CommPart(const CommPartRequest&);
  const std::vector<CommPartRequest>& GetCommPart() const;
  void ClearCommPart();

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
  /// Sets implementation of GetMpiRankFromBlockId
  void SetHandlerMpiRankFromId(std::function<int(int)>);
  /// Returns MPI rank that owns block for which `GetId()` returns `id`
  int GetMpiRankFromId(int id) const;
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
  void Reduce(std::string* buf, ReductionType::Concat) {
    Reduce(std::make_unique<typename UReduce<Scal>::OpCatString>(buf));
  }
  template <class T>
  void Reduce(std::vector<std::vector<T>>* buf, ReductionType::Concat) {
    Reduce(std::make_unique<typename UReduce<Scal>::template OpCatVT<T>>(buf));
  }
  void Reduce(std::vector<std::string>* buf, ReductionType::Concat) {
    Reduce(std::make_unique<typename UReduce<Scal>::OpCatVectorString>(buf));
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
  void Bcast(T* elem) {
    Bcast(std::make_unique<typename UReduce<Scal>::template OpCatRaw<T>>(elem));
  }
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
  void BcastFromLead(std::unique_ptr<Op>&& o);
  template <class T>
  void BcastFromLead(T* elem) {
    BcastFromLead(
        std::make_unique<typename UReduce<Scal>::template OpCatRaw<T>>(elem));
  }
  const std::vector<std::unique_ptr<Op>>& GetBcastFromLead() const;
  void ClearBcastFromLead();
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

  M* mshared_ = nullptr;
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
  std::array<size_t, kCellNumNeighborNodes> cell_node_;
  std::array<size_t, kFaceNumNeighborCells * dim> face_cell_;
  std::array<size_t, kCellNumNeighborFaces * dim> face_face_;
  std::array<size_t, kFaceNumNeighborNodes * dim> face_node_;
  std::array<size_t, kNumStencil> stencil_; // 3x3x3 stencil
  std::array<size_t, kNumStencil5> stencil5_; // 5x5x5 stencil

  FieldCell<Vect> fc_center_;
  FieldFace<Vect> ff_center_;
  Suspender susp_;
  std::string timer_report_path_;
  struct Imp;
  std::unique_ptr<Imp> imp;
};
