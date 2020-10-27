// Created by Petr Karnakov on 25.04.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <mpi.h>
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
#include "notation.h"
#include "range.h"
#include "rangein.h"
#include "transform.h"
#include "util/suspender.h"
#include "vect.h"

// Returns column of cells cmm,cm,cp,cpp.
// nci: 0 or 1 such that m.GetCell(f, nci) == cp
template <class M>
inline void GetCellColumn(
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
  using Direction = generic::Direction<Scal, dim>;
  using MIdx = GMIdx<dim>;
  using BlockCells = GBlockCells<dim>;
  using BlockFaces = GBlockFaces<dim>;
  using BlockNodes = GBlockNodes<dim>;
  using IndexCells = GIndex<IdxCell, dim>;
  using IndexFaces = GIndex<IdxFace, dim>;
  using IndexNodes = GIndex<IdxNode, dim>;
  using M = MeshStructured;
  enum class Type { regular, cut, excluded };

  static constexpr size_t kCellNumNeighbourFaces = 6;
  static constexpr size_t kCellNumNeighbourNodes = 8;
  static constexpr size_t kFaceNumNeighbourNodes = 4;
  static constexpr size_t kFaceNumNeighbourCells = 2;
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
    std::array<bool, dim> is_periodic = {false, false, false};
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
  MeshStructured(
      MIdx b, MIdx cs, Rect<Vect> domain, int halos, bool isroot, bool islead,
      MIdx gs, int id);
  MeshStructured(const MeshStructured&) = delete;
  MeshStructured(MeshStructured&&);
  MeshStructured& operator=(const MeshStructured&) = delete;
  MeshStructured& operator=(MeshStructured&&) = delete;
  ~MeshStructured();
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
    return (domain_.low + half_cell_size_) +
           Vect(indexc_.GetMIdx(c) - incells_begin_) * cell_size_;
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
    assert(q < kCellNumNeighbourFaces);
    return IdxCell(size_t(c) + cell_cell_[q]);
  }
  IdxFace GetFace(IdxCell c, size_t q) const {
    assert(q < kCellNumNeighbourFaces);
    return IdxFace(size_t(c) + cell_face_[q]);
  }
  IdxFace GetFace(IdxFace f, size_t q) const {
    const auto qm = kCellNumNeighbourFaces;
    assert(q < qm);
    return IdxFace(size_t(f) + face_face_[size_t(indexf_.GetDir(f)) * qm + q]);
  }
  Scal GetOutwardFactor(IdxCell, size_t q) const {
    assert(q < kCellNumNeighbourFaces);
    return cell_outward_[q];
  }
  Vect GetOutwardSurface(IdxCell c, size_t n) const {
    assert(n < kCellNumNeighbourFaces);
    return GetSurface(GetFace(c, n)) * GetOutwardFactor(c, n);
  }
  IdxNode GetNode(IdxCell c, size_t q) const {
    assert(q < kCellNumNeighbourNodes);
    return IdxNode(size_t(c) + cell_node_[q]);
  }
  Dir GetDir(IdxFace f) const {
    return indexf_.GetDir(f);
  }
  IdxCell GetCell(IdxFace f, size_t q) const {
    auto qm = kFaceNumNeighbourCells;
    assert(q < qm);
    return IdxCell(size_t(f) + face_cell_[size_t(indexf_.GetDir(f)) * qm + q]);
  }
  Vect GetVectToCell(IdxFace f, size_t n) const {
    assert(n < kFaceNumNeighbourCells);
    return GetCenter(GetCell(f, n)) - GetCenter(f);
  }
  IdxNode GetNode(IdxFace f, size_t q) const {
    auto qm = kFaceNumNeighbourNodes;
    assert(q < qm);
    size_t d(indexf_.GetDir(f));
    return IdxNode(size_t(f) + face_node_[d * qm + q]);
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
  // Cell indices over 5x5x5 stencil
  auto Stencil5(IdxCell c) const {
    return MakeTransformIterator<IdxCell>(
        GRange<size_t>(0, kNumStencil5),
        [this, c](size_t i) { return IdxCell(size_t(c) + stencil5_[i]); });
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
  MeshStructured& GetMesh() {
    return *this;
  }
  const MeshStructured& GetMesh() const {
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
  std::array<size_t, kCellNumNeighbourFaces> cell_cell_;
  std::array<size_t, kCellNumNeighbourFaces> cell_face_;
  std::array<Scal, kCellNumNeighbourFaces> cell_outward_;
  std::array<size_t, kCellNumNeighbourNodes> cell_node_;
  std::array<size_t, kFaceNumNeighbourCells * dim> face_cell_;
  std::array<size_t, kCellNumNeighbourFaces * dim> face_face_;
  std::array<size_t, kFaceNumNeighbourNodes * dim> face_node_;
  std::array<size_t, kNumStencil> stencil_; // 3x3x3 stencil
  std::array<size_t, kNumStencil5> stencil5_; // 5x5x5 stencil

  Suspender susp_;
  std::string timer_report_path_;
  struct Imp;
  std::unique_ptr<Imp> imp;
};

// Create uniform mesh
// domain: rectangle covering inner cells
// begin: index of lower inner cells
// s: number of inner cells in each direction
// halos: number of halo layers
// isroot: root block
// islead: lead block
// gs: global mesh size
// id: unique id
template <class M>
M InitUniformMesh(
    Rect<typename M::Vect> domain, typename M::MIdx begin, typename M::MIdx s,
    int halos, bool isroot, bool islead, typename M::MIdx gs, int id);
