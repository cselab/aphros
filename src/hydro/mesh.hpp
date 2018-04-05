#pragma once

#include <vector>
#include <array>
#include <cassert>
#include <utility>

#include "suspender.h"
#include "idx.h"
#include "vect.hpp"
#include "range.h"
#include "field.h"
#include "map.h"
#include "blockface.h"
#include "rangein.h"

namespace geom {

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
  using BlockNodes = GBlockNodes<dim>;
  using BlockCells = GBlockCells<dim>;
  using BlockFaces = GBlockFaces<dim>;
  static constexpr size_t kCellNumNeighbourFaces = 2 * dim;
  static constexpr size_t kCellNumNeighbourNodes = std::pow(2, dim);
  static constexpr size_t kFaceNumNeighbourNodes = std::pow(2, dim - 1);
  static constexpr size_t kFaceNumNeighbourCells = 2;

 private:
  // b:Block, fc:FieldCell, ff:FieldFace, fn:FieldNode
  BlockNodes b_nodes_;
  FieldNode<Vect> fn_node_;
  BlockCells b_cells_;
  BlockFaces b_faces_;
  // inner 
  BlockCells b_incells_;
  BlockFaces b_infaces_;
  BlockNodes b_innodes_;
  // support
  BlockCells b_sucells_;
  BlockFaces b_sufaces_;
  BlockNodes b_sunodes_;


  FieldCell<Vect> fc_center_;
  FieldCell<Scal> fc_volume_;
  FieldFace<Vect> ff_center_;
  FieldFace<Scal> ff_area_;
  FieldFace<Vect> ff_surface_;
  std::array<IntIdx, kCellNumNeighbourFaces> cell_neighbour_cell_offset_;
  FieldCell<std::array<IdxFace, kCellNumNeighbourFaces>> fc_neighbour_face_;
  FieldCell<IdxNode> fc_neighbour_node_base_;
  std::array<IntIdx, kCellNumNeighbourNodes> cell_neighbour_node_offset_;
  FieldFace<Dir> ff_direction_;
  FieldFace<std::array<IdxCell, kFaceNumNeighbourCells>> ff_neighbour_cell_;
  FieldFace<std::array<IdxNode, kFaceNumNeighbourNodes>> ff_neighbour_node_;
  FieldFace<std::array<Vect, kFaceNumNeighbourCells>> ff_to_cell_;
  FieldCell<bool> fc_is_inner_;
  FieldFace<bool> ff_is_inner_;

 public:
  MeshStructured(const BlockNodes& b_nodes, 
                 const FieldNode<Vect>& fn_node, 
                 int hl /*halos*/);
  // TODO: rename to GetBlockAllCells()
  const BlockCells& GetBlockCells() const {
    return b_cells_;
  }
  const BlockFaces& GetBlockFaces() const {
    return b_faces_;
  }
  const BlockNodes& GetBlockNodes() const {
    return b_nodes_;
  }
  // TODO: rename to GetBlockCells()
  const BlockCells& GetInBlockCells() const {
    return b_incells_;
  }
  const BlockFaces& GetInBlockFaces() const {
    return b_infaces_;
  }
  const BlockNodes& GetInBlockNodes() const {
    return b_innodes_;
  }
  // TODO: rename to GetBlockCells()
  const BlockCells& GetSuBlockCells() const {
    return b_sucells_;
  }
  const BlockFaces& GetSuBlockFaces() const {
    return b_sufaces_;
  }
  const BlockNodes& GetSuBlockNodes() const {
    return b_sunodes_;
  }
  Vect GetCenter(IdxCell idx) const {
    return fc_center_[idx];
  }
  Vect GetCenter(IdxFace idx) const {
    return ff_center_[idx];
  }
  Vect GetSurface(IdxFace idx) const {
    return ff_surface_[idx];
  }
  Vect GetNode(IdxNode idx) const {
    return fn_node_[idx];
  }
  Scal GetVolume(IdxCell idx) const {
    return fc_volume_[idx];
  }
  Scal GetArea(IdxFace idx) const {
    return ff_area_[idx];
  }
  size_t GetNumCells() const {
    return b_cells_.size();
  }
  size_t GetNumFaces() const {
    return b_faces_.size();
  }
  size_t GetNumNodes() const {
    return b_nodes_.size();
  }
  IdxCell GetNeighbourCell(IdxCell idx, size_t n) const {
    // TODO: 3d specific
    assert(n < kCellNumNeighbourFaces);
    if (fc_is_inner_[idx]) {
      idx.AddRaw(cell_neighbour_cell_offset_[n]);
    } else {
      MIdx base = b_cells_.GetMIdx(idx);
      MIdx offset = (std::vector<MIdx>{
        MIdx{-1, 0, 0}, MIdx{1, 0, 0},
        MIdx{0, -1, 0}, MIdx{0, 1, 0},
        MIdx{0, 0, -1}, MIdx{0, 0, 1}})[n];
      if (b_cells_.IsInside(base + offset)) {
        idx.AddRaw(cell_neighbour_cell_offset_[n]);
      } else {
        idx = IdxCell::None();
      }
    }
    return idx;
  }
  IdxFace GetNeighbourFace(IdxCell idxcell, size_t n) const {
    assert(n < kCellNumNeighbourFaces);
    return fc_neighbour_face_[idxcell][n];
  }
  Scal GetOutwardFactor(IdxCell, size_t n) const {
    // TODO: <= 3d specific, maybe replace with odd/even convention
    assert(n < kCellNumNeighbourFaces);
    Scal factor;
    switch (n) {
      case 0:
      case 2:
      case 4:
        factor = -1.;
        break;
      default: // n == 1, 3, 5
        factor = 1.;
    }
    return factor;
  }
  Vect GetOutwardSurface(IdxCell idxcell, size_t n) const {
    assert(n < kCellNumNeighbourFaces);
    return GetSurface(GetNeighbourFace(idxcell, n)) *
        GetOutwardFactor(idxcell, n);
  }
  IdxNode GetNeighbourNode(IdxCell idxcell, size_t n) const {
    assert(n < kCellNumNeighbourNodes);
    IdxNode idxnode = fc_neighbour_node_base_[idxcell];
    idxnode.AddRaw(cell_neighbour_node_offset_[n]);
    return idxnode;
  }
  Dir GetDir(IdxFace idx) const {
    return ff_direction_[idx];
  }
  IdxCell GetNeighbourCell(IdxFace idx, size_t n) const {
    assert(n < kFaceNumNeighbourCells);
    return ff_neighbour_cell_[idx][n];
  }
  Vect GetVectToCell(IdxFace idx, size_t n) const {
    assert(n < kFaceNumNeighbourCells);
    return ff_to_cell_[idx][n];
  }
  IdxNode GetNeighbourNode(IdxFace idx, size_t n) const {
    assert(n < kFaceNumNeighbourNodes);
    return ff_neighbour_node_[idx][n];
  }
  size_t GetNumNeighbourFaces(IdxCell) const {
    return kCellNumNeighbourFaces;
  }
  // Neighbour cell indices
  GRange<size_t> Nci(IdxCell c) const {
    return GRange<size_t>(0, GetNumNeighbourFaces(c));
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
  bool IsInner(IdxCell idxcell) const {
    return fc_is_inner_[idxcell];
  }
  bool IsInner(IdxFace idxface) const {
    return ff_is_inner_[idxface];
  }

  // Returns range of all indices
  template <class Idx>
  GRange<Idx> GetAll() const { 
    return GetAll((Idx*)0);
  }
  GRange<IdxCell> GetAll(IdxCell*) const { 
    return GRange<IdxCell>(0, GetNumCells());
  }
  GRange<IdxFace> GetAll(IdxFace*) const { 
    return GRange<IdxFace>(0, GetNumFaces());
  }
  GRange<IdxNode> GetAll(IdxNode*) const { 
    return GRange<IdxNode>(0, GetNumNodes());
  }
  GRange<IdxCell> AllCells() const {
    return GetAll<IdxCell>();
  }
  GRange<IdxFace> AllFaces() const {
    return GetAll<IdxFace>();
  }
  GRange<IdxNode> AllNodes() const {
    return GetAll<IdxNode>();
  }
  // Type-cast to GRange required for GField initialization
  template <class Idx>
  operator GRange<Idx>() const {
    return GetAll<Idx>(); 
  }

  // Returns range of inner indices
  template <class Idx>
  GRangeIn<Idx, dim> Get() const { 
    return Get((Idx*)0);
  }
  GRangeIn<IdxCell, dim> Get(IdxCell*) const {
    return GRangeIn<IdxCell, dim>(GetBlockCells(), GetInBlockCells());
  }
  GRangeIn<IdxFace, dim> Get(IdxFace*) const {
    return GRangeIn<IdxFace, dim>(GetBlockFaces(), GetInBlockFaces());
  }
  GRangeIn<IdxNode, dim> Get(IdxNode*) const {
    return GRangeIn<IdxNode, dim>(GetBlockNodes(), GetInBlockNodes());
  }
  GRangeIn<IdxCell, dim> Cells() const {
    return Get<IdxCell>();
  }
  GRangeIn<IdxFace, dim> Faces() const {
    return Get<IdxFace>();
  }
  GRangeIn<IdxNode, dim> Nodes() const {
    return Get<IdxNode>();
  }

  // Returns range of support indices
  template <class Idx>
  GRangeIn<Idx, dim> GetSu() const { 
    return GetSu((Idx*)0);
  }
  GRangeIn<IdxCell, dim> GetSu(IdxCell*) const {
    return GRangeIn<IdxCell, dim>(GetBlockCells(), GetSuBlockCells());
  }
  GRangeIn<IdxFace, dim> GetSu(IdxFace*) const {
    return GRangeIn<IdxFace, dim>(GetBlockFaces(), GetSuBlockFaces());
  }
  GRangeIn<IdxNode, dim> GetSu(IdxNode*) const {
    return GRangeIn<IdxNode, dim>(GetBlockNodes(), GetSuBlockNodes());
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

  bool IsInside(IdxCell idxcell, Vect vect) const {
    for (size_t i = 0; i < GetNumNeighbourFaces(idxcell); ++i) {
      IdxFace idxface = GetNeighbourFace(idxcell, i);
      if (GetOutwardSurface(idxcell, i).dot(vect - GetCenter(idxface)) > 0) {
        return false;
      }
    }
    return true;
  }
  IdxCell FindNearestCell(Vect point) const {
    IdxCell idxcell_nearest = IdxCell(0);
    for (auto idxcell : Cells()) {
      if (point.dist(GetCenter(idxcell)) <
          point.dist(GetCenter(idxcell_nearest))) {
        idxcell_nearest = idxcell;
      }
    }
    return idxcell_nearest;
  }
  Vect GetNormal(IdxFace idxface) const {
    return GetSurface(idxface) / GetArea(idxface);
  }

 private:
   Vect CalcCenter(IdxCell idxcell) const {
     Vect res = Vect::kZero;
     for (size_t i = 0; i < GetNumNeighbourNodes(idxcell); ++i) {
       res += GetNode(GetNeighbourNode(idxcell, i));
     }
     return res / Scal(GetNumNeighbourNodes(idxcell));
   }
   // TODO: Face center using weighted sum
   Vect CalcCenter(IdxFace idxface) const {
     Vect res = Vect::kZero;
     for (size_t i = 0; i < GetNumNeighbourNodes(idxface); ++i) {
       res += GetNode(GetNeighbourNode(idxface, i));
     }
     return res / Scal(GetNumNeighbourNodes(idxface));
   }
   static Vect CalcTriangleSurface(Vect a, Vect b, Vect c) {
     return (b - a).cross(c - a) * 0.5;
   }
   static Vect CalcTriangleCenter(Vect a, Vect b, Vect c) {
     return (a + b + c) / 3.;
   }
   Vect CalcSurface(IdxFace idxface) const {
     Vect res = Vect::kZero;

     for (size_t i = 2; i < GetNumNeighbourNodes(idxface); ++i) {
       Vect a = GetNode(GetNeighbourNode(idxface, 0));
       Vect b = GetNode(GetNeighbourNode(idxface, i - 1));
       Vect c = GetNode(GetNeighbourNode(idxface, i));
       res += CalcTriangleSurface(a, b, c);
     }

     return res;
   }
   Scal CalcArea(IdxFace idxface) const {
     return CalcSurface(idxface).norm();
   }
   Scal CalcVolume(IdxCell idxcell) const {
     Scal res = 0.;
     for (size_t i = 0; i < GetNumNeighbourFaces(idxcell); ++i) {
       IdxFace idxface = GetNeighbourFace(idxcell, i);
       res += GetCenter(idxface)[0] * GetOutwardSurface(idxcell, i)[0];
     }
     return res;
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
  void Comm(FieldCell<Scal>* u) {
    vcm_.push_back(u);
  }
  void Dump(FieldCell<Scal>* u, std::string name) {
    vd_.push_back(std::make_pair(u, name));
  }
  void Reduce(Scal* u, std::string op) {
    vrd_.push_back(std::make_pair(u, op));
  }
  void Solve(const LS& ls) {
    vls_.push_back(ls);
  }
  const std::vector<FieldCell<Scal>*>& GetComm() const {
    return vcm_;
  }
  const std::vector<std::pair<FieldCell<Scal>*, std::string>>& GetDump() const {
    return vd_;
  }
  void ClearComm() {
    vcm_.clear();
  }
  void ClearDump() {
    vd_.clear();
  }
  const std::vector<std::pair<Scal*, std::string>>& GetReduce() const {
    return vrd_;
  }
  void ClearReduce() {
    vrd_.clear();
  }
  const std::vector<LS>& GetSolve() const {
    return vls_;
  }
  void ClearSolve() {
    vls_.clear();
  }

 private:
  Suspender susp_;
  std::vector<FieldCell<Scal>*> vcm_; // fields for [c]o[m]munication
  std::vector<std::pair<FieldCell<Scal>*, std::string>> vd_; // fields for dump
  std::vector<std::pair<Scal*, std::string>> vrd_; // scalars for reduce
  std::vector<LS> vls_; // linear system
  // END DISTR
};


template <class _Scal, size_t _dim>
MeshStructured<_Scal, _dim>::MeshStructured(
    const BlockNodes& b_nodes,
    const FieldNode<Vect>& fn_node,
    int hl /*halos*/)
    : b_nodes_(b_nodes)
    , fn_node_(fn_node)
    , b_cells_(b_nodes_.GetBegin(), b_nodes_.GetDimensions() - MIdx(1))
    , b_faces_(b_nodes_.GetBegin(), b_cells_.GetDimensions())
    , b_incells_(b_cells_.GetBegin() + MIdx(hl), 
                 b_cells_.GetDimensions() - MIdx(2*hl))
    , b_infaces_(b_cells_.GetBegin() + MIdx(hl), 
                 b_cells_.GetDimensions() - MIdx(2*hl))
    , b_innodes_(b_nodes_.GetBegin() + MIdx(hl), 
                 b_nodes_.GetDimensions() - MIdx(2*hl))
    , b_sucells_(b_cells_.GetBegin() + MIdx(1), 
                 b_cells_.GetDimensions() - MIdx(2))
    , b_sufaces_(b_cells_.GetBegin() + MIdx(1), 
                 b_cells_.GetDimensions() - MIdx(2))
    , b_sunodes_(b_nodes_.GetBegin() + MIdx(1), 
                 b_nodes_.GetDimensions() - MIdx(2))
    , fc_center_(b_cells_)
    , fc_volume_(b_cells_)
    , ff_center_(b_faces_)
    , ff_area_(b_faces_)
    , ff_surface_(b_faces_)
    , fc_neighbour_face_(b_cells_)
    , fc_neighbour_node_base_(b_cells_)
    , ff_direction_(b_faces_)
    , ff_neighbour_cell_(b_faces_)
    , ff_neighbour_node_(b_faces_)
    , ff_to_cell_(b_faces_)
    , fc_is_inner_(b_cells_, false)
    , ff_is_inner_(b_faces_, false)
{
  static_assert(dim == 3, "Not implemented for dim != 3");
  MIdx mb = b_cells_.GetBegin(), me = b_cells_.GetEnd();

  // Base for cell neighbours
  for (auto midx : b_cells_) {
    IdxCell idxcell(b_cells_.GetIdx(midx));
    fc_neighbour_node_base_[idxcell] = b_nodes_.GetIdx(midx);

    auto nface = [this, midx](IntIdx i, IntIdx j, IntIdx k, Dir dir) {
      return b_faces_.GetIdx(midx + MIdx(i, j, k), dir);
    };
    fc_neighbour_face_[idxcell] = {{
        nface(0, 0, 0, Dir::i),
        nface(1, 0, 0, Dir::i),
        nface(0, 0, 0, Dir::j),
        nface(0, 1, 0, Dir::j),
        nface(0, 0, 0, Dir::k),
        nface(0, 0, 1, Dir::k)}};
  }

  // Mark inner cells
  for (auto i : this->Cells()) {
    fc_is_inner_[i] = true;
  }

  // Offset for cell neighbour cells
  {
    auto offset = [this](IntIdx i, IntIdx j, IntIdx k) {
      return static_cast<IntIdx>(b_cells_.GetIdx(MIdx(i, j, k)).GetRaw() -
          b_cells_.GetIdx(MIdx(0, 0, 0)).GetRaw());
    };

    cell_neighbour_cell_offset_ = {{offset(-1, 0, 0), offset(1, 0, 0),
                                  offset(0, -1, 0), offset(0, 1, 0),
                                  offset(0, 0, -1), offset(0, 0, 1)}};
  }

  // Offset for cell neighbour nodes
  {
    auto offset = [this](IntIdx i, IntIdx j, IntIdx k) {
      return static_cast<IntIdx>(b_nodes_.GetIdx(MIdx(i, j, k)).GetRaw() -
          b_nodes_.GetIdx(MIdx(0, 0, 0)).GetRaw());
    };

    cell_neighbour_node_offset_ = {{offset(0, 0, 0), offset(1, 0, 0),
                                  offset(0, 1, 0), offset(1, 1, 0),
                                  offset(0, 0, 1), offset(1, 0, 1),
                                  offset(0, 1, 1), offset(1, 1, 1)}};
  }

  // Base for i-face and j-face neighbours
  for (Dir dir : {Dir::i, Dir::j, Dir::k}) {
    for (auto midx : BlockCells(mb, me - mb + MIdx(dir))) {
      IdxFace idxface = b_faces_.GetIdx(midx, dir);
      ff_direction_[idxface] = dir;
      ff_neighbour_cell_[idxface] = {{
          b_cells_.GetIdx(midx - MIdx(dir)),
          b_cells_.GetIdx(midx)}};

      auto l = [this, &midx](IntIdx i, IntIdx j, IntIdx k) {
        return b_nodes_.GetIdx(midx + MIdx(i, j, k));
      };
      if (dir == Dir::i) {
        ff_neighbour_node_[idxface] = {{
            l(0,0,0), l(0,1,0), l(0,1,1), l(0,0,1)}};
      } else if (dir == Dir::j) {
        ff_neighbour_node_[idxface] = {{
            l(0,0,0), l(0,0,1), l(1,0,1), l(1,0,0)}};
      } else {
        // dir == Dir::k
        ff_neighbour_node_[idxface] = {{
            l(0,0,0), l(1,0,0), l(1,1,0), l(0,1,0)}};
      }

      const size_t sd(dir);
      if (midx[sd] == mb[sd]) {
        ff_neighbour_cell_[idxface][0] = IdxCell::None();
      }
      if (midx[sd] == me[sd]) {
        ff_neighbour_cell_[idxface][1] = IdxCell::None();
      }
    }
  }

  // Face centers, area
  for (auto idx : this->AllFaces()) {
    ff_center_[idx] = CalcCenter(idx);
    ff_area_[idx] = CalcArea(idx);
    ff_surface_[idx] = CalcSurface(idx);
  }

  // Cell centers, volume
  for (auto idx : this->AllCells()) {
    fc_center_[idx] = CalcCenter(idx);
    fc_volume_[idx] = CalcVolume(idx);
  }

  // Vect to cell
  for (Dir dir : {Dir::i, Dir::j, Dir::k}) {
    for (auto midx : BlockCells(mb, me - mb + MIdx(dir))) {
      IdxFace idxface = b_faces_.GetIdx(midx, dir);
      IdxCell c;
      c = GetNeighbourCell(idxface, 0);
      if (!c.IsNone()) {
        ff_to_cell_[idxface][0] = GetCenter(c) - GetCenter(idxface);
      }
      c = GetNeighbourCell(idxface, 1);
      if (!c.IsNone()) {
        ff_to_cell_[idxface][1] = GetCenter(c) - GetCenter(idxface);
      }
    }
  }

  // Mark inner faces
  for (auto idxface : this->Faces()) {
    ff_is_inner_[idxface] = true; 
  }
}


// Create uniform mesh
// domain    - rectangle covering inner cells
// begin     - index of lower inner cells
// s        - number of inner cells in each direction
// hl        - number of halo layers
template <class Mesh>
Mesh InitUniformMesh(Rect<typename Mesh::Vect> domain,
                     typename Mesh::MIdx begin,
                     typename Mesh::MIdx s, int hl) {
  using Vect = typename Mesh::Vect;
  using MIdx = typename Mesh::MIdx;
  typename Mesh::BlockNodes b_nodes(begin - MIdx(hl), s + MIdx(1 + 2*hl));
  FieldNode<Vect> fn_node(b_nodes);
  for (auto midx : b_nodes) {
    IdxNode idxnode = b_nodes.GetIdx(midx);
    Vect unit = Vect(midx - b_nodes.GetBegin() - MIdx(hl)) / Vect(s);
    fn_node[idxnode] = domain.lb + unit * domain.GetDimensions();
  }
  return Mesh(b_nodes, fn_node, hl);
}

template <class Mesh>
Mesh InitUniformMesh(const Rect<typename Mesh::Vect>& domain,
                     typename Mesh::MIdx mesh_size) {
  using MIdx = typename Mesh::MIdx;
  return InitUniformMesh<Mesh>(domain, MIdx(0), mesh_size);
}


} // namespace geom
