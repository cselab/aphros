/*
 * mesh3d.hpp
 *
 *  Created on: Apr 22, 2016
 *      Author: Petr Karnakov (petr.karnakov@gmail.com)
 */

#pragma once

#include "mesh.hpp"

#define PERX
#define PERZ

namespace geom {

namespace geom3d {


template <class Scal>
class MeshStructured : public MeshGeneric<Scal, 3> {
 public:
  static constexpr size_t dim = 3;
  using Vect = geom::Vect<Scal, dim>;
  using Direction = geom::Direction<dim>;
  using MIdx = geom::Vect<IntIdx, dim>;
  using BlockNodes = geom::BlockNodes<dim>;
  using BlockCells = geom::BlockCells<dim>;
  using BlockFaces = geom::BlockFaces<dim>;
  static constexpr size_t kCellNumNeighbourFaces = 6;
  static constexpr size_t kCellNumNeighbourNodes = 8;
  static constexpr size_t kFaceNumNeighbourNodes = 4;
  static constexpr size_t kFaceNumNeighbourCells = 2;

 private:
  // b:Block, fc:FieldCell, ff:FieldFace, fn:FieldNode
  BlockNodes b_nodes_;
  FieldNode<Vect> fn_node_;
  BlockCells b_cells_;
  BlockFaces b_faces_;
  FieldCell<Vect> fc_center_;
  FieldCell<Scal> fc_volume_;
  FieldFace<Vect> ff_center_;
  FieldFace<Scal> ff_area_;
  FieldFace<Vect> ff_surface_;
  std::array<IntIdx, kCellNumNeighbourFaces> cell_neighbour_cell_offset_;
  FieldCell<std::array<IdxFace, kCellNumNeighbourFaces>> fc_neighbour_face_;
  FieldCell<IdxNode> fc_neighbour_node_base_;
  std::array<IntIdx, kCellNumNeighbourNodes> cell_neighbour_node_offset_;
  FieldFace<Direction> ff_direction_;
  FieldFace<std::array<IdxCell, kFaceNumNeighbourCells>> ff_neighbour_cell_;
  FieldFace<std::array<IdxNode, kFaceNumNeighbourNodes>> ff_neighbour_node_;
  FieldFace<std::array<Vect, kFaceNumNeighbourCells>> ff_to_cell_;
  FieldCell<bool> fc_is_inner_;
  FieldFace<bool> ff_is_inner_;
  FieldCell<bool> fc_is_excluded_;
  FieldFace<bool> ff_is_excluded_;

 public:
  MeshStructured() {}
  MeshStructured(const BlockNodes& b_nodes, const FieldNode<Vect>& fn_node);
  const BlockCells& GetBlockCells() const {
    return b_cells_;
  }
  const BlockFaces& GetBlockFaces() const {
    return b_faces_;
  }
  const BlockNodes& GetBlockNodes() const {
    return b_nodes_;
  }
  Vect GetCenter(IdxCell idx) const override {
    return fc_center_[idx];
  }
  Vect GetCenter(IdxFace idx) const override {
    return ff_center_[idx];
  }
  Vect GetSurface(IdxFace idx) const override {
    return ff_surface_[idx];
  }
  Vect GetNode(IdxNode idx) const override {
    return fn_node_[idx];
  }
  Scal GetVolume(IdxCell idx) const override {
    return fc_volume_[idx];
  }
  Scal GetArea(IdxFace idx) const override {
    return ff_area_[idx];
  }
  size_t GetNumCells() const override {
    return b_cells_.size();
  }
  size_t GetNumFaces() const override {
    return b_faces_.size();
  }
  size_t GetNumNodes() const override {
    return b_nodes_.size();
  }
  // n = 0 .. 7
  IdxCell GetNeighbourCell(IdxCell idx, size_t n) const override {
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
  // n same as for GetNeighbourCell()
  IdxFace GetNeighbourFace(IdxCell idxcell, size_t n) const override {
    return fc_neighbour_face_[idxcell][n];
  }
  Scal GetOutwardFactor(IdxCell, size_t n) const override {
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
  // n same as for GetNeighbourCell()
  Vect GetOutwardSurface(IdxCell idxcell, size_t n) const override {
    return GetSurface(GetNeighbourFace(idxcell, n)) *
        GetOutwardFactor(idxcell, n);
  }
  // n = 0 .. 7
  IdxNode GetNeighbourNode(IdxCell idxcell, size_t n) const override {
    IdxNode idxnode = fc_neighbour_node_base_[idxcell];
    idxnode.AddRaw(cell_neighbour_node_offset_[n]);
    return idxnode;
  }
  Direction GetDirection(IdxFace idx) const {
    return ff_direction_[idx];
  }
  IdxCell GetNeighbourCell(IdxFace idx, size_t n) const override {
    return ff_neighbour_cell_[idx][n];
  }
  Vect GetVectToCell(IdxFace idx, size_t n) const {
    return ff_to_cell_[idx][n];
  }
  IdxNode GetNeighbourNode(IdxFace idx, size_t n) const override {
    return ff_neighbour_node_[idx][n];
  }
  size_t GetNumNeighbourFaces(IdxCell) const override {
    return kCellNumNeighbourFaces;
  }
  size_t GetNumNeighbourNodes(IdxCell) const override {
    return kCellNumNeighbourNodes;
  }
  size_t GetNumNeighbourCells(IdxFace) const override {
    return kFaceNumNeighbourCells;
  }
  size_t GetNumNeighbourNodes(IdxFace) const override {
    return kFaceNumNeighbourNodes;
  }
  bool IsInner(IdxCell idxcell) const override {
    return fc_is_inner_[idxcell];
  }
  bool IsInner(IdxFace idxface) const override {
    return ff_is_inner_[idxface];
  }
  bool IsInside(IdxCell idxcell, Vect vect) const override {
    for (size_t i = 0; i < GetNumNeighbourFaces(idxcell); ++i) {
      IdxFace idxface = GetNeighbourFace(idxcell, i);
      if (GetOutwardSurface(idxcell, i).dot(vect - GetCenter(idxface)) > 0) {
        return false;
      }
    }
    return true;
  }
  void ExcludeCells(const std::vector<IdxCell>& list) {
    for (auto idxcell : list) {
      fc_is_excluded_[idxcell] = true;
      fc_is_inner_[idxcell] = false;
    }
    for (auto idxface : this->Faces()) {
      for (size_t i = 0; i < GetNumNeighbourCells(idxface); ++i) {
        IdxCell idxcell = GetNeighbourCell(idxface, i);
        if (!idxcell.IsNone() && fc_is_excluded_[idxcell]) {
          ff_is_inner_[idxface] = false;
          ff_neighbour_cell_[idxface][i] = IdxCell::None();
        }
      }
    }
    for (auto idxface : this->Faces()) {
      bool found_valid = false;
      for (size_t i = 0; i < GetNumNeighbourCells(idxface); ++i) {
        IdxCell idxcell = GetNeighbourCell(idxface, i);
        if (!idxcell.IsNone()) {
          found_valid = true;
        }
      }
      if (!found_valid) {
        ff_is_excluded_[idxface] = true;
      }
    }
  }
  bool IsExcluded(IdxCell idxcell) const {
    return fc_is_excluded_[idxcell];
  }
  bool IsExcluded(IdxFace idxface) const {
    return ff_is_excluded_[idxface];
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
};

template <class Scal>
MeshStructured<Scal>::MeshStructured(const BlockNodes& b_nodes,
                                     const FieldNode<Vect>& fn_node)
    : b_nodes_(b_nodes)
    , fn_node_(fn_node)
    , b_cells_(b_nodes_.GetBegin(), b_nodes_.GetDimensions() - MIdx(1))
    , b_faces_(b_cells_.GetDimensions())
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
    , fc_is_inner_(b_cells_)
    , ff_is_inner_(b_faces_)
    , fc_is_excluded_(b_cells_, false)
    , ff_is_excluded_(b_faces_, false)
{
  MIdx mb = b_cells_.GetBegin(), me = b_cells_.GetEnd();

  // Base for cell neighbours
  for (auto midx : b_cells_) {
    IdxCell idxcell(b_cells_.GetIdx(midx));
    fc_neighbour_node_base_[idxcell] = b_nodes_.GetIdx(midx);

    auto nface = [this, midx](IntIdx i, IntIdx j, IntIdx k, Direction dir) {
      return b_faces_.GetIdx(midx + MIdx(i, j, k), dir);
    };
    fc_neighbour_face_[idxcell] = {{
        nface(0, 0, 0, Direction::i),
        nface(1, 0, 0, Direction::i),
        nface(0, 0, 0, Direction::j),
        nface(0, 1, 0, Direction::j),
        nface(0, 0, 0, Direction::k),
        nface(0, 0, 1, Direction::k)}};

    fc_is_inner_[idxcell] = (mb < midx && midx + MIdx(1) < me);
    #ifdef PERX
    // adhoc for periodic in x
    if (midx[0] == mb[0] || midx[0] + 1 == me[0]) {
      fc_is_inner_[idxcell] = true;
    }
    #endif
    #ifdef PERZ
    // adhoc for periodic in z
    if (midx[2] == mb[2] || midx[2] + 1 == me[2]) {
      fc_is_inner_[idxcell] = true;
    }
    #endif
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
  for (Direction dir : {Direction::i, Direction::j, Direction::k}) {
    for (auto midx : BlockCells(mb, me - mb + dir)) {
      IdxFace idxface = b_faces_.GetIdx(midx, dir);
      ff_direction_[idxface] = dir;
      ff_neighbour_cell_[idxface] = {{
          b_cells_.GetIdx(midx - dir),
          b_cells_.GetIdx(midx)}};

      auto l = [this, &midx](IntIdx i, IntIdx j, IntIdx k) {
        return b_nodes_.GetIdx(midx + MIdx(i, j, k));
      };
      if (dir == Direction::i) {
        ff_neighbour_node_[idxface] = {{
            l(0,0,0), l(0,1,0), l(0,1,1), l(0,0,1)}};
      } else if (dir == Direction::j) {
        ff_neighbour_node_[idxface] = {{
            l(0,0,0), l(0,0,1), l(1,0,1), l(1,0,0)}};
      } else {
        // dir == Direction::k
        ff_neighbour_node_[idxface] = {{
            l(0,0,0), l(1,0,0), l(1,1,0), l(0,1,0)}};
      }

      if (midx[dir] == mb[dir]) {
        ff_neighbour_cell_[idxface][0] = IdxCell::None();
      }
      if (midx[dir] == me[dir]) {
        ff_neighbour_cell_[idxface][1] = IdxCell::None();
      }
      #ifdef PERX
      // adhoc for periodic in x
      if (midx[0] == mb[0] && dir == Direction::i) {
        ff_neighbour_cell_[idxface][0] = 
            b_cells_.GetIdx(MIdx(me[0] - 1, midx[1], midx[2]));
      }
      if (midx[0] == me[0] && dir == Direction::i) {
        ff_neighbour_cell_[idxface][1] = 
            b_cells_.GetIdx(MIdx(mb[0], midx[1], midx[2]));
      }
      #endif
      #ifdef PERZ
      // adhoc for periodic in x
      if (midx[2] == mb[2] && dir == Direction::k) {
        ff_neighbour_cell_[idxface][0] = 
            b_cells_.GetIdx(MIdx(midx[0], midx[1], me[2] - 1));
      }
      if (midx[2] == me[2] && dir == Direction::k) {
        ff_neighbour_cell_[idxface][1] = 
            b_cells_.GetIdx(MIdx(midx[0], midx[1], mb[2]));
      }
      #endif
    }
  }

  // Face centers, area
  for (auto idx : this->Faces()) {
    ff_center_[idx] = CalcCenter(idx);
    ff_area_[idx] = CalcArea(idx);
    ff_surface_[idx] = CalcSurface(idx);
  }

  // Cell centers, volume
  for (auto idx : this->Cells()) {
    fc_center_[idx] = CalcCenter(idx);
    fc_volume_[idx] = CalcVolume(idx);
  }

  // Vect to cell
  for (Direction dir : {Direction::i, Direction::j, Direction::k}) {
    for (auto midx : BlockCells(mb, me - mb + dir)) {
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
      #ifdef PERX
      // adhoc for periodic in x
      if (midx[0] == mb[0] && dir == Direction::i) {
        auto pf = b_faces_.GetIdx(MIdx(me[0], midx[1], midx[2]), dir);
        auto pc = b_cells_.GetIdx(MIdx(me[0] - 1, midx[1], midx[2]));
        ff_to_cell_[idxface][0] = GetCenter(pc) - GetCenter(pf);
      }
      if (midx[0] == me[0] && dir == Direction::i) {
        auto pf = b_faces_.GetIdx(MIdx(mb[0], midx[1], midx[2]), dir);
        auto pc = b_cells_.GetIdx(MIdx(mb[0], midx[1], midx[2]));
        ff_to_cell_[idxface][1] = GetCenter(pc) - GetCenter(pf);
      }
      #endif
      #ifdef PERZ
      // adhoc for periodic in z
      if (midx[2] == mb[2] && dir == Direction::k) {
        auto pf = b_faces_.GetIdx(MIdx(midx[0], midx[1], me[2]), dir);
        auto pc = b_cells_.GetIdx(MIdx(midx[0], midx[1], me[2] - 1));
        ff_to_cell_[idxface][0] = GetCenter(pc) - GetCenter(pf);
      }
      if (midx[2] == me[2] && dir == Direction::k) {
        auto pf = b_faces_.GetIdx(MIdx(midx[0], midx[1], mb[2]), dir);
        auto pc = b_cells_.GetIdx(MIdx(midx[0], midx[1], mb[2]));
        ff_to_cell_[idxface][1] = GetCenter(pc) - GetCenter(pf);
      }
      #endif
    }
  }

  // Mark inner faces
  for (auto idxface : this->Faces()) {
    ff_is_inner_[idxface] = (!GetNeighbourCell(idxface, 0).IsNone() &&
        !GetNeighbourCell(idxface, 1).IsNone());
  }
}

} // namespace geom2d

} // namespace geom
