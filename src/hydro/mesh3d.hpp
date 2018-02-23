/*
 * mesh3d.hpp
 *
 *  Created on: Apr 22, 2016
 *      Author: Petr Karnakov (petr.karnakov@gmail.com)
 */

#pragma once

#include "mesh.hpp"

//#define PERX
//#define PERZ

namespace geom3d {

template <class Scal>
using MeshStructured = geom::MeshStructured<Scal, 3>;

} // namespace geom3d

namespace geom {

// specialization 
//template <class _Scal>
template <class _Scal>
MeshStructured<_Scal, 3>::MeshStructured(
    const geom::BlockNodes<_Scal>& b_nodes,
    const geom::FieldNode<Vect>& fn_node,
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
    for (auto midx : BlockCells(mb, me - mb + dir)) {
      IdxFace idxface = b_faces_.GetIdx(midx, dir);
      ff_direction_[idxface] = dir;
      ff_neighbour_cell_[idxface] = {{
          b_cells_.GetIdx(midx - dir),
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

      if (midx[dir] == mb[dir]) {
        ff_neighbour_cell_[idxface][0] = IdxCell::None();
      }
      if (midx[dir] == me[dir]) {
        ff_neighbour_cell_[idxface][1] = IdxCell::None();
      }
      #ifdef PERX
      // adhoc for periodic in x
      if (midx[0] == mb[0] && dir == Dir::i) {
        ff_neighbour_cell_[idxface][0] = 
            b_cells_.GetIdx(MIdx(me[0] - 1, midx[1], midx[2]));
      }
      if (midx[0] == me[0] && dir == Dir::i) {
        ff_neighbour_cell_[idxface][1] = 
            b_cells_.GetIdx(MIdx(mb[0], midx[1], midx[2]));
      }
      #endif
      #ifdef PERZ
      // adhoc for periodic in x
      if (midx[2] == mb[2] && dir == Dir::k) {
        ff_neighbour_cell_[idxface][0] = 
            b_cells_.GetIdx(MIdx(midx[0], midx[1], me[2] - 1));
      }
      if (midx[2] == me[2] && dir == Dir::k) {
        ff_neighbour_cell_[idxface][1] = 
            b_cells_.GetIdx(MIdx(midx[0], midx[1], mb[2]));
      }
      #endif
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
      if (midx[0] == mb[0] && dir == Dir::i) {
        auto pf = b_faces_.GetIdx(MIdx(me[0], midx[1], midx[2]), dir);
        auto pc = b_cells_.GetIdx(MIdx(me[0] - 1, midx[1], midx[2]));
        ff_to_cell_[idxface][0] = GetCenter(pc) - GetCenter(pf);
      }
      if (midx[0] == me[0] && dir == Dir::i) {
        auto pf = b_faces_.GetIdx(MIdx(mb[0], midx[1], midx[2]), dir);
        auto pc = b_cells_.GetIdx(MIdx(mb[0], midx[1], midx[2]));
        ff_to_cell_[idxface][1] = GetCenter(pc) - GetCenter(pf);
      }
      #endif
      
      #ifdef PERZ
      // adhoc for periodic in z
      if (midx[2] == mb[2] && dir == Dir::k) {
        auto pf = b_faces_.GetIdx(MIdx(midx[0], midx[1], me[2]), dir);
        auto pc = b_cells_.GetIdx(MIdx(midx[0], midx[1], me[2] - 1));
        ff_to_cell_[idxface][0] = GetCenter(pc) - GetCenter(pf);
      }
      if (midx[2] == me[2] && dir == Dir::k) {
        auto pf = b_faces_.GetIdx(MIdx(midx[0], midx[1], mb[2]), dir);
        auto pc = b_cells_.GetIdx(MIdx(midx[0], midx[1], mb[2]));
        ff_to_cell_[idxface][1] = GetCenter(pc) - GetCenter(pf);
      }
      #endif
    }
  }

  // Mark inner faces
  for (auto idxface : this->Faces()) {
    ff_is_inner_[idxface] = true; 
  }
}

} // namespace geom3d
