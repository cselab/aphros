// Created by Petr Karnakov on 27.10.2020
// Copyright 2020 ETH Zurich

#include "mesh.h"

using M = MeshStructured<double, 3>;

template <class _Scal, size_t _dim>
struct MeshStructured<_Scal, _dim>::Imp {
  using Owner = MeshStructured<_Scal, _dim>;
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
MeshStructured<_Scal, _dim>::MeshStructured(
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
    , imp(new Imp(this)) {
  static_assert(dim == 3, "Not implemented for dim != 3");

  // surface area
  face_area_[0] = cell_size_[1] * cell_size_[2];
  face_area_[1] = cell_size_[2] * cell_size_[0];
  face_area_[2] = cell_size_[0] * cell_size_[1];

  // surface vectors
  face_surface_[0] = Vect(face_area_[0], 0., 0.);
  face_surface_[1] = Vect(0., face_area_[1], 0.);
  face_surface_[2] = Vect(0., 0., face_area_[2]);

  { // cell neighbour cell offset
    MIdx w = indexc_.GetBegin(); // any cell
    IdxCell c = indexc_.GetIdx(w);
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
      IdxCell cn = indexc_.GetIdx(w + wo);
      cell_cell_[q] = size_t(cn) - size_t(c);
    }
  }

  { // cell neighbour face offset
    MIdx w = indexc_.GetBegin(); // any cell
    IdxCell c = indexc_.GetIdx(w);
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
      IdxFace fn = indexf_.GetIdx(std::make_pair(w + wo, d));
      cell_face_[q] = size_t(fn) - size_t(c);
    }
  }

  { // cell outward factor
    for (size_t q = 0; q < kCellNumNeighbourFaces; ++q) {
      cell_outward_[q] = (q % 2 == 0 ? -1. : 1.);
    }
  }

  { // cell neighbour node offset
    MIdx w = indexc_.GetBegin(); // any cell
    IdxCell c = indexc_.GetIdx(w);
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
      IdxNode nn = indexn_.GetIdx(w + wo);
      cell_node_[q] = size_t(nn) - size_t(c);
    }
  }

  { // face neighbour cell offset
    const MIdx w = indexc_.GetBegin(); // any cell
    for (size_t d = 0; d < dim; ++d) {
      IdxFace f = indexf_.GetIdx(w, Dir(d));
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
        IdxCell cn = indexc_.GetIdx(w + wo);
        face_cell_[d * qm + q] = size_t(cn) - size_t(f);
      }
    }
  }

  { // face neighbour face offset
    const MIdx w = indexc_.GetBegin(); // any cell
    for (size_t d = 0; d < dim; ++d) {
      const IdxFace f = indexf_.GetIdx(w, Dir(d));
      auto qm = kCellNumNeighbourFaces;
      for (size_t q = 0; q < qm; ++q) {
        MIdx wo(0);
        wo[q / 2] = (q % 2 == 0 ? -1 : 1);
        const IdxFace fn = indexf_.GetIdx(w + wo, Dir(d));
        face_face_[d * qm + q] = size_t(fn) - size_t(f);
      }
    }
  }

  { // face neighbour cell offset
    const MIdx w = indexc_.GetBegin(); // any cell
    for (size_t d = 0; d < dim; ++d) {
      IdxFace f = indexf_.GetIdx(w, Dir(d));
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
        IdxNode nn = indexn_.GetIdx(w + wo);
        face_node_[d * qm + q] = size_t(nn) - size_t(f);
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
MeshStructured<Scal, dim>::~MeshStructured() = default;

template <class Scal, size_t dim>
MeshStructured<Scal, dim>::MeshStructured(MeshStructured&&) = default;

template <class M>
M InitUniformMesh(
    Rect<typename M::Vect> domain, typename M::MIdx begin, typename M::MIdx s,
    int halos, bool isroot, bool islead, typename M::MIdx gs, int id) {
  return {begin, s, domain, halos, isroot, islead, gs, id};
}

template <class Scal, size_t dim>
void MeshStructured<Scal, dim>::Comm(std::unique_ptr<CommRequest>&& r) {
  imp->commreq.emplace_back(std::move(r));
}
template <class Scal, size_t dim>
void MeshStructured<Scal, dim>::Comm(FieldCell<Scal>* f) {
  Comm(std::make_unique<CommRequestScal>(f));
}
template <class Scal, size_t dim>
void MeshStructured<Scal, dim>::Comm(FieldCell<Vect>* f, int d) {
  Comm(std::make_unique<CommRequestVect>(f, d));
}
template <class Scal, size_t dim>
void MeshStructured<Scal, dim>::Comm(FieldCell<Vect>* f) {
  Comm(f, -1);
}
template <class Scal, size_t dim>
void MeshStructured<Scal, dim>::Dump(const FieldCell<Scal>* f, std::string n) {
  auto ff = const_cast<FieldCell<Scal>*>(f);
  imp->dump.emplace_back(std::make_unique<CommRequestScal>(ff), n);
}
template <class Scal, size_t dim>
void MeshStructured<Scal, dim>::Dump(
    const FieldCell<Vect>* f, int d, std::string n) {
  auto ff = const_cast<FieldCell<Vect>*>(f);
  imp->dump.emplace_back(std::make_unique<CommRequestVect>(ff, d), n);
}
template <class Scal, size_t dim>
void MeshStructured<Scal, dim>::Dump(
    std::unique_ptr<CommRequest>&& o, std::string name) {
  imp->dump.emplace_back(std::move(o), name);
}
template <class Scal, size_t dim>
auto MeshStructured<Scal, dim>::GetComm() const
    -> const std::vector<std::unique_ptr<CommRequest>>& {
  return imp->commreq;
}
template <class Scal, size_t dim>
auto MeshStructured<Scal, dim>::GetDump() const -> const
    std::vector<std::pair<std::unique_ptr<CommRequest>, std::string>>& {
  return imp->dump;
}
template <class Scal, size_t dim>
void MeshStructured<Scal, dim>::ClearComm() {
  imp->commreq.clear();
}
template <class Scal, size_t dim>
void MeshStructured<Scal, dim>::ClearDump() {
  imp->dump.clear();
}
template <class Scal, size_t dim>
void MeshStructured<Scal, dim>::Reduce(
    std::unique_ptr<typename UReduce<Scal>::Op>&& o) {
  imp->reduce.Add(std::move(o));
}
template <class Scal, size_t dim>
void MeshStructured<Scal, dim>::Reduce(Scal* u, std::string o) {
  imp->reduce.Add(u, o);
}
template <class Scal, size_t dim>
const std::vector<std::unique_ptr<typename UReduce<Scal>::Op>>&
MeshStructured<Scal, dim>::GetReduce() const {
  return imp->reduce.Get();
}
template <class Scal, size_t dim>
void MeshStructured<Scal, dim>::Reduce(Scal* buf, ReductionType::Sum) {
  Reduce(buf, "sum");
}
template <class Scal, size_t dim>
void MeshStructured<Scal, dim>::Reduce(Scal* buf, ReductionType::Prod) {
  Reduce(buf, "prod");
}
template <class Scal, size_t dim>
void MeshStructured<Scal, dim>::Reduce(Scal* buf, ReductionType::Max) {
  Reduce(buf, "max");
}
template <class Scal, size_t dim>
void MeshStructured<Scal, dim>::Reduce(Scal* buf, ReductionType::Min) {
  Reduce(buf, "min");
}
template <class Scal, size_t dim>
void MeshStructured<Scal, dim>::Reduce(
    std::pair<Scal, int>* buf, ReductionType::MaxLoc) {
  Reduce(std::make_unique<typename UReduce<Scal>::OpMaxloc>(buf));
}
template <class Scal, size_t dim>
void MeshStructured<Scal, dim>::Reduce(
    std::pair<Scal, int>* buf, ReductionType::MinLoc) {
  Reduce(std::make_unique<typename UReduce<Scal>::OpMinloc>(buf));
}
template <class Scal, size_t dim>
void MeshStructured<Scal, dim>::ClearReduce() {
  imp->reduce.Clear();
}
template <class Scal, size_t dim>
void MeshStructured<Scal, dim>::ReduceToLead(
    std::unique_ptr<typename UReduce<Scal>::Op>&& o) {
  imp->reduce_lead.Add(std::move(o));
}
template <class Scal, size_t dim>
const std::vector<std::unique_ptr<typename UReduce<Scal>::Op>>&
MeshStructured<Scal, dim>::GetReduceToLead() const {
  return imp->reduce_lead.Get();
}
template <class Scal, size_t dim>
void MeshStructured<Scal, dim>::ClearReduceToLead() {
  imp->reduce_lead.Clear();
}

template <class Scal, size_t dim>
void MeshStructured<Scal, dim>::Scatter(const ScatterRequest& req) {
  imp->scatter.push_back(req);
}
template <class Scal, size_t dim>
auto MeshStructured<Scal, dim>::GetScatter() const
    -> const std::vector<ScatterRequest>& {
  return imp->scatter;
}
template <class Scal, size_t dim>
void MeshStructured<Scal, dim>::ClearScatter() {
  imp->scatter.clear();
}
template <class Scal, size_t dim>
void MeshStructured<Scal, dim>::Bcast(
    std::unique_ptr<typename UReduce<Scal>::Op>&& o) {
  imp->bcast.emplace_back(std::move(o));
}
template <class Scal, size_t dim>
auto MeshStructured<Scal, dim>::GetBcast() const
    -> const std::vector<std::unique_ptr<Op>>& {
  return imp->bcast;
}
template <class Scal, size_t dim>
void MeshStructured<Scal, dim>::ClearBcast() {
  imp->bcast.clear();
}
template <class Scal, size_t dim>
auto MeshStructured<Scal, dim>::GetNanFaces() const
    -> const std::vector<std::pair<IdxFace, size_t>>& {
  return imp->vfnan;
}
template <class Scal, size_t dim>
void MeshStructured<Scal, dim>::SetNanFaces(
    const std::vector<std::pair<IdxFace, size_t>>& vfnan) {
  imp->vfnan = vfnan;
}
template <class Scal, size_t dim>
template <class T>
void MeshStructured<Scal, dim>::ApplyNanFaces(FieldCell<T>& fc) {
  for (auto p : imp->vfnan) {
    IdxFace f = p.first;
    size_t nci = p.second;
    auto cc = GetCellColumn(f, nci);
    fc[cc[0]] = T(flags.nan_faces_value);
    fc[cc[1]] = T(flags.nan_faces_value);
  }
}
