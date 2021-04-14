// Created by Petr Karnakov on 27.10.2020
// Copyright 2020 ETH Zurich

#include "mesh.h"

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
  std::vector<std::unique_ptr<typename UReduce<Scal>::Op>> bcast_lead;
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
    , blockcs_(b - MIdx(1), cs + MIdx(2))
    , blockfs_(blockcs_.GetBegin(), blockcs_.GetSize())
    , blockns_(blockcs_.GetBegin(), blockcs_.GetSize() + MIdx(1))
    // index
    , indexc_(blockca_.GetBegin(), blockca_.GetSize() + MIdx(1))
    , indexf_(indexc_.GetBegin(), indexc_.GetSize())
    , indexn_(indexc_.GetBegin(), indexc_.GetSize())
    , mshared_(this)
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
    using HyperMIdx = generic::MIdx<dim - 1>;
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

  { // face centers
    ff_center_.Reinit(*this);
    for (auto f : AllFaces()) {
      auto p = indexf_.GetMIdxDir(f);
      const MIdx& w = p.first;
      size_t d(p.second);
      Vect r = (domain_.low + half_cell_size_) +
               Vect(w - incells_begin_) * cell_size_;
      r[d] -= half_cell_size_[d];
      ff_center_[f] = r;
    }
  }
}

template <class Scal, size_t dim>
MeshCartesian<Scal, dim>::~MeshCartesian() = default;

template <class Scal, size_t dim>
MeshCartesian<Scal, dim>::MeshCartesian(MeshCartesian&&) = default;

template <class Scal, size_t dim>
void MeshCartesian<Scal, dim>::Comm(std::unique_ptr<CommRequest>&& r) {
  imp->commreq.emplace_back(std::move(r));
}
template <class Scal, size_t dim>
void MeshCartesian<Scal, dim>::Comm(FieldCell<Scal>* f, CommStencil stencil) {
  Comm(std::make_unique<CommRequestScal>(f, stencil));
}
template <class Scal, size_t dim>
void MeshCartesian<Scal, dim>::Comm(
    FieldCell<Vect>* f, int d, CommStencil stencil) {
  Comm(std::make_unique<CommRequestVect>(f, d, stencil));
}
template <class Scal, size_t dim>
void MeshCartesian<Scal, dim>::Comm(FieldCell<Vect>* f, CommStencil stencil) {
  Comm(f, -1, stencil);
}
template <class Scal, size_t dim>
void MeshCartesian<Scal, dim>::Dump(const FieldCell<Scal>* f, std::string n) {
  auto ff = const_cast<FieldCell<Scal>*>(f);
  imp->dump.emplace_back(
      std::make_unique<CommRequestScal>(ff, CommStencil::none), n);
}
template <class Scal, size_t dim>
void MeshCartesian<Scal, dim>::Dump(
    const FieldCell<Vect>* f, int d, std::string n) {
  auto ff = const_cast<FieldCell<Vect>*>(f);
  imp->dump.emplace_back(
      std::make_unique<CommRequestVect>(ff, d, CommStencil::none), n);
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
void MeshCartesian<Scal, dim>::BcastFromLead(
    std::unique_ptr<typename UReduce<Scal>::Op>&& o) {
  imp->bcast_lead.emplace_back(std::move(o));
}
template <class Scal, size_t dim>
auto MeshCartesian<Scal, dim>::GetBcastFromLead() const
    -> const std::vector<std::unique_ptr<Op>>& {
  return imp->bcast_lead;
}
template <class Scal, size_t dim>
void MeshCartesian<Scal, dim>::ClearBcastFromLead() {
  imp->bcast_lead.clear();
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
