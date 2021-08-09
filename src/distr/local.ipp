// Created by Petr Karnakov on 24.11.2019
// Copyright 2019 ETH Zurich

#pragma once

#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <vector>

#include "distr.h"
#include "dump/output.h"
#include "dump/output_paraview.h"
#include "util/mpi.h"

template <class M_>
class Local : public DistrMesh<M_> {
 public:
  using M = M_;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;

  Local(MPI_Comm comm, const KernelMeshFactory<M>& kf, Vars& var);

 private:
  using MIdx = typename M::MIdx;
  using P = DistrMesh<M>;
  using P::dim;
  using RedOp = typename M::Op;
  using BlockInfoProxy = generic::BlockInfoProxy<dim>;

  using P::comm_;
  using P::domain_;
  using P::frame_;
  using P::isroot_;
  using P::kernelfactory_;
  using P::kernels_;
  using P::mshared_;
  using P::stage_;
  using P::var;

  std::vector<FieldCell<Scal>> buf_; // buffer on mesh
  M globalmesh_; // global mesh
  std::unique_ptr<output::Ser> output_; // output series
  std::vector<BlockInfoProxy> proxies_;
  generic::Vect<bool, dim> is_periodic_; // periodic in direction

  size_t WriteBuffer(const FieldCell<Scal>& fc, size_t e, M& m);
  size_t WriteBuffer(const FieldCell<Vect>& f, size_t d, size_t e, M& m);
  size_t WriteBuffer(const FieldCell<Vect>& fc, size_t e, M& m);
  size_t WriteBuffer(typename M::CommRequest* o, size_t e, M& m);
  void WriteBuffer(M& m);

  size_t ReadBuffer(FieldCell<Scal>& fc, size_t e, M& m);
  size_t ReadBuffer(FieldCell<Vect>& fc, size_t d, size_t e, M& m);
  size_t ReadBuffer(FieldCell<Vect>& fc, size_t e, M& m);
  size_t ReadBuffer(typename M::CommRequest* o, size_t e, M& m);
  void ReadBuffer(M& m);

  // bs: inner block size
  // b: blocks per rank
  // p: ranks
  // ext: extent
  static M CreateGlobalMesh(MIdx bs, MIdx b, MIdx p, Scal ext);

  std::vector<size_t> TransferHalos(bool inner) override;
  void ReduceSingleRequest(const std::vector<RedOp*>& blocks) override;
  void Scatter(const std::vector<size_t>& bb) override;
  void Bcast(const std::vector<size_t>& bb) override;
  void DumpWrite(const std::vector<size_t>& bb) override;
  void TransferParticles(const std::vector<size_t>& bb) override;
};

template <class M>
auto Local<M>::CreateGlobalMesh(MIdx bs, MIdx b, MIdx p, Scal ext) -> M {
  const MIdx gs = bs * b * p; // global size in cells (without halos)
  const Scal h = ext / gs.max(); // cell size
  const Rect<Vect> domain(Vect(0), Vect(gs) * h); // bounding box
  const MIdx begin(0);

  return {begin, gs, domain, 0, true, true, gs, 0};
}

template <class M>
Local<M>::Local(MPI_Comm comm, const KernelMeshFactory<M>& kf, Vars& var_)
    : DistrMesh<M>(comm, kf, var_)
    , buf_(var.Int["loc_maxcomm"])
    , globalmesh_(CreateGlobalMesh(
          domain_.blocksize, domain_.nblocks, domain_.nprocs, domain_.extent)) {
  // Resize buffer for mesh
  for (auto& u : buf_) {
    u.Reinit(globalmesh_);
  }

  GBlock<size_t, dim> procs(domain_.nprocs);
  GBlock<size_t, dim> blocks(domain_.nblocks);
  for (auto wproc : procs) { // same ordering as in Cubism
    for (auto wblock : blocks) {
      BlockInfoProxy p;
      p.index = domain_.nblocks * wproc + wblock;
      p.cellsize = globalmesh_.GetCellSize();
      p.blocksize = domain_.blocksize;
      p.halos = domain_.halos;
      p.isroot = (p.index == MIdx(0));
      p.islead = (p.index == MIdx(0));
      p.globalsize = domain_.nprocs * domain_.nblocks * domain_.blocksize;
      proxies_.push_back(p);
    }
  }

  isroot_ = true; // XXX: overwrite isroot_
  is_periodic_[0] = var.Int["loc_periodic_x"];
  is_periodic_[1] = var.Int["loc_periodic_y"];
  is_periodic_[2] = var.Int["loc_periodic_z"];

  this->MakeKernels(proxies_);

  for (auto& kernel : kernels_) {
    auto& m = kernel->GetMesh();
    m.SetHandlerMpiRankFromId([](int) -> int { //
      return 0;
    });
  }
}

template <class M>
auto Local<M>::TransferHalos(bool inner) -> std::vector<size_t> {
  fassert_equal(
      mshared_->GetComm().size(), 0,
      ". Comm() on shared mesh is not implemented");

  std::vector<size_t> bb;
  if (inner) {
    bb.resize(proxies_.size());
    std::iota(bb.begin(), bb.end(), 0);
    for (auto b : bb) {
      WriteBuffer(kernels_[b]->GetMesh());
    }
    for (auto b : bb) {
      ReadBuffer(kernels_[b]->GetMesh());
    }
  }
  return bb;
}

template <class M>
void Local<M>::ReduceSingleRequest(const std::vector<RedOp*>& blocks) {
  using OpScal = typename UReduce<Scal>::OpS;
  using OpScalInt = typename UReduce<Scal>::OpSI;
  using OpConcat = typename UReduce<Scal>::OpCat;

  auto* firstbase = blocks.front();

  if (auto* first = dynamic_cast<OpScal*>(firstbase)) {
    auto buf = first->Neutral();

    for (auto otherbase : blocks) {
      auto* other = dynamic_cast<OpScal*>(otherbase);
      other->Append(buf);
    }

    for (auto otherbase : blocks) {
      auto* other = dynamic_cast<OpScal*>(otherbase);
      other->Set(buf);
    }
    return;
  }
  if (auto* first = dynamic_cast<OpScalInt*>(firstbase)) {
    auto buf = first->Neutral();

    for (auto otherbase : blocks) {
      auto* other = dynamic_cast<OpScalInt*>(otherbase);
      other->Append(buf);
    }

    for (auto otherbase : blocks) {
      auto* other = dynamic_cast<OpScalInt*>(otherbase);
      other->Set(buf);
    }
    return;
  }
  if (auto* first = dynamic_cast<OpConcat*>(firstbase)) {
    // Concatenation of std::vector<T>
    auto buf = first->Neutral();

    for (auto otherbase : blocks) {
      auto* other = dynamic_cast<OpConcat*>(otherbase);
      other->Append(buf);
    }

    // Write results to root block only (assume first)
    // FIXME get IsRoot() from kernel_ after using std::vector for kernels
    first->Set(buf);
    return;
  }
  fassert(false, "Reduce: Unknown M::Op implementation");
}

template <class M>
void Local<M>::Scatter(const std::vector<size_t>& bb) {
  const size_t nreq = kernels_.front()->GetMesh().GetScatter().size();

  for (auto b : bb) {
    fassert_equal(kernels_[b]->GetMesh().GetScatter().size(), nreq);
  }

  for (size_t q = 0; q < nreq; ++q) {
    for (auto broot : bb) {
      auto& mroot = kernels_[broot]->GetMesh();
      if (mroot.IsRoot()) {
        auto& req = mroot.GetScatter()[q];
        size_t i = 0;
        for (auto b : bb) {
          auto& v = *kernels_[b]->GetMesh().GetScatter()[q].second;
          v = (*req.first)[i++];
        }
      }
    }
  }

  for (auto b : bb) {
    kernels_[b]->GetMesh().ClearScatter();
  }
}

template <class M>
void Local<M>::Bcast(const std::vector<size_t>& bb) {
  using OpConcat = typename UReduce<Scal>::OpCat;
  auto& vf = kernels_.front()->GetMesh().GetBcast(); // pointers to broadcast

  for (auto b : bb) {
    fassert_equal(kernels_[b]->GetMesh().GetBcast().size(), vf.size());
  }

  for (size_t i = 0; i < vf.size(); ++i) {
    if (OpConcat* first = dynamic_cast<OpConcat*>(vf[i].get())) {
      std::vector<char> buf = first->Neutral();

      // read from root block
      for (auto b : bb) {
        auto& m = kernels_[b]->GetMesh();
        if (m.IsRoot()) {
          auto& v = m.GetBcast();
          OpConcat* other = dynamic_cast<OpConcat*>(v[i].get());
          other->Append(buf);
        }
      }

      // write to all blocks
      for (auto b : bb) {
        auto& m = kernels_[b]->GetMesh();
        auto& v = m.GetBcast();
        OpConcat* other = dynamic_cast<OpConcat*>(v[i].get());
        other->Set(buf);
      }
    } else {
      fassert(false, "Bcast: Unknown M::Op instance");
    }
  }

  for (auto b : bb) {
    kernels_[b]->GetMesh().ClearBcast();
  }
}

template <class M>
void Local<M>::DumpWrite(const std::vector<size_t>& bb) {
  auto& m = kernels_.front()->GetMesh();
  if (m.GetDump().size()) {
    std::string dumpformat = var.String["dumpformat"];
    if (dumpformat == "default") {
      dumpformat = "vtk";
    }

    if (dumpformat == "vtk") {
      // Initialize on first call
      if (!output_) {
        // TODO: check all blocks are same as first
        output::VOut v;
        size_t k = 0; // offset in buffer
        // Skip comm
        for (auto& o : m.GetComm()) {
          k += o->GetSize();
        }
        // Write dump
        for (auto& on : m.GetDump()) {
          v.emplace_back(new output::OutFldFunc<Scal, IdxCell, M>(
              on.second, globalmesh_,
              [this, k](IdxCell i) { return buf_[k][i]; }));
          k += on.first->GetSize();
          fassert_equal(
              on.first->GetSize(), 1, ". DumpWrite(): Support only size 1");
        }

        output_.reset(
            new output::SerVtkStruct<M>(v, "title", "p", globalmesh_));
      }

      // TODO: Check no change in comm list between time steps
      //       (otherwise output_ needs reinitialization)

      output_->Write(frame_, "title"); // TODO: t instead of stage_
      std::cerr << "Dump " << frame_ << ": format=" << dumpformat << std::endl;
      ++frame_;
    } else {
      P::DumpWrite(bb);
    }
  }
}

// Reads from buffer to scalar field [a]
// fc: field
// l: lab
// e: offset in buffer
// Returns:
// number of scalar fields read
template <class M>
size_t Local<M>::ReadBuffer(FieldCell<Scal>& fc, size_t e, M& m) {
  fassert(e < buf_.size(), "ReadBuffer: Too many fields for Comm()");
  auto& ndc = m.GetIndexCells();
  auto& gndc = globalmesh_.GetIndexCells();
  MIdx gs = globalmesh_.GetInBlockCells().GetSize();
  for (auto c : m.AllCells()) {
    auto w = ndc.GetMIdx(c);
    // periodic
    for (auto d : M::dirs) {
      if (is_periodic_[d]) {
        w[d] = (w[d] + gs[d]) % gs[d];
      }
    }
    // XXX: addhoc nan in halos
    if (MIdx(0) <= w && w < gs) {
      auto gc = gndc.GetIdx(w);
      fc[c] = buf_[e][gc];
    } else {
      fc[c] = std::numeric_limits<Scal>::quiet_NaN();
    }
  }
  return 1;
}

// Reads from buffer to component of vector field [a]
// fc: field
// comp: component (0,1,2)
// e: offset in buffer
// Returns:
// number of scalar fields read
template <class M>
size_t Local<M>::ReadBuffer(FieldCell<Vect>& fc, size_t comp, size_t e, M& m) {
  fassert(e < buf_.size(), "ReadBuffer: Too many fields for Comm()");
  auto& indexc = m.GetIndexCells();
  auto& indexc_global = globalmesh_.GetIndexCells();
  MIdx gs = globalmesh_.GetInBlockCells().GetSize();
  for (auto c : m.AllCells()) {
    auto w = indexc.GetMIdx(c);
    // periodic
    for (auto d : M::dirs) {
      if (is_periodic_[d]) {
        w[d] = (w[d] + gs[d]) % gs[d];
      }
    }
    // XXX: addhoc nan in halos
    if (MIdx(0) <= w && w < gs) {
      auto gc = indexc_global.GetIdx(w);
      fc[c][comp] = buf_[e][gc];
    } else {
      fc[c][comp] = std::numeric_limits<Scal>::quiet_NaN();
    }
  }
  return 1;
}

// Reads from buffer to all components of vector field [a]
// fc: field
// e: offset in buffer
// Returns:
// number of scalar fields read
template <class M>
size_t Local<M>::ReadBuffer(FieldCell<Vect>& fc, size_t e, M& m) {
  for (auto d : M::dirs) {
    e += ReadBuffer(fc, d, e, m);
  }
  return Vect::dim;
}
// Reads from buffer to CommRequest
// o: instance of CommRequest
// e: offset in buffer
// Returns:
// number of scalar fields written
template <class M>
size_t Local<M>::ReadBuffer(typename M::CommRequest* o, size_t e, M& m) {
  if (auto od = dynamic_cast<typename M::CommRequestScal*>(o)) {
    return ReadBuffer(*od->field, e, m);
  }
  if (auto od = dynamic_cast<typename M::CommRequestVect*>(o)) {
    if (od->d == -1) {
      return ReadBuffer(*od->field, e, m);
    }
    return ReadBuffer(*od->field, od->d, e, m);
  }
  fassert(false, "ReadBuffer: Unknown CommRequest instance");
  return 0;
}
template <class M>
void Local<M>::ReadBuffer(M& m) {
  size_t e = 0;
  for (auto& o : m.GetComm()) {
    e += ReadBuffer(o.get(), e, m);
  }
}

// Writes scalar field to buffer [i].
// fc: scalar field
// e: offset in buffer
// Returns:
// number of scalar fields written
template <class M>
size_t Local<M>::WriteBuffer(const FieldCell<Scal>& fc, size_t e, M& m) {
  fassert(e < buf_.size(), "WriteBuffer: Too many fields for Comm()");
  auto& indexc = m.GetIndexCells();
  auto& indexc_global = globalmesh_.GetIndexCells();
  for (auto c : m.Cells()) {
    auto w = indexc.GetMIdx(c);
    auto gc = indexc_global.GetIdx(w);
    buf_[e][gc] = fc[c];
  }
  return 1;
}
// Writes component of vector field to buffer [i].
// fc: vector field
// d: component (0,1,2)
// e: offset in buffer
// Returns:
// number of scalar fields written
template <class M>
size_t Local<M>::WriteBuffer(
    const FieldCell<Vect>& fc, size_t d, size_t e, M& m) {
  fassert(e < buf_.size(), "WriteBuffer: Too many fields for Comm()");
  auto& indexc = m.GetIndexCells();
  auto& indexc_global = globalmesh_.GetIndexCells();
  for (auto c : m.Cells()) {
    auto w = indexc.GetMIdx(c);
    auto gc = indexc_global.GetIdx(w);
    buf_[e][gc] = fc[c][d];
  }
  return 1;
}
// Writes all components of vector field to buffer [i].
// fc: vector field
// e: offset in buffer
// Returns:
// number of scalar fields written
template <class M>
size_t Local<M>::WriteBuffer(const FieldCell<Vect>& fc, size_t e, M& m) {
  for (auto d : M::dirs) {
    e += WriteBuffer(fc, d, e, m);
  }
  return Vect::dim;
}
// Writes CommRequest to buffer.
// o: instance of CommRequest
// e: offset in buffer
// Returns:
// number of scalar fields written
template <class M>
size_t Local<M>::WriteBuffer(typename M::CommRequest* o, size_t e, M& m) {
  if (auto od = dynamic_cast<typename M::CommRequestScal*>(o)) {
    return WriteBuffer(*od->field, e, m);
  }
  if (auto od = dynamic_cast<typename M::CommRequestVect*>(o)) {
    if (od->d == -1) {
      return WriteBuffer(*od->field, e, m);
    }
    return WriteBuffer(*od->field, od->d, e, m);
  }
  fassert(false, "WriteBuffer: Unknown CommRequest instance");
}
template <class M>
void Local<M>::WriteBuffer(M& m) {
  size_t e = 0;
  for (auto& o : m.GetComm()) {
    e += WriteBuffer(o.get(), e, m);
  }
  for (auto& on : m.GetDump()) {
    e += WriteBuffer(on.first.get(), e, m);
  }
}

template <class M>
void Local<M>::TransferParticles(const std::vector<size_t>& bb) {
  const size_t nreq = kernels_.front()->GetMesh().GetCommPart().size();

  for (auto b : bb) {
    fassert_equal(kernels_[b]->GetMesh().GetCommPart().size(), nreq);
  }

  for (size_t q = 0; q < nreq; ++q) {
    auto& mroot = kernels_.front()->GetMesh();
    const auto halorad =
        mroot.flags.particles_halo_radius * mroot.GetCellSize();
    for (size_t b : bb) {
      auto& m = kernels_[b]->GetMesh();
      auto& req = m.GetCommPart()[q];
      const auto& bbox = m.GetBoundingBox();
      // TODO: consider periodic
      const Rect<Vect> halobox(bbox.low - halorad, bbox.high + halorad);
      std::vector<Vect> xx;
      std::vector<std::vector<Scal>> attr_scal(req.attr_scal.size());
      std::vector<std::vector<Vect>> attr_vect(req.attr_vect.size());
      for (size_t bn : bb) {
        auto& reqn = kernels_[bn]->GetMesh().GetCommPart()[q];
        for (size_t i = 0; i < reqn.x->size(); ++i) {
          const auto& x = (*reqn.x)[i];
          if (halobox.low <= x && x < halobox.high) {
            xx.push_back(x);
            for (size_t a = 0; a < attr_scal.size(); ++a) {
              attr_scal[a].push_back((*reqn.attr_scal[a])[i]);
            }
            for (size_t a = 0; a < attr_vect.size(); ++a) {
              attr_vect[a].push_back((*reqn.attr_vect[a])[i]);
            }
          }
        }
      }
    }
  }

  for (auto b : bb) {
    kernels_[b]->GetMesh().ClearCommPart();
  }
}
