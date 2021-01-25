// Created by Petr Karnakov on 24.11.2019
// Copyright 2019 ETH Zurich

#pragma once

#include <mpi.h>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <vector>

#include "distr.h"
#include "dump/output.h"
#include "dump/output_paraview.h"
#include "local.h"

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
  using RedOp = typename M::Op;
  using BlockInfoProxy = generic::BlockInfoProxy<3>;

  using P::blocksize_;
  using P::comm_;
  using P::dim;
  using P::extent_;
  using P::frame_;
  using P::halos_;
  using P::isroot_;
  using P::kernelfactory_;
  using P::kernels_;
  using P::nblocks_;
  using P::nprocs_;
  using P::stage_;
  using P::var;

  std::vector<FieldCell<Scal>> buf_; // buffer on mesh
  M globalmesh; // global mesh
  std::unique_ptr<output::Ser> oser_; // output series
  std::vector<BlockInfoProxy> proxies_;
  generic::Vect<bool, 3> per_; // periodic in direction

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
};

template <class M>
auto Local<M>::CreateGlobalMesh(MIdx bs, MIdx b, MIdx p, Scal ext) -> M {
  MIdx gs = bs * b * p; // global size in cells (without halos)
  Scal h = ext / gs.max(); // cell size
  Rect<Vect> d(Vect(0), Vect(gs) * h); // bounding box
  MIdx o(0); // origin index

  return InitUniformMesh<M>(d, o, gs, 0, true, true, gs, 0);
}

template <class M>
Local<M>::Local(MPI_Comm comm, const KernelMeshFactory<M>& kf, Vars& var_)
    : DistrMesh<M>(comm, kf, var_)
    , buf_(var.Int["loc_maxcomm"])
    , globalmesh(CreateGlobalMesh(blocksize_, nblocks_, nprocs_, extent_)) {
  // Resize buffer for mesh
  for (auto& u : buf_) {
    u.Reinit(globalmesh);
  }

  GBlock<size_t, dim> procs(nprocs_);
  GBlock<size_t, dim> blocks(nblocks_);
  for (auto wproc : procs) { // same ordering as in Cubism
    for (auto wblock : blocks) {
      BlockInfoProxy p;
      p.index = nblocks_ * wproc + wblock;
      p.cellsize = globalmesh.GetCellSize();
      p.blocksize = blocksize_;
      p.halos = halos_;
      p.isroot = (p.index == MIdx(0));
      p.islead = (p.index == MIdx(0));
      p.globalsize = nprocs_ * nblocks_ * blocksize_;
      proxies_.push_back(p);
    }
  }

  isroot_ = true; // XXX: overwrite isroot_

  // periodic
  per_[0] = var.Int["loc_periodic_x"];
  per_[1] = var.Int["loc_periodic_y"];
  per_[2] = var.Int["loc_periodic_z"];

  this->MakeKernels(proxies_);
}

template <class M>
auto Local<M>::TransferHalos(bool inner) -> std::vector<size_t> {
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
      if (!oser_) {
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
              on.second, globalmesh,
              [this, k](IdxCell i) { return buf_[k][i]; }));
          k += on.first->GetSize();
          if (on.first->GetSize() != 1) {
            throw std::runtime_error("DumpWrite(): Support only size 1");
          }
        }

        oser_.reset(new output::SerVtkStruct<M>(v, "title", "p", globalmesh));
      }

      // TODO: Check no change in comm list between time steps
      //       (otherwise oser_ needs reinitialization)

      oser_->Write(frame_, "title"); // TODO: t instead of stage_
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
  if (e >= buf_.size()) {
    throw std::runtime_error("ReadBuffer: Too many fields for Comm()");
  }
  auto& ndc = m.GetIndexCells();
  auto& gndc = globalmesh.GetIndexCells();
  MIdx gs = globalmesh.GetInBlockCells().GetSize();
  for (auto c : m.AllCells()) {
    auto w = ndc.GetMIdx(c);
    // periodic
    for (int d = 0; d < 3; ++d) {
      if (per_[d]) {
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
  if (e >= buf_.size()) {
    throw std::runtime_error("ReadBuffer: Too many fields for Comm()");
  }
  auto& indexc = m.GetIndexCells();
  auto& indexc_global = globalmesh.GetIndexCells();
  MIdx gs = globalmesh.GetInBlockCells().GetSize();
  for (auto c : m.AllCells()) {
    auto w = indexc.GetMIdx(c);
    // periodic
    for (int d = 0; d < 3; ++d) {
      if (per_[d]) {
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
  for (size_t d = 0; d < Vect::dim; ++d) {
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
  throw std::runtime_error("ReadBuffer: Unknown CommRequest instance");
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
  if (e >= buf_.size()) {
    throw std::runtime_error("WriteBuffer: Too many fields for Comm()");
  }
  auto& indexc = m.GetIndexCells();
  auto& indexc_global = globalmesh.GetIndexCells();
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
  if (e >= buf_.size()) {
    throw std::runtime_error("WriteBuffer: Too many fields for Comm()");
  }
  auto& indexc = m.GetIndexCells();
  auto& indexc_global = globalmesh.GetIndexCells();
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
  for (size_t d = 0; d < Vect::dim; ++d) {
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
  throw std::runtime_error("WriteBuffer: Unknown CommRequest instance");
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
