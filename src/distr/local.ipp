// Created by Petr Karnakov on 24.11.2019
// Copyright 2019 ETH Zurich

#pragma once

#include <mpi.h>
#include <limits>
#include <map>
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

  using P::b_;
  using P::bs_;
  using P::comm_;
  using P::dim;
  using P::ext_;
  using P::frame_;
  using P::hl_;
  using P::isroot_;
  using P::kernels_;
  using P::kf_;
  using P::p_;
  using P::stage_;
  using P::var;

  std::vector<FieldCell<Scal>> buf_; // buffer on mesh
  M gm; // global mesh
  std::unique_ptr<output::Ser> oser_; // output series
  std::vector<MyBlockInfo> bb_;
  generic::Vect<bool, 3> per_; // periodic in direction

  size_t WriteBuffer(const FieldCell<Scal>& fc, size_t e, M& m);
  size_t WriteBuffer(const FieldCell<Vect>& f, size_t d, size_t e, M& m);
  size_t WriteBuffer(const FieldCell<Vect>& fc, size_t e, M& m);
  size_t WriteBuffer(typename M::Co* o, size_t e, M& m);
  void WriteBuffer(M& m);

  size_t ReadBuffer(FieldCell<Scal>& fc, size_t e, M& m);
  size_t ReadBuffer(FieldCell<Vect>& fc, size_t d, size_t e, M& m);
  size_t ReadBuffer(FieldCell<Vect>& fc, size_t e, M& m);
  size_t ReadBuffer(typename M::Co* o, size_t e, M& m);
  void ReadBuffer(M& m);

  // bs: inner block size
  // b: blocks per rank
  // p: ranks
  // ext: extent
  static M CreateGlobalMesh(MIdx bs, MIdx b, MIdx p, Scal ext);

  std::vector<MIdx> GetBlocks(bool inner) override;
  void ReadBuffer(const std::vector<MIdx>& bb) override;
  void WriteBuffer(const std::vector<MIdx>& bb) override;
  void ReduceSingleRequest(
      const std::vector<std::shared_ptr<RedOp>>& blocks) override;
  void Scatter(const std::vector<MIdx>& bb) override;
  void Bcast(const std::vector<MIdx>& bb) override;
  void DumpWrite(const std::vector<MIdx>& bb) override;
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
    , gm(CreateGlobalMesh(bs_, b_, p_, ext_)) {
  // Resize buffer for mesh
  for (auto& u : buf_) {
    u.Reinit(gm);
  }

  // Fill block info
  MIdx ms(bs_); // block size
  MIdx mb(b_[0], b_[1], b_[2]); // number of blocks
  MIdx mp(p_[0], p_[1], p_[2]); // number of PEs
  GBlock<size_t, dim> procs(p_);
  GBlock<size_t, dim> blocks(b_);
  Scal h = (gm.GetNode(IdxNode(1)) - gm.GetNode(IdxNode(0)))[0];
  assert(h > 0);
  if (var.Int["verbose"]) {
    std::cerr << "h from gm = " << h << std::endl;
  }
  for (auto wp : procs) { // same ordering as in Cubism
    for (auto wb : blocks) {
      auto w = b_ * wp + wb;
      MyBlockInfo b;
      IdxNode n = gm.GetIndexNodes().GetIdx(w * ms);
      Vect o = gm.GetNode(n);
      if (var.Int["verbose"]) {
        std::cerr << "o=" << o << " n=" << n.GetRaw() << " w=" << w
                  << std::endl;
      }
      for (int q = 0; q < 3; ++q) {
        b.index[q] = w[q];
        b.origin[q] = o[q];
        b.bs[q] = bs_[q];
      }
      b.h_gridpoint = h;
      b.hl = hl_;
      b.maxcomm = buf_.size();
      MIdx gs = p_ * b_ * bs_; // global size
      for (int j = 0; j < 3; ++j) {
        b.gs[j] = gs[j];
      }
      bb_.push_back(b);
    }
  }

  isroot_ = true; // XXX: overwrite isroot_

  // periodic
  per_[0] = var.Int["loc_periodic_x"];
  per_[1] = var.Int["loc_periodic_y"];
  per_[2] = var.Int["loc_periodic_z"];

  bool islead = true;
  for (auto& b : bb_) {
    MIdx d(b.index);
    b.isroot = (d == MIdx(0));
    b.islead = islead;
    islead = false;
  }

  this->MakeKernels(bb_);
}

template <class M>
auto Local<M>::GetBlocks(bool inner) -> std::vector<MIdx> {
  std::vector<MIdx> bb;
  if (inner) {
    for (auto e : bb_) {
      bb.emplace_back(e.index);
    }
  }
  return bb;
}

template <class M>
void Local<M>::ReadBuffer(const std::vector<MIdx>& bb) {
  for (auto& b : bb) {
    auto& m = kernels_.at(b)->GetMesh();
    ReadBuffer(m);
  }
}

template <class M>
void Local<M>::WriteBuffer(const std::vector<MIdx>& bb) {
  for (auto& b : bb) {
    auto& m = kernels_.at(b)->GetMesh();
    WriteBuffer(m);
  }
}

template <class M>
void Local<M>::ReduceSingleRequest(
    const std::vector<std::shared_ptr<RedOp>>& blocks) {
  using OpScal = typename M::OpS;
  using OpScalInt = typename M::OpSI;
  using OpConcat = typename M::OpCat;

  auto* firstbase = blocks.front().get();

  if (auto* first = dynamic_cast<OpScal*>(firstbase)) {
    auto buf = first->Neutral();

    for (auto otherbase : blocks) {
      auto* other = dynamic_cast<OpScal*>(otherbase.get());
      other->Append(buf);
    }

    for (auto otherbase : blocks) {
      auto* other = dynamic_cast<OpScal*>(otherbase.get());
      other->Set(buf);
    }
    return;
  }
  if (auto* first = dynamic_cast<OpScalInt*>(firstbase)) {
    auto buf = first->Neutral();

    for (auto otherbase : blocks) {
      auto* other = dynamic_cast<OpScalInt*>(otherbase.get());
      other->Append(buf);
    }

    for (auto otherbase : blocks) {
      auto* other = dynamic_cast<OpScalInt*>(otherbase.get());
      other->Set(buf);
    }
    return;
  }
  if (auto* first = dynamic_cast<OpConcat*>(firstbase)) {
    // Concatenation of std::vector<T>
    auto buf = first->Neutral();

    for (auto otherbase : blocks) {
      auto* other = dynamic_cast<OpConcat*>(otherbase.get());
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
void Local<M>::Scatter(const std::vector<MIdx>& bb) {
  auto& vreq0 =
      kernels_.at(bb[0])->GetMesh().GetScatter(); // requests on first block

  // Check size is the same for all blocks
  for (auto& b : bb) {
    auto& vreq = kernels_.at(b)->GetMesh().GetScatter();
    if (vreq.size() != vreq0.size()) {
      throw std::runtime_error("Scatter: vreq.size() != vreq0.size()");
    }
  }

  for (size_t q = 0; q < vreq0.size(); ++q) {
    // find root block
    for (auto& b : bb) {
      auto& m = kernels_.at(b)->GetMesh();
      if (m.IsRoot()) {
        auto& req = m.GetScatter()[q];

        GBlock<size_t, dim> qp(p_);
        GBlock<size_t, dim> qb(b_);

        // write to blocks on current rank
        size_t i = 0;
        for (auto wp : qp) { // same ordering as with Cubism
          for (auto wb : qb) {
            auto w = b_ * wp + wb;
            auto& v = *kernels_.at(w)->GetMesh().GetScatter()[q].second;
            v = (*req.first)[i++];
          }
        }
      }
    }
  }

  // Clear requests
  for (auto& b : bb) {
    kernels_.at(b)->GetMesh().ClearScatter();
  }
}

template <class M>
void Local<M>::Bcast(const std::vector<MIdx>& bb) {
  using OpConcat = typename M::OpCat;
  auto& vf = kernels_.at(bb[0])->GetMesh().GetBcast(); // pointers to broadcast

  // Check size is the same for all kernels
  for (auto& b : bb) {
    auto& v = kernels_.at(b)->GetMesh().GetBcast(); // pointers to broadcast
    if (v.size() != vf.size()) {
      throw std::runtime_error(
          "Bcast: v.size()=" + std::to_string(v.size()) +
          " != vf.size()=" + std::to_string(vf.size()));
    }
  }

  for (size_t i = 0; i < vf.size(); ++i) {
    if (OpConcat* o = dynamic_cast<OpConcat*>(vf[i].get())) {
      std::vector<char> buf = o->Neutral(); // buffer

      // read from root block
      for (auto& b : bb) {
        auto& m = kernels_.at(b)->GetMesh();
        if (m.IsRoot()) {
          auto& v = m.GetBcast();
          OpConcat* ob = dynamic_cast<OpConcat*>(v[i].get());
          ob->Append(buf);
        }
      }

      // write to all blocks
      for (auto& b : bb) {
        auto& m = kernels_.at(b)->GetMesh();
        auto& v = m.GetBcast();
        OpConcat* ob = dynamic_cast<OpConcat*>(v[i].get());
        ob->Set(buf);
      }
    } else {
      throw std::runtime_error("Bcast: Unknown M::Op instance");
    }
  }

  // Clear bcast requests
  for (auto& b : bb) {
    auto& k = *kernels_.at(b);
    auto& m = k.GetMesh();
    m.ClearBcast();
  }
}

template <class M>
void Local<M>::DumpWrite(const std::vector<MIdx>& bb) {
  auto& m = kernels_.at(bb[0])->GetMesh();
  if (m.GetDump().size()) {
    std::string df = var.String["dumpformat"];
    if (df == "default") {
      df = "vtk";
    }

    if (df == "vtk") {
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
              on.second, gm, [this, k](IdxCell i) { return buf_[k][i]; }));
          k += on.first->GetSize();
          if (on.first->GetSize() != 1) {
            throw std::runtime_error("DumpWrite(): Support only size 1");
          }
        }

        oser_.reset(new output::SerVtkStruct<M>(v, "title", "p", gm));
      }

      // TODO: Check no change in comm list between time steps
      //       (otherwise oser_ needs reinitialization)

      oser_->Write(frame_, "title"); // TODO: t instead of stage_
      std::cerr << "Dump " << frame_ << ": format=" << df << std::endl;
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
  auto& gndc = gm.GetIndexCells();
  MIdx gs = gm.GetInBlockCells().GetSize();
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
  auto& bc = m.GetIndexCells();
  auto& gbc = gm.GetIndexCells();
  MIdx gs = gm.GetInBlockCells().GetSize();
  for (auto c : m.AllCells()) {
    auto w = bc.GetMIdx(c);
    // periodic
    for (int d = 0; d < 3; ++d) {
      if (per_[d]) {
        w[d] = (w[d] + gs[d]) % gs[d];
      }
    }
    // XXX: addhoc nan in halos
    if (MIdx(0) <= w && w < gs) {
      auto gc = gbc.GetIdx(w);
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
// Reads from buffer to Co
// o: instance of Co
// e: offset in buffer
// Returns:
// number of scalar fields written
template <class M>
size_t Local<M>::ReadBuffer(typename M::Co* o, size_t e, M& m) {
  if (auto od = dynamic_cast<typename M::CoFcs*>(o)) {
    return ReadBuffer(*od->field, e, m);
  }
  if (auto od = dynamic_cast<typename M::CoFcv*>(o)) {
    if (od->d == -1) {
      return ReadBuffer(*od->field, e, m);
    }
    return ReadBuffer(*od->field, od->d, e, m);
  }
  throw std::runtime_error("ReadBuffer: Unknown Co instance");
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
  auto& bc = m.GetIndexCells();
  auto& gbc = gm.GetIndexCells();
  for (auto c : m.Cells()) {
    auto w = bc.GetMIdx(c);
    auto gc = gbc.GetIdx(w);
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
  auto& bc = m.GetIndexCells();
  auto& gbc = gm.GetIndexCells();
  for (auto c : m.Cells()) {
    auto w = bc.GetMIdx(c);
    auto gc = gbc.GetIdx(w);
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
// Writes Co to buffer.
// o: instance of Co
// e: offset in buffer
// Returns:
// number of scalar fields written
template <class M>
size_t Local<M>::WriteBuffer(typename M::Co* o, size_t e, M& m) {
  if (auto od = dynamic_cast<typename M::CoFcs*>(o)) {
    return WriteBuffer(*od->field, e, m);
  }
  if (auto od = dynamic_cast<typename M::CoFcv*>(o)) {
    if (od->d == -1) {
      return WriteBuffer(*od->field, e, m);
    }
    return WriteBuffer(*od->field, od->d, e, m);
  }
  throw std::runtime_error("WriteBuffer: Unknown Co instance");
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
