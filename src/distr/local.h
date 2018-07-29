#pragma once

#include <vector>
#include <limits>
#include <map>
#include <mpi.h>
#include <stdexcept>

#include "dump/output.h"
#include "dump/output_paraview.h"
#include "distr.h"
#include "ilocal.h"

template <class KF>
class Local : public DistrMesh<KF> {
 public:
  using K = typename KF::K;
  using M = typename KF::M;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;

  Local(MPI_Comm comm, KF& kf, Vars& par);
  typename M::BlockCells GetGlobalBlock() const override;
  typename M::IndexCells GetGlobalIndex() const override;
  // Returns data field i from buffer defined on global mesh
  FieldCell<Scal> GetGlobalField(size_t i) override; 

 private:
  using MIdx = typename M::MIdx;
  using P = DistrMesh<KF>;

  using P::mk;
  using P::kf_;
  using P::par;
  using P::bs_;
  using P::es_;
  using P::hl_;
  using P::p_;
  using P::b_; 
  using P::stage_;
  using P::isroot_;
  using P::comm_;
  using P::ext_;
  using P::frame_;
  using P::dim;

  std::vector<FieldCell<Scal>> buf_; // buffer on mesh
  M gm; // global mesh
  std::unique_ptr<output::Ser> oser_; // output series
  std::vector<MyBlockInfo> bb_;

  size_t WriteBuffer(const FieldCell<Scal>& fc, size_t e, M& m);
  size_t WriteBuffer(const FieldCell<Vect>& f, size_t d, size_t e, M& m);
  size_t WriteBuffer(const FieldCell<Vect>& fc, size_t e, M& m);
  size_t WriteBuffer(typename M::Co* o, size_t e, M& m);
  void WriteBuffer(M& m);

  size_t ReadBuffer(FieldCell<Scal>& fc, size_t e,  M& m);
  size_t ReadBuffer(FieldCell<Vect>& fc, size_t d, size_t e,  M& m);
  size_t ReadBuffer(FieldCell<Vect>& fc, size_t e,  M& m);
  size_t ReadBuffer(typename M::Co* o, size_t e, M& m);
  void ReadBuffer(M& m);


  static M CreateMesh(MIdx bs, MIdx b, MIdx p, int es, Scal ext_);

  std::vector<MIdx> GetBlocks() override;
  void ReadBuffer(const std::vector<MIdx>& bb) override;
  void WriteBuffer(const std::vector<MIdx>& bb) override;
  void Reduce(const std::vector<MIdx>& bb) override;
  void DumpWrite(const std::vector<MIdx>& bb) override;
};

template <class KF>
auto Local<KF>::CreateMesh(MIdx bs, MIdx b, MIdx p, int /*es*/, Scal ext) -> M {
  // Init global mesh
  MIdx ms(bs); // block size 
  MIdx mb(b); // number of blocks
  MIdx mp(p); // number of PEs
  MIdx mm = mp * mb * ms; // total size in cells (without halos)

  Scal h = ext / std::max(std::max(mm[0], mm[1]), mm[2]);
  Vect d0(0); // origin coord
  Vect d1 = d0 + Vect(mm) * h;      // end coord
  Rect<Vect> d(d0, d1);

  MIdx o(0); // origin index
  std::cout 
    << "o=" << o 
    << " dom=" << d0 << "," << d1 
    << " h=" << h
    << std::endl;

  return InitUniformMesh<M>(d, o, mm, 0, true, mm);
}

template <class KF>
Local<KF>::Local(MPI_Comm comm, KF& kf, Vars& par) 
  : DistrMesh<KF>(comm, kf, par)
  , buf_(es_)
  , gm(CreateMesh(bs_, b_, p_, es_, ext_))
{

  // Resize buffer for mesh
  for (auto& u : buf_) {
    u.Reinit(gm);
  }

  // Fill block info
  MIdx ms(bs_); // block size 
  MIdx mb(b_[0], b_[1], b_[2]); // number of blocks
  MIdx mp(p_[0], p_[1], p_[2]); // number of PEs
  GBlockCells<3> bc(mb * mp);
  Scal h = (gm.GetNode(IdxNode(1)) - gm.GetNode(IdxNode(0)))[0];
  assert(h > 0);
  std::cerr << "h from gm = " << h << std::endl;
  for (MIdx i : bc) {
    MyBlockInfo b;
    IdxNode n = gm.GetIndexNodes().GetIdx(i * ms);
    Vect o = gm.GetNode(n);
    std::cerr << "o=" << o << " n=" << n.GetRaw() <<  " i=" << i << std::endl;
    for (int q = 0; q < 3; ++q) {
      b.index[q] = i[q];
      b.origin[q] = o[q];
      b.bs[q] = bs_[q];
    }
    b.h_gridpoint = h;
    b.ptrBlock = nullptr;
    b.hl = hl_;
    MIdx gs = p_ * b_ * bs_; // global size
    for (int j = 0; j < 3; ++j) {
      b.gs[j] = gs[j];
    }
    bb_.push_back(b);
  }

  isroot_ = true; // XXX: overwrite isroot_

  bool islead = true;
  for (auto& b : bb_) {
    MIdx d(b.index);
    b.isroot = (d == MIdx(0));
    b.islead = islead;
    islead = false;
  }

  this->MakeKernels(bb_);
}

template <class KF>
auto Local<KF>::GetBlocks() -> std::vector<MIdx> {
  // Put blocks to map by index 
  std::vector<MIdx> bb;
  for (auto e : bb_) {
    MIdx b(e.index);
    bb.push_back(b);
  }

  return bb;
}

template <class KF>
void Local<KF>::ReadBuffer(const std::vector<MIdx>& bb) {
  for (auto& b : bb) {
    auto& m = mk.at(b)->GetMesh();
    ReadBuffer(m);
  }
}

template <class KF>
void Local<KF>::WriteBuffer(const std::vector<MIdx>& bb) {
  for (auto& b : bb) {
    auto& m = mk.at(b)->GetMesh();
    WriteBuffer(m);
  }
}

template <class KF>
void Local<KF>::Reduce(const std::vector<MIdx>& bb) {
  using OpS = typename M::OpS;
  using OpSI = typename M::OpSI;
  using OpCat = typename M::OpCat;
  auto& f = *mk.at(bb[0]); // first kernel
  auto& mf = f.GetMesh();
  auto& vf = mf.GetReduce();  // reduce requests

  // Check size is the same for all kernels
  for (auto& b : bb) {
    auto& v = mk.at(b)->GetMesh().GetReduce();  // reduce requests
    if (v.size() != vf.size()) {
      throw std::runtime_error("Reduce: v.size() != vf.size()");
    }
  }

  // TODO: Check operation is the same for all kernels
  // TODO: avoid code duplication

  for (size_t i = 0; i < vf.size(); ++i) {
    if (OpS* o = dynamic_cast<OpS*>(vf[i].get())) {
      // Reduction on Scal
      
      auto r = o->Neut(); // result

      // Reduce over all blocks 
      for (auto& b : bb) {
        auto& v = mk.at(b)->GetMesh().GetReduce(); 
        OpS* ob = dynamic_cast<OpS*>(v[i].get());
        ob->Append(r);
      }

      // Write results to all blocks 
      for (auto& b : bb) {
        auto& v = mk.at(b)->GetMesh().GetReduce(); 
        OpS* ob = dynamic_cast<OpS*>(v[i].get());
        ob->Set(r);
      }
    } else if (OpSI* o = dynamic_cast<OpSI*>(vf[i].get())) {
      // Reduction on std::pair<Scal, int>
      
      auto r = o->Neut(); // result

      // Reduce over all blocks 
      for (auto& b : bb) {
        auto& v = mk.at(b)->GetMesh().GetReduce(); 
        OpSI* ob = dynamic_cast<OpSI*>(v[i].get());
        ob->Append(r);
      }

      // Write results to all blocks 
      for (auto& b : bb) {
        auto& v = mk.at(b)->GetMesh().GetReduce(); 
        OpSI* ob = dynamic_cast<OpSI*>(v[i].get());
        ob->Set(r);
      }
    } else if (OpCat* o = dynamic_cast<OpCat*>(vf[i].get())) {
      // Concatenation of std::vector<T>
      
      auto r = o->Neut(); // result

      GBlock<size_t, dim> qp(p_); 
      GBlock<size_t, dim> qb(b_); 

      // Reduce over all blocks 
      for (auto wp : qp) { // same ordering as with Cubism
        for (auto wb : qb) {
          auto w = b_ * wp + wb;
          auto& v = mk.at(w)->GetMesh().GetReduce(); 
          OpCat* ob = dynamic_cast<OpCat*>(v[i].get());
          ob->Append(r);
        }
      }

      // Write results to root block 
      for (auto& b : bb) {
        auto& m = mk.at(b)->GetMesh();
        if (m.IsRoot()) {
          auto& v = m.GetReduce(); 
          OpCat* ob = dynamic_cast<OpCat*>(v[i].get());
          ob->Set(r); 
        }
      }
    } else {
      throw std::runtime_error("Reduce: Unknown M::Op implementation");
    }
  }

  // Clear reduce requests
  for (auto& b : bb) {
    auto& k = *mk.at(b); 
    auto& m = k.GetMesh();
    m.ClearReduce();
  }
}

template <class KF>
void Local<KF>::DumpWrite(const std::vector<MIdx>& bb) {
  auto& m = mk.at(bb[0])->GetMesh();
  if (m.GetDump().size()) {
    std::string df = par.String["dumpformat"];
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
          v.emplace_back(
              new output::OutFldFunc<Scal, IdxCell, M>(
                  on.second, gm, [this,k](IdxCell i) { return buf_[k][i]; }));
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
template <class KF>
size_t Local<KF>::ReadBuffer(FieldCell<Scal>& fc, size_t e,  M& m) {
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
      w[d] = (w[d] + gs[d]) % gs[d];
    }
    auto gc = gndc.GetIdx(w);
    fc[c] = buf_[e][gc];
  }
  return 1;
}

// Reads from buffer to component of vector field [a]
// fc: field
// d: component (0,1,2)
// e: offset in buffer
// Returns:
// number of scalar fields read
template <class KF>
size_t Local<KF>::ReadBuffer(FieldCell<Vect>& fc, size_t d, size_t e,  M& m) {
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
      w[d] = (w[d] + gs[d]) % gs[d];
    }
    auto gi = gbc.GetIdx(w);
    fc[c][d] = buf_[e][gi];
  }
  return 1;
}

// Reads from buffer to all components of vector field [a]
// fc: field
// e: offset in buffer
// Returns:
// number of scalar fields read
template <class KF>
size_t Local<KF>::ReadBuffer(FieldCell<Vect>& fc, size_t e,  M& m) {
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
template <class KF>
size_t Local<KF>::ReadBuffer(typename M::Co* o, size_t e, M& m) {
  if (auto od = dynamic_cast<typename M::CoFcs*>(o)) {
    return ReadBuffer(*od->f, e, m);
  } else if (auto od = dynamic_cast<typename M::CoFcv*>(o)) {
    if (od->d == -1) {
      return ReadBuffer(*od->f, e, m);
    } 
    return ReadBuffer(*od->f, od->d, e, m);
  } 
  throw std::runtime_error("ReadBuffer: Unknown Co instance");
  return 0;
}
template <class KF>
void Local<KF>::ReadBuffer(M& m) {
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
template <class KF>
size_t Local<KF>::WriteBuffer(const FieldCell<Scal>& fc, size_t e, M& m) {
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
template <class KF>
size_t Local<KF>::WriteBuffer(const FieldCell<Vect>& fc, 
                              size_t d, size_t e, M& m) {
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
template <class KF>
size_t Local<KF>::WriteBuffer(const FieldCell<Vect>& fc, size_t e, M& m) {
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
template <class KF>
size_t Local<KF>::WriteBuffer(typename M::Co* o, size_t e, M& m) {
  if (auto od = dynamic_cast<typename M::CoFcs*>(o)) {
    return WriteBuffer(*od->f, e, m);
  } else if (auto od = dynamic_cast<typename M::CoFcv*>(o)) {
    if (od->d == -1) {
      return WriteBuffer(*od->f, e, m);
    } 
    return WriteBuffer(*od->f, od->d, e, m);
  }
  throw std::runtime_error("WriteBuffer: Unknown Co instance");
  return 0;
}
template <class KF>
void Local<KF>::WriteBuffer(M& m) {
  size_t e = 0;
  for (auto& o : m.GetComm()) {
    e += WriteBuffer(o.get(), e, m);
  }
  for (auto& on : m.GetDump()) {
    e += WriteBuffer(on.first.get(), e, m);
  }
}

template <class KF>
auto Local<KF>::GetGlobalBlock() const -> typename M::BlockCells {
  return gm.GetInBlockCells();
}

template <class KF>
auto Local<KF>::GetGlobalIndex() const -> typename M::IndexCells {
  return gm.GetIndexCells();
}

template <class KF>
auto Local<KF>::GetGlobalField(size_t i) -> FieldCell<Scal> {
  return buf_[i];
}

