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

  Local(MPI_Comm comm, KF& kf, Vars& par);
  typename M::BlockCells GetGlobalBlock() const override;
  // Returns data field i from buffer defined on global mesh
  FieldCell<Scal> GetGlobalField(size_t i) const override; 

 private:
  using MIdx = typename M::MIdx;
  using Vect = typename M::Vect;
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

  std::vector<FieldCell<Scal>> buf_; // buffer on mesh
  M gm; // global mesh
  std::unique_ptr<output::Session> session_;
  std::vector<MyBlockInfo> bb_;

  void ReadBuffer(M& m);
  void WriteBuffer(M& m);
  static M CreateMesh(MIdx bs, MIdx b, MIdx p, int es, Scal ext_);
  std::vector<MIdx> GetBlocks() override;
  void ReadBuffer(const std::vector<MIdx>& bb) override;
  void WriteBuffer(const std::vector<MIdx>& bb) override;
  void Reduce(const std::vector<MIdx>& bb) override;
  void DumpWrite(const std::vector<MIdx>& bb) override;
};

template <class KF>
auto Local<KF>::CreateMesh(MIdx bs, MIdx b, MIdx p, int es, Scal ext) -> M {
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

  return InitUniformMesh<M>(d, o, mm, 0);
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
    IdxNode n = gm.GetBlockNodes().GetIdx(i * ms);
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
    auto& k = *mk.at(b); // kernel
    auto& m = k.GetMesh();

    ReadBuffer(m);
  }
}

template <class KF>
void Local<KF>::WriteBuffer(const std::vector<MIdx>& bb) {
  for (auto& b : bb) {
    auto& k = *mk.at(b); // kernel
    auto& m = k.GetMesh();
    WriteBuffer(m);
  }
}

template <class KF>
void Local<KF>::Reduce(const std::vector<MIdx>& bb) {
  using Op = typename M::Op;
  using OpS = typename M::OpS;
  using OpSI = typename M::OpSI;
  auto& f = *mk.at(bb[0]); // first kernel
  auto& mf = f.GetMesh();
  auto& vf = mf.GetReduce();  // pointers to reduce

  // Check size is the same for all kernels
  for (auto& b : bb) {
    // TODO: collapse next 2 lines
    auto& k = *mk.at(b); // kernel
    auto& m = k.GetMesh();
    auto& v = m.GetReduce();  // pointers to reduce
  }

  // TODO: Check operation is the same for all kernels

  for (size_t i = 0; i < vf.size(); ++i) {
    if (OpS* o = dynamic_cast<OpS*>(vf[i].get())) {
      auto r = o->Neut(); // result

      // Reduce over all blocks on current rank
      for (auto& b : bb) {
        auto& v = mk.at(b)->GetMesh().GetReduce(); 
        OpS* ob = dynamic_cast<OpS*>(v[i].get());
        ob->Append(r);
      }

      // Write results to all blocks on current rank
      for (auto& b : bb) {
        auto& v = mk.at(b)->GetMesh().GetReduce(); 
        OpS* ob = dynamic_cast<OpS*>(v[i].get());
        ob->Set(r);
      }
    } else if (OpSI* o = dynamic_cast<OpSI*>(vf[i].get())) {
      auto r = o->Neut(); // result

      // Reduce over all blocks on current rank
      for (auto& b : bb) {
        auto& v = mk.at(b)->GetMesh().GetReduce(); 
        OpSI* ob = dynamic_cast<OpSI*>(v[i].get());
        ob->Append(r);
      }

      // Write results to all blocks on current rank
      for (auto& b : bb) {
        auto& v = mk.at(b)->GetMesh().GetReduce(); 
        OpSI* ob = dynamic_cast<OpSI*>(v[i].get());
        ob->Set(r);
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
      if (!session_) {
        // TODO: check all blocks are same as first
        output::Content c;
        size_t k = 0; // offset in buffer
        auto& vcm = m.GetComm(); // comm and dump added by DumpComm
        auto& vd = m.GetDump(); // dump
        // Skip comm 
        for (size_t i = 0; i < vcm.size() - vd.size(); ++i) {
          k += vcm[i]->GetSize();
        }
        // Write dump
        for (auto& o : m.GetDump()) {
          c.emplace_back(
              new output::EntryFunction<Scal, IdxCell, M>(
                  o.second, gm, [this,k](IdxCell i) { return buf_[k][i]; }));
          k += o.first->GetSize();
          if (o.first->GetSize() != 1) {
            throw std::runtime_error("DumpWrite(): Support only size 1");
          }
        }

        session_.reset(new output::SessionParaviewStructured<M>(
              c, "title", "p" /*filename*/, gm));
      }

      // TODO: Check no change in comm list between time steps
      //       (otherwise session_ needs reinitialization)

      std::cerr << "Dump " << frame_ << ": format=" << df << std::endl;
      ++frame_;
      session_->Write(stage_ * 1., "title"); // TODO: t instead of stage_
    } else {
      P::DumpWrite(bb);
    }
  }
}


template <class KF>
void Local<KF>::ReadBuffer(M& m) {
  int e = 0; // buffer field idx

  auto& bc = m.GetBlockCells();
  auto& gbc = gm.GetBlockCells();
  MIdx gs = gbc.GetDimensions();
  for (auto& u : m.GetComm()) {
    if (auto ud = dynamic_cast<typename M::CoFcs*>(u.get())) {
      for (auto c : m.AllCells()) {
        auto w = bc.GetMIdx(c);
        // periodic
        for (int d = 0; d < 3; ++d) {
          w[d] = (w[d] + gs[d]) % gs[d];
        }
        auto gi = gbc.GetIdx(w);
        (*ud->f)[c] = buf_[e][gi];
      }
    }
    ++e;
  }
}

template <class KF>
void Local<KF>::WriteBuffer(M& m) {
  using MIdx = typename M::MIdx;

  // Check buffer has enough space for all fields
  assert(m.GetComm().size() <= buf_.size() && "Too many fields for Comm()");

  int e = 0; // buffer field idx

  auto& bc = m.GetBlockCells();
  auto& gbc = gm.GetBlockCells();
  for (auto u : m.GetComm()) {
    if (auto ud = dynamic_cast<typename M::CoFcs*>(u.get())) {
      for (auto c : m.Cells()) {
        auto w = bc.GetMIdx(c); 
        auto gc = gbc.GetIdx(w); 
        buf_[e][gc] = (*ud->f)[c];
      }
    }
    ++e;
  }
}

template <class KF>
auto Local<KF>::GetGlobalBlock() const -> typename M::BlockCells {
  return gm.GetBlockCells();
}

template <class KF>
auto Local<KF>::GetGlobalField(size_t i) const -> FieldCell<Scal> {
  return buf_[i];
}

