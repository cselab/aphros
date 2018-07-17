#pragma once

#include <memory>
#include <limits>
#include <map>
#include <stdexcept>
#include <mpi.h>

#include "distr.h"
#include "icubism.h"
#include "dump/dumper.h"

#include "Cubism/BlockInfo.h"
#include "Cubism/Grid.h"
#include "Cubism/GridMPI.h"
#include "Cubism/BlockLab.h"
#include "Cubism/BlockLabMPI.h"
#include "Cubism/StencilInfo.h"
#include "Cubism/HDF5Dumper_MPI.h"

// Hide implementation and avoid collision with GBlock
// TODO: rename cubism_impl::GBlock
namespace cubism_impl {

// Static parameters for Cubism
// bx, by, bz - block size
// es - Elem size
template <class Scal_, size_t bx_, size_t by_, size_t bz_, size_t es_>
struct GPar {
  using Scal = Scal_;
  static const size_t bx = bx_;
  static const size_t by = by_;
  static const size_t bz = bz_;
  static const size_t es = es_;
};

template <class Scal_, int es_>
struct GElem {
  using Scal = Scal_;

  static const size_t es = es_;

  Scal a[es];

  void init(Scal val) { 
    for (size_t i = 0; i < es; ++i) {
      a[i] = val;
    }
  }
  void clear() {
    init(0);
  }
  GElem& operator=(const GElem&) = default;
};

template <class Par_>
struct GBlock {
  using Par = Par_;
  using Scal = typename Par::Scal;
  static const size_t bx = Par::bx;
  static const size_t by = Par::by;
  static const size_t bz = Par::bz;
  static const size_t es = Par::es;

  static const int n = bx * by * bz;

  // required by Cubism
  static const int sizeX = bx;
  static const int sizeY = by;
  static const int sizeZ = bz;

  using Elem = GElem<Scal, es>;
  using ElementType = Elem;
  using element_type = Elem;

  // floats per element
  static const int fe = sizeof(Elem) / sizeof(Scal);
  static_assert(fe == es, "Block: fe != es");

  Elem __attribute__((__aligned__(_ALIGNBYTES_))) data[bz][by][bx];

  Scal __attribute__((__aligned__(_ALIGNBYTES_))) tmp[bz][by][bx][fe];

  void clear_data() {
    Elem* e = &data[0][0][0];
    for(int i = 0; i < n; ++i) {
      e[i].clear();
    }
  }

  void clear_tmp() {
    Scal* t = &tmp[0][0][0][0];
    for(int i = 0; i < n * fe; ++i) {
      t[i] = 0;
    }
  }

  void clear() {
    clear_data();
    clear_tmp();
  }

  inline Elem& operator()(int x, int y=0, int z=0) {
    assert(0 <= x && x < bx);
    assert(0 <= y && y < by);
    assert(0 <= z && z < bz);

    return data[z][y][x];
  }
};


// Par - instance of GPar
template<class Par, template<typename X> class A=std::allocator>
class LabPer : public BlockLab<GBlock<Par>, A>
{
  using Block = GBlock<Par>;
  using ElementTypeBlock = typename Block::Elem;

 public:
  virtual inline std::string name() const { return "LabPer"; }
  bool is_xperiodic() {return true;}
  bool is_yperiodic() {return true;}
  bool is_zperiodic() {return true;}

  // TODO: remove
  LabPer() : BlockLab<Block,A>() {}
};


// Par - instance of GPar
template <class Par>
using GLab = BlockLabMPI<LabPer<Par, std::allocator>>;

// B - block type
template <class B>
using GGrid = GridMPI<Grid<B, std::allocator>>;

// Par - instance of GPar
// KF - instance of KernelMeshFactory
template <class Par, class KF>
class Cubism : public DistrMesh<KF> {
 public:
  using K = typename KF::K;
  using M = typename KF::M;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;

  Cubism(MPI_Comm comm, KF& kf, Vars& par);
  typename M::BlockCells GetGlobalBlock() const override;
  FieldCell<Scal> GetGlobalField(size_t i) override; 

 private:
  using Lab = GLab<Par>;
  using Block = GBlock<Par>;
  using Grid = GGrid<Block>;
  using Elem = typename Block::Elem;
  using Synch = typename Grid::Synch;

  using P = DistrMesh<KF>; // parent
  using MIdx = typename M::MIdx;

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

  Grid g_;
  struct S { // cubism [s]tate
    Synch* s;
    std::unique_ptr<Lab> l;
    std::map<MIdx, BlockInfo, typename MIdx::LexLess> mb;
  };
  S s_;

  // hl: number of halo cells from each side
  // cs: number of fields for communication
  static StencilInfo GetStencil(int hl, int cs) {
    const int a = -hl;
    const int b = hl + 1;
    assert(cs <= Elem::es);
    if (cs == 0) {
      return StencilInfo(a,a,a,b,b,b, true, cs);
    } else if (cs == 1) {
      return StencilInfo(a,a,a,b,b,b, true, cs, 0);
    } else if (cs == 2) {
      return StencilInfo(a,a,a,b,b,b, true, cs, 0,1);
    } else if (cs == 3) {
      return StencilInfo(a,a,a,b,b,b, true, cs, 0,1,2);
    } else if (cs == 4) {
      return StencilInfo(a,a,a,b,b,b, true, cs, 0,1,2,3);
    } else if (cs == 5) {
      return StencilInfo(a,a,a,b,b,b, true, cs, 0,1,2,3,4);
    } else if (cs == 6) {
      return StencilInfo(a,a,a,b,b,b, true, cs, 0,1,2,3,4,5);
    } else if (cs == 7) {
      return StencilInfo(a,a,a,b,b,b, true, cs, 0,1,2,3,4,5,6);
    } else if (cs == 8) {
      return StencilInfo(a,a,a,b,b,b, true, cs, 0,1,2,3,4,5,6,7);
    } else {
      std::cerr << "GetStencil(): Unknown cs=" << cs << std::endl;
      assert(false);
    }
  }
  // Convert Cubism BlockInfo to MyBlockInfo
  static std::vector<MyBlockInfo> GetBlocks(
      const std::vector<BlockInfo>&, MIdx bs, size_t hl);

  // Reads from buffer to scalar field [a]
  // fc: field
  // l: lab
  // e: offset in buffer, 0 <= e < Elem::es
  // Returns:
  // number of scalar fields read
  size_t ReadBuffer(FieldCell<Scal>& fc, Lab& l, size_t e,  M& m) {
    if (e >= Elem::es) {
      throw std::runtime_error("ReadBuffer: Too many fields for Comm()");
    }
    auto& bc = m.GetBlockCells();
    for (auto c : m.AllCells()) {
      auto w = bc.GetMIdx(c) - MIdx(hl_) - bc.GetBegin();
      fc[c] = l(w[0], w[1], w[2]).a[e];
    }
    return 1;
  }
  // Reads from buffer to component of vector field [a]
  // fc: field
  // d: component (0,1,2)
  // l: lab
  // e: offset in buffer, 0 <= e < Elem::es
  // Returns:
  // number of scalar fields read
  size_t ReadBuffer(FieldCell<Vect>& fc, size_t d, Lab& l, size_t e,  M& m) {
    if (e >= Elem::es) {
      throw std::runtime_error("ReadBuffer: Too many fields for Comm()");
    }
    if (d >= Vect::dim) {
      throw std::runtime_error("ReadBuffer: d >= Vect::dim");
    }
    auto& bc = m.GetBlockCells();
    for (auto c : m.AllCells()) {
      auto w = bc.GetMIdx(c) - MIdx(hl_) - bc.GetBegin();
      fc[c][d] = l(w[0], w[1], w[2]).a[e];
    }
    return 1;
  }
  // Reads from buffer to all components of vector field [a]
  // fc: field
  // l: lab
  // e: offset in buffer, 0 <= e < Elem::es
  // Returns:
  // number of scalar fields read
  size_t ReadBuffer(FieldCell<Vect>& fc, Lab& l, size_t e,  M& m) {
    for (size_t d = 0; d < Vect::dim; ++d) {
      e += ReadBuffer(fc, d, l, e, m);
    }
    return Vect::dim;
  }
  // Reads from buffer to Co
  // o: instance of Co
  // l: lab
  // e: offset in buffer, 0 <= e < Elem::es
  // Returns:
  // number of scalar fields written
  size_t ReadBuffer(typename M::Co* o, Lab& l, size_t e, M& m) {
    if (auto od = dynamic_cast<typename M::CoFcs*>(o)) {
      return ReadBuffer(*od->f, l, e, m);
    } else if (auto od = dynamic_cast<typename M::CoFcv*>(o)) {
      if (od->d == -1) {
        return ReadBuffer(*od->f, l, e, m);
      } 
      return ReadBuffer(*od->f, od->d, l, e, m);
    } 
    throw std::runtime_error("ReadBuffer: Unknown Co instance");
    return 0;
  }
  void ReadBuffer(M& m, Lab& l) {
    size_t e = 0;
    for (auto& o : m.GetComm()) {
      e += ReadBuffer(o.get(), l, e, m);
    }
  }
  // Writes scalar field to buffer [i].
  // fc: scalar field
  // b: block 
  // e: offset in buffer, 0 <= e < Elem::es
  // Returns:
  // number of scalar fields written
  size_t WriteBuffer(const FieldCell<Scal>& fc, Block& b, size_t e, M& m) {
    if (e >= Elem::es) {
      throw std::runtime_error("WriteBuffer: Too many fields for Comm()");
    }
    auto& bc = m.GetBlockCells();
    for (auto c : m.Cells()) {
      auto w = m.GetBlockCells().GetMIdx(c) - MIdx(hl_) - bc.GetBegin();
      b.data[w[2]][w[1]][w[0]].a[e] = fc[c];
    }
    return 1;
  }
  // Writes component of vector field to buffer [i].
  // fc: vector field
  // d: component (0,1,2)
  // b: block 
  // e: offset in buffer, 0 <= e < Elem::es
  // Returns:
  // number of scalar fields written
  size_t WriteBuffer(const FieldCell<Vect>& fc, 
                     size_t d, Block& b, size_t e, M& m) {
    if (e >= Elem::es) {
      throw std::runtime_error("WriteBuffer: Too many fields for Comm()");
    }
    if (d >= Vect::dim) {
      throw std::runtime_error("ReadBuffer: d >= Vect::dim");
    }
    auto& bc = m.GetBlockCells();
    for (auto c : m.Cells()) {
      auto w = m.GetBlockCells().GetMIdx(c) - MIdx(hl_) - bc.GetBegin();
      b.data[w[2]][w[1]][w[0]].a[e] = fc[c][d];
    }
    return 1;
  }
  // Writes all components of vector field to buffer [i].
  // fc: vector field
  // b: block 
  // e: offset in buffer, 0 <= e < Elem::es
  // Returns:
  // number of scalar fields written
  size_t WriteBuffer(const FieldCell<Vect>& fc, Block& b, size_t e, M& m) {
    for (size_t d = 0; d < Vect::dim; ++d) {
      e += WriteBuffer(fc, d, b, e, m);
    }
    return Vect::dim;
  }
  // Writes Co to buffer.
  // o: instance of Co
  // b: block 
  // e: offset in buffer, 0 <= e < Elem::es
  // Returns:
  // number of scalar fields written
  size_t WriteBuffer(typename M::Co* o, Block& b, size_t e, M& m) {
    if (auto od = dynamic_cast<typename M::CoFcs*>(o)) {
      return WriteBuffer(*od->f, b, e, m);
    } else if (auto od = dynamic_cast<typename M::CoFcv*>(o)) {
      if (od->d == -1) {
        return WriteBuffer(*od->f, b, e, m);
      } 
      return WriteBuffer(*od->f, od->d, b, e, m);
    }
    throw std::runtime_error("WriteBuffer: Unknown Co instance");
    return 0;
  }
  void WriteBuffer(M& m, Block& b) {
    size_t e = 0; 
    for (auto& o : m.GetComm()) {
      e += WriteBuffer(o.get(), b, e, m);
    }
    for (auto& on : m.GetDump()) {
      e += WriteBuffer(on.first.get(), b, e, m);
    }
  }

  std::vector<MIdx> GetBlocks() override;
  void ReadBuffer(const std::vector<MIdx>& bb) override;
  void WriteBuffer(const std::vector<MIdx>& bb) override;
  void Reduce(const std::vector<MIdx>& bb) override;
  void DumpWrite(const std::vector<MIdx>& bb) override;
};


// B_ - instance of GBlock
template <class B_, int ID>
struct StreamHdf {
  using B = B_;
  using Scal = typename B::Scal;
  using T = typename B::Elem;

  // Required by Cubism
  static const std::string NAME;
  static const std::string EXT;
  static const int NCHANNELS = 1;
  static const int CLASS = 0;

  B& b;

  StreamHdf(B& b): b(b) {}

  // write
  void operate(const int ix, const int iy, const int iz, Scal out[0]) const
  {
    const T& in = *((T*)&b.data[iz][iy][ix].a[0]);
    out[0] = in.a[ID];
  }

  // read
  void operate(const Scal out[0], const int ix, const int iy, const int iz) const
  {
    T& in = *((T*)&b.data[iz][iy][ix].a[0]);
    in.a[ID] = out[0];
  }

  static const char * getAttributeName() { return "Scalar"; }
};

// TODO: somehow remove global parameter ID
template <class B_>
struct StreamHdfDyn {
  using B = B_;
  using Scal = typename B::Scal;
  using T = typename B::Elem;

  // Required by Cubism
  static std::string NAME;
  static const std::string EXT;
  static const int NCHANNELS = 1;
  static const int CLASS = 0;

  // Used as global variable
  static int ID;

  B& b;

  StreamHdfDyn(B& b): b(b) {}

  // write
  void operate(const int ix, const int iy, const int iz, Scal out[0]) const
  {
    const T& in = *((T*)&b.data[iz][iy][ix].a[0]);
    out[0] = in.a[ID];
  }

  // read
  void operate(const Scal out[0], const int ix, const int iy, const int iz) const
  {
    T& in = *((T*)&b.data[iz][iy][ix].a[0]);
    in.a[ID] = out[0];
  }

  static const char * getAttributeName() { return "Scalar"; }
};

// Class with field 'stencil' needed for SynchronizerMPI::sync(Processing)
struct FakeProc {
  StencilInfo stencil;
  explicit FakeProc(StencilInfo si) 
    : stencil(si)
  {}
};

template <class Par, class KF>
std::vector<MyBlockInfo> Cubism<Par, KF>::GetBlocks(
    const std::vector<BlockInfo>& cc, MIdx bs, size_t hl) {
  std::vector<MyBlockInfo> bb;
  for(size_t i = 0; i < cc.size(); i++) {
    const BlockInfo& c = cc[i];
    MyBlockInfo b;
    for (int j = 0; j < 3; ++j) {
      b.index[j] = c.index[j];
      b.origin[j] = c.origin[j];
      b.bs[j] = bs[j];
    }
    b.h_gridpoint = c.h_gridpoint;
    b.ptrBlock = c.ptrBlock;
    b.hl = hl;
    bb.push_back(b);
  }
  return bb;
}


template <class Par, class KF>
Cubism<Par, KF>::Cubism(MPI_Comm comm, KF& kf, Vars& par) 
  : DistrMesh<KF>(comm, kf, par)
  , g_(p_[0], p_[1], p_[2], b_[0], b_[1], b_[2], ext_, comm)
{
  assert(bs_ == MIdx(Block::bx, Block::by, Block::bz));

  int r;
  MPI_Comm_rank(comm, &r);
  isroot_ = (0 == r);  // XXX: overwrite isroot_

  std::vector<BlockInfo> cc = g_.getBlocksInfo(); // [c]ubism block info
  std::vector<MyBlockInfo> ee = GetBlocks(cc, bs_, hl_);

  bool islead = true;
  for (auto& e : ee) {
    MIdx d(e.index);
    e.isroot = (d == MIdx(0));
    e.islead = islead;
    islead = false;
  }

  this->MakeKernels(ee);

  comm_ = g_.getCartComm(); // XXX: overwrite comm_
}

template <class Par, class KF>
auto Cubism<Par, KF>::GetBlocks() -> std::vector<MIdx> {
  MPI_Barrier(comm_);

  // Get all blocks
  std::vector<BlockInfo> cc = g_.getBlocksInfo(); // all blocks
  auto& m = mk.at(MIdx(cc[0].index))->GetMesh();

  size_t cs = 0; // comm size
  for (auto& o : m.GetComm()) {
    cs += o->GetSize();
  }

  std::vector<BlockInfo> aa;
  // Perform communication if necessary or s_.l not initialized
  if (cs > 0 || !s_.l) { 
    // 1. Exchange halos in buffer mesh.
    // max(cs, 1) to prevent forbidden call with zero components
    FakeProc fp(GetStencil(hl_, std::max<size_t>(cs, 1))); 
    Synch& s = g_.sync(fp); 

    s_.l.reset(new Lab);
    s_.l->prepare(g_, s);   // allocate memory for lab cache

    MPI_Barrier(comm_);

    // Do communication and get all blocks
    aa = s.avail(cc.size());
  } else {
    aa = cc;
  }

  if (aa.size() != cc.size()) {
    std::cerr 
        << "aa.size()=" << aa.size()
        << " != "
        << "cc.size()=" << cc.size()
        << std::endl;
    assert(false);
  }

  // Create vector of indices and save block info to map 
  std::vector<MIdx> bb;
  s_.mb.clear();
  for (auto a : aa) {
    MIdx b(a.index);
    s_.mb.emplace(b, a);
    bb.push_back(b);
  }

  return bb;
}

template <class Par, class KF>
void Cubism<Par, KF>::ReadBuffer(const std::vector<MIdx>& bb) {
  for (auto& b : bb) {
    auto& k = *mk.at(b); // kernel
    auto& m = k.GetMesh();
    s_.l->load(s_.mb[b], stage_);
    ReadBuffer(m, *s_.l);
  }
}

template <class Par, class KF>
void Cubism<Par, KF>::WriteBuffer(const std::vector<MIdx>& bb) {
  for (auto& b : bb) {
    auto& k = *mk.at(b); // kernel
    auto& m = k.GetMesh();
    WriteBuffer(m, *(typename Grid::BlockType*)s_.mb[b].ptrBlock);
  }
}

template <class Par, class KF>
void Cubism<Par, KF>::Reduce(const std::vector<MIdx>& bb) {
  using Op = typename M::Op;
  using OpS = typename M::OpS;
  using OpSI = typename M::OpSI;
  using OpCat = typename M::OpCat;
  auto& f = *mk.at(bb[0]); // first kernel
  auto& mf = f.GetMesh();
  auto& vf = mf.GetReduce();  // pointers to reduce

  // Check size is the same for all kernels
  for (auto& b : bb) {
    auto& k = *mk.at(b); // kernel
    auto& m = k.GetMesh();
    auto& v = m.GetReduce();  // pointers to reduce
    assert(v.size() == vf.size());
  }

  // TODO: Check operation is the same for all kernels
  // TODO: avoid code duplication

  for (size_t i = 0; i < vf.size(); ++i) {
    if (OpS* o = dynamic_cast<OpS*>(vf[i].get())) {
      auto r = o->Neut(); // result
      
      // Reduce over all blocks on current rank
      for (auto& b : bb) {
        auto& v = mk.at(b)->GetMesh().GetReduce(); 
        OpS* ob = dynamic_cast<OpS*>(v[i].get());
        ob->Append(r);
      }

      MPI_Op mo;
      if (dynamic_cast<typename M::OpSum*>(o)) {
        mo = MPI_SUM;
      } else if (dynamic_cast<typename M::OpProd*>(o)) {
        mo = MPI_PROD;
      } else if (dynamic_cast<typename M::OpMax*>(o)) {
        mo = MPI_MAX;
      } else if (dynamic_cast<typename M::OpMin*>(o)) {
        mo = MPI_MIN;
      } else {
        throw std::runtime_error("Reduce(): Can't find MPI_Op");
      }
      MPI_Datatype mt = (sizeof(Scal) == 8 ? MPI_DOUBLE : MPI_FLOAT);

      // Reduce over all ranks
      MPI_Allreduce(MPI_IN_PLACE, &r, 1, mt, mo, comm_); 

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

      MPI_Op mo;
      if (dynamic_cast<typename M::OpMinloc*>(o)) {
        mo = MPI_MINLOC;
      } else if (dynamic_cast<typename M::OpMaxloc*>(o)) {
        mo = MPI_MAXLOC;
      } else {
        throw std::runtime_error("Reduce(): Can't find MPI_Op");
      }

      MPI_Datatype mt = (sizeof(Scal) == 8 ? MPI_DOUBLE_INT : MPI_FLOAT_INT);

      // Reduce over all ranks
      MPI_Allreduce(MPI_IN_PLACE, &r, 1, mt, mo, comm_); 

      // Write results to all blocks on current rank
      for (auto& b : bb) {
        auto& v = mk.at(b)->GetMesh().GetReduce(); 
        OpSI* ob = dynamic_cast<OpSI*>(v[i].get());
        ob->Set(r);
      }
    } else if (OpCat* o = dynamic_cast<OpCat*>(vf[i].get())) {
      std::vector<char> r = o->Neut(); // result local
      
      // Reduce over all local blocks
      for (auto& b : bb) {
        auto& v = mk.at(b)->GetMesh().GetReduce(); 
        OpCat* ob = dynamic_cast<OpCat*>(v[i].get());
        ob->Append(r);
      }

      int s = r.size(); // size local

      if (isroot_) {
        int sc; // size of communicator
        MPI_Comm_size(comm_, &sc);

        std::vector<int> ss(sc); // size of r on all ranks

        // Gather ss
        MPI_Gather(&s, 1, MPI_INT, 
                   ss.data(), 1, MPI_INT, 
                   0, comm_);


        int sa = 0; // size all
        std::vector<int> oo = {0}; // offsets
        for (auto& q : ss) {
          sa += q;
          oo.push_back(oo.back() + q);
        }
        oo.pop_back();
        assert(ss.size() == oo.size());

        std::vector<char> ra(sa); // result all

        // Gather ra
        MPI_Gatherv(r.data(), r.size(), MPI_CHAR,
                    ra.data(), ss.data(), oo.data(), MPI_CHAR,
                    0, comm_);
       
        // Write results to root block 
        size_t cnt = 0;
        for (auto& b : bb) {
          auto& m = mk.at(b)->GetMesh();
          if (m.IsRoot()) {
            auto& v = m.GetReduce(); 
            OpCat* ob = dynamic_cast<OpCat*>(v[i].get());
            ob->Set(ra);
            ++cnt;
          }
        }
        assert(cnt == 1);
      } else {
        // Send s to root
        MPI_Gather(&s, 1, MPI_INT, nullptr, 0, MPI_INT, 0, comm_);

        // Send r to root
        MPI_Gatherv(r.data(), r.size(), MPI_CHAR,
                    nullptr, nullptr, nullptr, MPI_CHAR,
                    0, comm_);
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

template <class Par, class KF>
void Cubism<Par, KF>::DumpWrite(const std::vector<MIdx>& bb) {
  auto& m = mk.at(bb[0])->GetMesh();
  if (m.GetDump().size()) {
    std::string df = par.String["dumpformat"];
    if (df == "default") {
      df = "hdf";
    }

    if (df == "hdf") {
      size_t k = 0; // offset in buffer
      // Skip comm 
      for (auto& o : m.GetComm()) {
        k += o->GetSize();
      }
      // Write dump
      for (auto& on : m.GetDump()) {
        StreamHdfDyn<Block>::ID = k;
        StreamHdfDyn<Block>::NAME = on.second;
        auto fn = GetDumpName(on.second, "", frame_);
        DumpHDF5_MPI<Grid, StreamHdfDyn<Block>>(g_, frame_, frame_, fn);
        k += on.first->GetSize();
        if (on.first->GetSize() != 1) {
          throw std::runtime_error("DumpWrite(): Support only size 1");
        }
      }
      if (isroot_) {
        std::cerr << "Dump " << frame_ << ": format=" << df << std::endl;
      }
      ++frame_;
    } else {
      P::DumpWrite(bb);
    }
  }
}

template <class Par, class KF>
auto Cubism<Par, KF>::GetGlobalBlock() const -> typename M::BlockCells {
  return typename M::BlockCells(p_ * b_ * bs_);
}

template <class Par, class KF>
auto Cubism<Par, KF>::GetGlobalField(size_t e) -> FieldCell<Scal> {
  using BC = typename M::BlockCells;
  auto gbc = GetGlobalBlock();
  // collective, does actual communication
  auto bb = GetBlocks();
  FieldCell<Scal> fc; // tmp
  std::vector<Scal> v(bs_.prod()); // tmp
  MPI_Datatype mt = (sizeof(Scal) == 8 ? MPI_DOUBLE : MPI_FLOAT);
  BC bc(bs_); // cells
  BC gbb(p_ * b_);  // all blocks
  if (isroot_) {
    FieldCell<Scal> gfc(gbc); // result
    // Copy from blocks on root
    for (auto& b : bb) {
      // block mesh
      auto& m = mk.at(b)->GetMesh();
      // block cells
      auto& mbc = m.GetBlockCells();
      // resize field for block mesh
      fc.Reinit(m);
      // load from grid to lab
      s_.l->load(s_.mb[b], stage_);
      // read from lab to fc
      ReadBuffer(fc, *s_.l, e, m);
      // get corner of inner cells block
      MIdx wb = m.GetInBlockCells().GetBegin();
      // copy from inner cells to global field
      for (auto w : bc) {
        gfc[gbc.GetIdx(wb + w)] = fc[mbc.GetIdx(wb + w)];
      }
    }
    // recv from other ranks
    for (auto b : gbb) {
      if (!s_.mb.count(b)) { // not local block
        MPI_Status st;
        MPI_Recv(v.data(), v.size(), mt, MPI_ANY_SOURCE, 
                 MPI_ANY_TAG, comm_, &st);

        size_t i = 0;
        MIdx wb = gbb.GetMIdx(IdxCell(st.MPI_TAG)) * bs_;
        for (auto w : bc) {
          gfc[gbc.GetIdx(wb + w)] = v[i++];
        }
      }
    }

    MPI_Barrier(comm_);
    return gfc;
  } else {
    // send to root
    for (auto b : bb) {
      // block mesh
      auto& m = mk.at(b)->GetMesh();
      // block cells
      auto& mbc = m.GetBlockCells();
      // resize field for block mesh
      fc.Reinit(m);
      // load from grid to lab
      s_.l->load(s_.mb[b], stage_);
      // read from lab to fc
      ReadBuffer(fc, *s_.l, e, m);
      // get corner of inner cells block
      MIdx wb = m.GetInBlockCells().GetBegin();
      // copy from inner cells to v
      size_t i = 0;
      // copy from inner cells to global field
      for (auto w : bc) {
        v[i++] = fc[mbc.GetIdx(wb + w)];
      }
      // XXX: assume same order of Recv on root
      MPI_Send(v.data(), v.size(), mt, 0, gbb.GetIdx(b).GetRaw(), comm_);
    }
    MPI_Barrier(comm_);
    return FieldCell<Scal>();
  }
}

template <class B, int i>
const std::string StreamHdf<B, i>::NAME = "alpha";
template <class B, int i>
const std::string StreamHdf<B, i>::EXT = "00";

template <class B>
std::string StreamHdfDyn<B>::NAME = "alpha";
template <class B>
const std::string StreamHdfDyn<B>::EXT = "";
template <class B>
int StreamHdfDyn<B>::ID = 0;

} // namespace cubism_impl

template <class Scal_, size_t bx_, size_t by_, size_t bz_, size_t es_>
using GPar = cubism_impl::GPar<Scal_, bx_, by_, bz_, es_>;

template <class Par, class KF>
using Cubism = cubism_impl::Cubism<Par, KF>;
