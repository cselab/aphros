#pragma once

#include <memory>
#include <limits>
#include <map>

#include "Cubism/BlockInfo.h"
#include "Cubism/Grid.h"
#include "Cubism/GridMPI.h"
#include "Cubism/BlockLab.h"
#include "Cubism/BlockLabMPI.h"
#include "Cubism/StencilInfo.h"
#include "Cubism/HDF5Dumper_MPI.h"
#include "ICubism.h"
#include "Vars.h"
#include "Distr.h"

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
  Cubism(MPI_Comm comm, KF& kf, Vars& par);

 private:
  using Lab = GLab<Par>;
  using Block = GBlock<Par>;
  using Grid = GGrid<Block>;
  using Elem = typename Block::Elem;
  using Synch = typename Grid::Synch;

  using K = typename KF::K;
  using M = typename KF::M;
  using P = DistrMesh<KF>; // parent
  using MIdx = typename M::MIdx;
  using Scal = typename M::Scal;

  using P::mk;
  using P::kf_;
  using P::par;
  using P::bs_;
  using P::es_;
  using P::hl_;
  using P::p_;
  using P::b_; 
  using P::stage_;
  using P::frame_;
  using P::isroot_;
  using P::comm_;
  using P::ext_;

  Grid g_;
  struct S { // cubism [s]tate
    Synch* s;
    std::unique_ptr<Lab> l;
    std::map<MIdx, BlockInfo, typename MIdx::LexLess> mb;
  };
  S s_;

  static StencilInfo GetStencil(int hl /*halo cells*/,
                                int cs /*comm size*/) {
    const int a = -hl;
    const int b = hl + 1;
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
  static std::vector<MyBlockInfo> GetBlocks(
      const std::vector<BlockInfo>&, MIdx bs, size_t hl);

  void ReadBuffer(M& m, Lab& l) {
    using MIdx = typename M::MIdx;

    // Check buffer has enough space for all fields
    assert(m.GetComm().size() <= Elem::es && "Too many fields for Comm()");

    int e = 0; // buffer field idx

    for (auto u : m.GetComm()) {
      for (auto i : m.AllCells()) {
        auto& bc = m.GetBlockCells();
        auto d = bc.GetMIdx(i) - MIdx(hl_) - bc.GetBegin();
        (*u)[i] = l(d[0], d[1], d[2]).a[e];
      }
      ++e;
    }
  }

  void WriteBuffer(M& m, Block& o) {
    using MIdx = typename M::MIdx;

    // Check buffer has enough space for all fields
    assert(m.GetComm().size() <= Elem::es && "Too many fields for Comm()");

    int e = 0; // buffer field idx

    for (auto u : m.GetComm()) {
      for (auto i : m.Cells()) {
        auto& bc = m.GetBlockCells();
        auto d = m.GetBlockCells().GetMIdx(i) - MIdx(hl_) - bc.GetBegin();
        o.data[d[2]][d[1]][d[0]].a[e] = (*u)[i];
      }
      ++e;
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

  std::vector<BlockInfo> cc = g_.getBlocksInfo(); // [c]ubism block info
  std::vector<MyBlockInfo> ee = GetBlocks(cc, bs_, hl_);

  for (auto& e : ee) {
    MIdx b(e.index);
    mk.emplace(b, std::unique_ptr<K>(kf_.Make(par, e)));
  }

  int r;
  MPI_Comm_rank(comm, &r);
  isroot_ = (0 == r);

  comm_ = g_.getCartComm(); // XXX: overwrite comm_
}

template <class Par, class KF>
auto Cubism<Par, KF>::GetBlocks() -> std::vector<MIdx> {
  MPI_Barrier(comm_);

  // Get all blocks
  std::vector<BlockInfo> cc = g_.getBlocksInfo(); // all blocks
  auto& m = mk.at(MIdx(cc[0].index))->GetMesh();
  size_t cs = m.GetComm().size(); // comm size

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

  for (size_t i = 0; i < vf.size(); ++i) {
    Scal r; // result
    std::string s = vf[i].second; // operation string

    enum class Op { sum, prod, max, min };
    Op o;
    if (s == "sum") {
      o = Op::sum;
    } else if (s == "prod") {
      o = Op::prod;
    } else if (s == "max") {
      o = Op::max;
    } else if (s == "min") {
      o = Op::min;
    } else {
      std::cerr << "Reduce(): unknown operation '" << s <<  "'" << std::endl;
      assert(false);
    }

    
    // Init result
    switch (o) {
      case Op::sum:  
        r = 0;  
        break;
      case Op::prod: 
        r = 1;  
        break;
      case Op::max:  
        r = -std::numeric_limits<double>::max();  
        break;
      case Op::min:  
        r = std::numeric_limits<double>::max();  
        break;
      default:
        assert(false);
    }

    // Reduce over all blocks on current rank
    for (auto& b : bb) {
      auto& v = mk.at(b)->GetMesh().GetReduce(); 
      Scal a = *v[i].first;
      switch (o) {
        case Op::sum:  
          r += a;  
          break;
        case Op::prod: 
          r *= a;  
          break;
        case Op::max:  
          r = std::max(r, a);
          break;
        case Op::min:  
          r = std::min(r, a);
          break;
        default:
          assert(false);
      }
    }

    // Reduce over all ranks
    std::map<Op, MPI_Op> mo = {
        {Op::sum, MPI_SUM},
        {Op::prod, MPI_PROD},
        {Op::max, MPI_MAX},
        {Op::min, MPI_MIN}
      };
    std::map<size_t, MPI_Datatype> mt = {
        {4, MPI_FLOAT},
        {8, MPI_DOUBLE}
      };
    MPI_Allreduce(MPI_IN_PLACE, &r, 1, mt.at(sizeof(Scal)), mo.at(o), comm_); 

    // Write results to all blocks on current rank
    for (auto& b : bb) {
      auto& v = mk.at(b)->GetMesh().GetReduce(); 
      *v[i].first = r;
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
    size_t k = m.GetComm().size() - m.GetDump().size();
    for (auto d : m.GetDump()) {
      auto suff = "_" + std::to_string(frame_);
      StreamHdfDyn<Block>::ID = k;
      StreamHdfDyn<Block>::NAME = d.second;
      DumpHDF5_MPI<Grid, StreamHdfDyn<Block>>(
          g_, frame_, frame_, d.second + suff);
      ++k;
    }
    ++frame_;
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
