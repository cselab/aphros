#pragma once

#include <memory>

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

using Real = double;

struct Elem {
  static const size_t s = 8;
  Real a[s];
  void init(Real val) { 
    for (size_t i = 0; i < s; ++i) {
      a[i] = val;
    }
  }
  void clear() {
    init(0);
  }

  Elem& operator=(const Elem&) = default;
};

struct Block {
  static const int bs = _BLOCKSIZE_;
  static const int sx = bs;
  static const int sy = bs;
  static const int sz = bs;
  static const int n = sx * sy * sz;

  // required by framework
  static const int sizeX = sx;
  static const int sizeY = sy;
  static const int sizeZ = sz;


  // floats per element
  static const int fe = sizeof(Elem) / sizeof(Real);

  using ElementType = Elem;
  using element_type = Elem;

  Elem __attribute__((__aligned__(_ALIGNBYTES_))) data[bs][bs][bs];

  Real __attribute__((__aligned__(_ALIGNBYTES_))) tmp[bs][bs][bs][fe];

  void clear_data() {
    Elem* e = &data[0][0][0];
    for(int i = 0; i < n; ++i) {
      e[i].clear();
    }
  }

  void clear_tmp() {
    Real* t = &tmp[0][0][0][0];
    for(int i = 0; i < n * fe; ++i) {
      t[i] = 0;
    }
  }

  void clear() {
    clear_data();
    clear_tmp();
  }

  inline Elem& operator()(int ix, int iy=0, int iz=0) {
    assert(ix>=0 && ix<sx);
    assert(iy>=0 && iy<sy);
    assert(iz>=0 && iz<sz);

    return data[iz][iy][ix];
  }
};


typedef Block Block_t;  
typedef Grid<Block_t, std::allocator> GridBase;
typedef GridBase Grid_t;

template<typename BlockType, template<typename X> class Alloc=std::allocator>
class LabPer: public BlockLab<BlockType,Alloc>
{
  typedef typename BlockType::ElementType ElementTypeBlock;

 public:
  virtual inline std::string name() const { return "name"; }
  bool is_xperiodic() {return true;}
  bool is_yperiodic() {return true;}
  bool is_zperiodic() {return true;}

  LabPer()
    : BlockLab<BlockType,Alloc>(){}
};



using Lab = LabPer<Block_t, std::allocator>;

typedef BlockLabMPI<Lab> LabMPI;
typedef GridMPI<Grid_t> GridMPI_t;
using TGrid = GridMPI_t;

using Scal = double;

template <class KF>
class Cubism : public DistrMesh<KF> {
 public:
  Cubism(MPI_Comm comm, KF& kf, int bs, int es, int h, Vars& par);

 private:
  using K = typename KF::K;
  using M = typename KF::M;
  using P = DistrMesh<KF>;
  using MIdx = typename M::MIdx;

  using P::mk;
  using P::kf_;
  using P::par;
  using P::bs_;
  using P::es_;
  using P::hl_;
  using P::p_;
  using P:: b_; 
  using P::step_;
  using P::stage_;
  using P::frame_;
  using P::isroot_;
  using P::comm_;

  TGrid g_;
  struct S { // cubism [s]tate
    SynchronizerMPI* s;
    std::unique_ptr<LabMPI> l;
    std::map<MIdx, BlockInfo, typename MIdx::LexLess> mb;
  };
  S s_;

  static StencilInfo GetStencil(int h) {
    return StencilInfo(-h,-h,-h,h+1,h+1,h+1, true, 8, 0,1,2,3,4,5,6,7);
  }
  static std::vector<MyBlockInfo> GetBlocks(
      const std::vector<BlockInfo>&, size_t hl);

  void ReadBuffer(M& m, LabMPI& l) {
    using MIdx = typename M::MIdx;
    int e = 0; // buffer field idx

    for (auto u : m.GetComm()) {
      for (auto i : m.AllCells()) {
        auto& bc = m.GetBlockCells();
        auto d = bc.GetMIdx(i) - MIdx(hl_) - bc.GetBegin();
        (*u)[i] = l(d[0], d[1], d[2]).a[e];
      }
      ++e;
    }

    m.ClearComm();
  }

  void WriteBuffer(M& m, Block_t& o) {
    using MIdx = typename M::MIdx;
    int bs = _BLOCKSIZE_;

    // Check buffer has enough space for all fields
    assert(m.GetComm().size() <= Elem::s && "Too many fields for Comm()");

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
  void Dump(int frame, int step) override;
};


template <int ID>
struct StreamHdf {
  static const std::string NAME;
  static const std::string EXT;
  static const int NCHANNELS = 1;
  static const int CLASS = 0;
  struct T { Real a[8]; };

  using B = Block;
  B& b;

  StreamHdf(B& b): b(b) {}

  // write
  void operate(const int ix, const int iy, const int iz, Real out[0]) const
  {
    const T& in = *((T*)&b.data[iz][iy][ix].a[0]);
    out[0] = in.a[ID];
  }

  // read
  void operate(const Real out[0], const int ix, const int iy, const int iz) const
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

template <class KF>
std::vector<MyBlockInfo> Cubism<KF>::GetBlocks(
    const std::vector<BlockInfo>& cc, size_t hl) {
  std::vector<MyBlockInfo> bb;
  for(size_t i = 0; i < cc.size(); i++) {
    const BlockInfo& c = cc[i];
    MyBlockInfo b;
    for (int j = 0; j < 3; ++j) {
      b.index[j] = c.index[j];
      b.origin[j] = c.origin[j];
    }
    b.h_gridpoint = c.h_gridpoint;
    b.ptrBlock = c.ptrBlock;
    b.hl = hl;
    bb.push_back(b);
  }
  return bb;
}


template <class KF>
Cubism<KF>::Cubism(MPI_Comm comm, KF& kf, 
    int bs, int es, int hl, Vars& par) 
  : DistrMesh<KF>(comm, kf, bs, es, hl, par)
  , g_(p_[0], p_[1], p_[2], b_[0], b_[1], b_[2], 1., comm)
{
  std::vector<BlockInfo> cc = g_.getBlocksInfo(); // [c]ubism block info
  std::vector<MyBlockInfo> ee = GetBlocks(cc, hl_);

  for (auto& e : ee) {
    auto d = e.index;
    // TODO: constructor
    MIdx b(d[0], d[1], d[2]);
    mk.emplace(b, std::unique_ptr<K>(kf_.Make(par, e)));
  }

  int r;
  MPI_Comm_rank(comm, &r);
  isroot_ = (0 == r);

  comm_ = g_.getCartComm(); // XXX: overwrite comm_
}

template <class KF>
auto Cubism<KF>::GetBlocks() -> std::vector<MIdx> {
  MPI_Barrier(comm_);

  // 1. Exchange halos in buffer mesh
  FakeProc fp(GetStencil(hl_));       // object with field 'stencil'
  SynchronizerMPI& s = g_.sync(fp); 

  s_.l.reset(new LabMPI);
  s_.l->prepare(g_, s);   // allocate memory for lab cache

  MPI_Barrier(comm_);

  // Do communication and get all blocks
  std::vector<BlockInfo> aa = s.avail();

  // Put blocks to map by index 
  std::vector<MIdx> bb;
  s_.mb.clear();
  for (auto a : aa) {
    auto d = a.index;
    MIdx b(d[0], d[1], d[2]);
    s_.mb.emplace(b, a);
    bb.push_back(b);
  }

  return bb;
}

template <class KF>
void Cubism<KF>::ReadBuffer(const std::vector<MIdx>& bb) {
  for (auto& b : bb) {
    auto& k = *mk.at(b); // kernel
    auto& m = k.GetMesh();
    s_.l->load(s_.mb[b], stage_);
    ReadBuffer(m, *s_.l);
  }
}

template <class KF>
void Cubism<KF>::WriteBuffer(const std::vector<MIdx>& bb) {
  for (auto& b : bb) {
    auto& k = *mk.at(b); // kernel
    auto& m = k.GetMesh();
    WriteBuffer(m, *(typename TGrid::BlockType*)s_.mb[b].ptrBlock);
  }
}

template <class KF>
void Cubism<KF>::Reduce(const std::vector<MIdx>& bb) {
  auto& f = *mk.at(bb[0]); // first kernel
  auto& mf = f.GetMesh();
  auto& vf = mf.GetReduce();  // pointers to reduce

  std::vector<Scal> r(vf.size(), 0); // results

  // Check size is the same for all kernels
  for (auto& b : bb) {
    auto& k = *mk.at(b); // kernel
    auto& m = k.GetMesh();
    auto& v = m.GetReduce();  // pointers to reduce
    assert(v.size() == r.size());
  }

  // Reduce over all kernels on current rank
  for (auto& b : bb) {
    auto& k = *mk.at(b); 
    auto& m = k.GetMesh();
    auto& v = m.GetReduce();  
    for (size_t i = 0; i < r.size(); ++i) {
      r[i] += *v[i];
    }
  }

  // Reduce over all ranks
  MPI_Allreduce(
      MPI_IN_PLACE, r.data(), r.size(), 
      MPI_DOUBLE, MPI_SUM, comm_); // TODO: type from Scal

  // Write results to all kernels on current rank
  for (auto& b : bb) {
    auto& k = *mk.at(b); 
    auto& m = k.GetMesh();
    auto& v = m.GetReduce();  
    for (size_t i = 0; i < r.size(); ++i) {
      *v[i] = r[i];
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
void Cubism<KF>::Dump(int frame, int step) {
  auto suff = "_" + std::to_string(frame_);
  DumpHDF5_MPI<TGrid, StreamHdf<0>>(g_, frame, step*1., "p" + suff);
}


template <int i>
const std::string StreamHdf<i>::NAME = "alpha";
template <int i>
const std::string StreamHdf<i>::EXT = "00";
