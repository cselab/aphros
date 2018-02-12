#include "BlockInfo.h"
#include "Grid.h"
#include "GridMPI.h"
#include "BlockLab.h"
#include "BlockLabMPI.h"
#include "StencilInfo.h"
#include "HDF5Dumper_MPI.h"
#include "ICubism.h"

#define NUMFR 10
#define TEND 100

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

using Scal = double;

class Cubism {
 public:
  using Idx = std::array<int, 3>;

  Cubism(MPI_Comm comm, KernelFactory<K>& kf, 
      int bs, Idx b, Idx p, int es, int h);

  bool IsDone() const;
  void Step();

 private:
  std::map<Idx, std::unique_ptr<K>> mk;

  static Idx GetIdx(const int* d) {
    return {d[0], d[1], d[2]};
  }

  int bs_; // block size
  int es_; // element size in Scal
  int h_; // number of halo cells (same in all directions)

  TGrid g_;

  int step_ = 0;
  int stage_ = 0;
  int frame_ = 0;

  static StencilInfo GetStencil(int h) {
    return StencilInfo(-h,-h,-h,h+1,h+1,h+1, true, 8, 0,1,2,3,4,5,6,7);
  }

  bool isroot_;

  void ReadBuffer(Hydro<MeshStructured>&);
};

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

struct FakeProc {
  StencilInfo stencil;
  explicit FakeProc(StencilInfo si) 
    : stencil(si)
  {}
};

