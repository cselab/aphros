#include "vect.hpp"

template <class _Scal, size_t _dim>
class IMesh {
 public:
  using Scal = _Scal;
  static constexpr size_t dim = _dim;
  using Vect = geom::Vect<Scal, dim>;
  MeshGeneric() = default;
  virtual ~MeshGeneric() {}
  virtual Vect GetCenter(IdxCell) const = 0;
  virtual Vect GetCenter(IdxFace) const = 0;
  virtual Vect GetSurface(IdxFace) const = 0;
  virtual Vect GetNode(IdxNode) const = 0;
  virtual Scal GetVolume(IdxCell) const = 0;
  virtual Scal GetArea(IdxFace) const = 0;
  virtual size_t GetNumCells() const = 0;
  virtual size_t GetNumFaces() const = 0;
  virtual size_t GetNumNodes() const = 0;
  virtual IdxCell GetNeighbourCell(IdxCell, size_t) const = 0;
  virtual IdxFace GetNeighbourFace(IdxCell, size_t) const = 0;
  virtual Scal GetOutwardFactor(IdxCell, size_t) const = 0;
  virtual Vect GetOutwardSurface(IdxCell, size_t) const = 0;
  virtual IdxNode GetNeighbourNode(IdxCell, size_t) const = 0;
  virtual IdxCell GetNeighbourCell(IdxFace, size_t) const = 0;
  virtual IdxNode GetNeighbourNode(IdxFace, size_t) const = 0;
  // For any cell, indices of its neighbour cells and faces are the same.
  // Thus, there's no GetNumNeighbourCells(IdxCell) function.
  virtual size_t GetNumNeighbourFaces(IdxCell) const = 0;
  virtual size_t GetNumNeighbourNodes(IdxCell) const = 0;
  virtual size_t GetNumNeighbourCells(IdxFace) const = 0;
  virtual size_t GetNumNeighbourNodes(IdxFace) const = 0;
  operator Range<IdxCell>() const {
    return Range<IdxCell>(0, GetNumCells());
  }
  operator Range<IdxFace>() const {
    return Range<IdxFace>(0, GetNumFaces());
  }
  operator Range<IdxNode>() const {
    return Range<IdxNode>(0, GetNumNodes());
  }
  Range<IdxCell> AllCells() const {
    return Range<IdxCell>(*this);
  }
  Range<IdxFace> AllFaces() const {
    return Range<IdxFace>(*this);
  }
  Range<IdxNode> AllNodes() const {
    return Range<IdxNode>(*this);
  }
  virtual bool IsInner(IdxCell) const = 0;
  virtual bool IsInner(IdxFace) const = 0;
  virtual bool IsInside(IdxCell, Vect) const = 0;
  virtual IdxCell FindNearestCell(Vect point) const {
    IdxCell idxcell_nearest = IdxCell(0);
    for (auto idxcell : AllCells()) {
      if (point.dist(GetCenter(idxcell)) <
          point.dist(GetCenter(idxcell_nearest))) {
        idxcell_nearest = idxcell;
      }
    }
    return idxcell_nearest;
  }
  virtual size_t GetValidNeighbourCellId(IdxFace idxface) const {
    size_t id = 0;
    while (id < GetNumNeighbourCells(idxface) &&
        GetNeighbourCell(idxface, id).IsNone()) {
      ++id;
    }
    return id;
  }
  virtual Vect GetNormal(IdxFace idxface) const {
    return GetSurface(idxface) / GetArea(idxface);
  }
};


template <class IdxType, size_t dim>
class BlockGeneric {
  using MIdx = geom::MIdxGeneral<dim>;
  MIdx begin_, size_, end_;

 public:
  class iterator {
    const BlockGeneric* owner_;
    MIdx midx_;

   public:
    explicit iterator(const BlockGeneric* owner, MIdx midx)
        : owner_(owner)
        , midx_(midx)
    {}
    iterator& operator++() {
      for (size_t n = 0; n < dim; ++n) {
        ++midx_[n];
        if (midx_[n] == owner_->end_[n] && n < dim - 1) {
          midx_[n] = owner_->begin_[n];
        } else {
          break;
        }
      }
      return *this;
    }
    iterator& operator--() {
      for (size_t n = 0; n < dim; ++n) {
        if (midx_[n] == owner_->begin_[n]) {
          midx_[n] = owner_->end_[n] - 1;
        } else {
          --midx_[n];
          break;
        }
      }
      return *this;
    }
    bool operator==(const iterator& other) const {
      return midx_ == other.midx_;
    }
    bool operator!=(const iterator& other) const {
      return !(*this == other);
    }
    MIdx operator*() const {
      return midx_;
    }
  };

  BlockGeneric()
      : begin_(MIdx::kZero), size_(MIdx::kZero), end_(MIdx::kZero)
  {}
  BlockGeneric(MIdx size)
      : begin_(MIdx::kZero), size_(size), end_(begin_ + size_)
  {}
  BlockGeneric(MIdx begin, MIdx size)
      : begin_(begin), size_(size), end_(begin_ + size_)
  {}
  MIdx GetDimensions() const {
    return size_;
  }
  MIdx GetBegin() const {
    return begin_;
  }
  MIdx GetEnd() const {
    return end_;
  }
  size_t size() const {
    size_t res = 1;
    for (size_t i = 0; i < dim; ++i) {
      res *= size_[i];
    }
    return res;
  }
  operator Range<IdxType>() const {
    return Range<IdxType>(0, size());
  }
/*  void Advance(IdxType& idx, Direction dir, int n) const {
    switch (dir) {
      case Direction::i: {
        idx.AddRaw(n);
        break;
      }
      case Direction::j: {
        idx.AddRaw(n*size_[0]);
        break;
      }
      case Direction::k: {
        idx.AddRaw(n*size_[0]*size_[1]);
        break;
      }
    }
  }
  void Advance(IdxType& idx, MIdx midx) const {
    // TODO: unwrap multiple Advance()
    Advance(idx, Direction::i, midx[0]);
    Advance(idx, Direction::j, midx[1]);
    Advance(idx, Direction::k, midx[2]);
  }
  IdxType GetAdvanced(IdxType idx, MIdx midx) const {
    Advance(idx, midx);
    return idx;
  }
  void Inc(IdxType& idx, Direction dir) const {
    Advance(idx, dir, 1);
  }
  IdxType GetInc(IdxType idx, Direction dir) const {
    Inc(idx, dir);
    return idx;
  }
  void Dec(IdxType& idx, Direction dir) const {
    Advance(idx, dir, -1);
  }
  IdxType GetDec(IdxType idx, Direction dir) const {
    Dec(idx, dir);
    return idx;
  }*/
  IdxType GetIdx(MIdx midx) const {
    midx -= begin_;
    size_t res = 0;
    for (size_t i = dim; i != 0; ) {
      --i;
      res *= size_[i];
      res += midx[i];
    }
    return IdxType(res);
  }
  MIdx GetMIdx(IdxType idx) const {
    MIdx midx;
    size_t raw = idx.GetRaw();
    for (size_t i = 0; i < dim; ++i) {
      midx[i] = raw % size_[i];
      raw /= size_[i];
    }
    return begin_ + midx;
  }
  bool IsInside(MIdx midx) const {
    return begin_ <= midx && midx < end_;
  }
  // TODO: rename to clip
  void Cut(MIdx& midx) const {
    for (size_t i = 0; i < midx.size(); ++i) {
      midx[i] = std::max(
          begin_[i],
          std::min(
              end_[i] - 1,
              midx[i]));
    }
  }
  iterator begin() const {
    return iterator(this, begin_);
  }
  iterator end() const {
    MIdx midx = begin_;
    midx[dim - 1] = end_[dim - 1];
    return iterator(this, midx);
  }
};

template <size_t dim>
using BlockCells = BlockGeneric<IdxCell, dim>;

template <size_t dim>
using BlockNodes = BlockGeneric<IdxNode, dim>;


template <size_t dim>
class Direction {
  using MIdx = Vect<IntIdx, dim>;
  size_t diridx;
 public:
  Direction() {}
  explicit Direction(size_t i)
    : diridx(i) {}
  char GetLetter() const {
    return std::string("xyz")[diridx];
  }
  operator size_t() const {
    return diridx;
  }
  operator MIdx() const {
    MIdx res = MIdx::kZero;
    ++res[diridx];
    return res;
  }
  bool operator==(const Direction& other) const {
    return diridx == other.diridx;
  }
  bool operator!=(const Direction& other) const {
    return !((*this) == other);
  }
  bool operator<(const Direction& other) const {
    return diridx < other.diridx;
  }
  static const Direction i;
  static const Direction j;
  static const Direction k;
};

template <size_t dim>
const Direction<dim> Direction<dim>::i(0);
template <size_t dim>
const Direction<dim> Direction<dim>::j(1);
template <size_t dim>
const Direction<dim> Direction<dim>::k(2);

template <size_t dim>
class BlockFaces {
  using MIdx = MIdxGeneral<dim>;
  using Direction = geom::Direction<dim>;
  MIdx begin_;
  MIdx block_cells_size_;
  size_t GetNumFaces(Direction dir) const {
    MIdx bcs = block_cells_size_;
    ++bcs[dir];
    size_t res = 1;
    for (size_t i = 0; i < dim; ++i) {
      res *= bcs[i];
    }
    return res;
  }
  size_t GetNumFaces(size_t i) const {
    return GetNumFaces(Direction(i));
  }
  size_t GetFlat(MIdx midx, Direction dir) const {
    midx -= begin_;
    MIdx bcs = block_cells_size_;
    ++bcs[dir];
    size_t res = 0;
    for (size_t i = dim; i != 0; ) {
      --i;
      res *= bcs[i];
      res += midx[i];
    }
    return res;
  }
  MIdx GetMIdxFromOffset(size_t raw, Direction dir) const {
    MIdx bcs = block_cells_size_;
    ++bcs[dir];
    MIdx midx;
    for (size_t i = 0; i < dim; ++i) {
      midx[i] = raw % bcs[i];
      raw /= bcs[i];
    }
    return begin_ + midx;
  }
  std::pair<MIdx, Direction> GetMIdxDirection(IdxFace idxface) const {
    size_t raw = idxface.GetRaw();
    size_t diridx = 0;
    while (raw >= GetNumFaces(diridx) && diridx < dim) {
      raw -= GetNumFaces(diridx);
      ++diridx;
    }
    Direction dir(diridx);
    return {GetMIdxFromOffset(raw, dir), dir};
  }

 public:
  BlockFaces()
      : begin_(MIdx::kZero), block_cells_size_(MIdx::kZero)
  {}
  BlockFaces(MIdx block_cells_size)
      : begin_(MIdx::kZero), block_cells_size_(block_cells_size)
  {}
  BlockFaces(MIdx begin, MIdx block_cells_size)
      : begin_(begin), block_cells_size_(block_cells_size)
  {}
  size_t size() const {
    size_t res = 0;
    for (size_t i = 0; i < dim; ++i) {
      res += GetNumFaces(i);
    }
    return res;
  }
  operator Range<IdxFace>() const {
    return Range<IdxFace>(0, size());
  }
  IdxFace GetIdx(MIdx midx, Direction dir) const {
    size_t raw = 0;
    for (size_t i = 0; i < dir; ++i) {
      raw += GetNumFaces(i);
    }
    raw += GetFlat(midx, dir);
    return IdxFace(raw);
  }
  MIdx GetMIdx(IdxFace idx) const {
    return GetMIdxDirection(idx).first;
  }
  Direction GetDirection(IdxFace idx) const {
    return GetMIdxDirection(idx).second;
  }
};

