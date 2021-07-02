// Created by Petr Karnakov on 01.02.2021
// Copyright 2021 ETH Zurich

#include <cstdint>
#include <fstream>
#include <limits>
#include <map>
#include <sstream>

#include "util/distr.h"
#include "util/format.h"
#include "util/logger.h"

#include "raw.h"

namespace dump {

#if USEFLAG(MPI)
template <class Type>
static MPI_Datatype GetMpiType(Type type) {
  switch (type) {
    case Type::UInt16:
      return MPI_UINT16_T;
    case Type::Float32:
      return MPI_FLOAT;
    case Type::Float64:
      return MPI_DOUBLE;
    default:
      fassert(false);
  }
}

template <class T>
MPI_Datatype GetMpiType() {
  return MPI_DATATYPE_NULL;
}

template <>
MPI_Datatype GetMpiType<double>() {
  return MPI_DOUBLE;
}

template <>
MPI_Datatype GetMpiType<float>() {
  return MPI_FLOAT;
}

template <>
MPI_Datatype GetMpiType<std::uint16_t>() {
  return MPI_UINT16_T;
}

template <class MIdx>
static void CreateSubarray(
    MIdx wsize, MIdx wsubsize, MIdx wstart, MPI_Datatype numbertype,
    /*out*/ MPI_Datatype* filetype) {
  const auto dim = MIdx::dim;
  int size[dim];
  int subsize[dim];
  int start[dim];
  for (size_t d = 0; d < dim; ++d) {
    size[d] = wsize[dim - d - 1];
    subsize[d] = wsubsize[dim - d - 1];
    start[d] = wstart[dim - d - 1];
  }
  MPI_Type_create_subarray(
      dim, size, subsize, start, MPI_ORDER_C, numbertype, filetype);
  MPI_Type_commit(filetype);
}
#endif

template <class M>
template <class T>
void Raw<M>::Write(
    const std::string& path, const std::vector<MIdx>& starts,
    const std::vector<MIdx>& sizes, const std::vector<std::vector<T>>& data,
    const MIdx global_size, Type type, const MpiWrapper& mpi) {
  const auto nblocks = data.size();
  fassert(nblocks > 0);
  MIdx comb_start(std::numeric_limits<IntIdx>::max());
  MIdx comb_end(std::numeric_limits<IntIdx>::min());
  for (size_t b = 0; b < nblocks; ++b) {
    comb_start = comb_start.min(starts[b]);
    comb_end = comb_start.max(starts[b] + sizes[b]);
  }
  // size of block combined from all blocks on current rank
  const MIdx comb_size = comb_end - comb_start;

  // XXX check the blocks on current rank constitute a rectangular subarray
  fassert_equal(comb_size % sizes[0], MIdx(0));
  fassert_equal((comb_size / sizes[0]).prod(), (int)nblocks);

  typename M::IndexCells comb_indexer(comb_start, comb_size);

  auto to_elemtype = [type](T value, void* ptr) {
    switch (type) {
      case Type::UInt16:
        (*reinterpret_cast<std::uint16_t*>(ptr)) = value;
        break;
      case Type::Float32:
        (*reinterpret_cast<float*>(ptr)) = value;
        break;
      case Type::Float64:
        (*reinterpret_cast<double*>(ptr)) = value;
        break;
      default:
        fassert(false);
    }
  };
  const size_t elemtypesize = GetPrecision(type);
  std::vector<char> buf(comb_size.prod() * elemtypesize);
  for (size_t b = 0; b < nblocks; ++b) {
    const typename M::BlockCells block(starts[b], sizes[b]);
    size_t i = 0;
    for (auto c : GRangeIn<IdxCell, dim>(comb_indexer, block)) {
      to_elemtype(data[b][i++], &buf[c.raw() * elemtypesize]);
    }
  }

#if USEFLAG(MPI)
  const MPI_Datatype elemtype = GetMpiType(type);
  MPI_Datatype filetype;
  CreateSubarray(global_size, comb_size, comb_start, elemtype, &filetype);

  MPI_Datatype memtype;
  CreateSubarray(comb_size, comb_size, MIdx(0), elemtype, &memtype);

  const MPI_Comm comm = mpi.GetComm();
  MPI_File fh;
  try {
    MPICALL(MPI_File_open(
        comm, path.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL,
        &fh));
  } catch (const std::exception& e) {
    fassert(
        false,
        std::string() + e.what() + ", while opening file '" + path + "'");
  }
  MPICALL(
      MPI_File_set_view(fh, 0, elemtype, filetype, "native", MPI_INFO_NULL));

  MPICALL(MPI_File_write_all(fh, buf.data(), 1, memtype, MPI_STATUS_IGNORE));
  MPI_File_close(&fh);
  MPI_Type_free(&memtype);
  MPI_Type_free(&filetype);
#else
  (void)mpi;
  (void)global_size;
  std::ofstream file(path, std::ios::binary);
  fassert(file.good(), "Can't open file '" + path + "' for write");
  file.write(buf.data(), buf.size());
#endif
}

template <class M>
template <class T>
void Raw<M>::Write(
    const FieldCell<T>& fc, const Meta& meta, std::string path, M& m) {
  auto sem = m.GetSem();
  struct {
    std::vector<MIdx> starts;
    std::vector<MIdx> sizes;
    std::vector<std::vector<T>> data;
  } * ctx(sem);
  auto& t = *ctx;

  if (sem("gather")) {
    const auto bc = m.GetInBlockCells();
    t.starts.push_back(bc.GetBegin());
    t.sizes.push_back(bc.GetSize());

    t.data.emplace_back();
    auto& data = t.data.back();
    data.reserve(bc.size());
    for (auto c : m.Cells()) {
      data.push_back(fc[c]);
    }

    m.GatherToLead(&t.starts);
    m.GatherToLead(&t.sizes);
    m.GatherToLead(&t.data);
  }
  if (sem("write") && m.IsLead()) {
    MpiWrapper mpi(m.GetMpiComm());
    Write(path, t.starts, t.sizes, t.data, m.GetGlobalSize(), meta.type, mpi);
  }
  if (sem()) { // XXX empty stage
  }
}

template <class M>
template <class T>
void Raw<M>::Read(FieldCell<T>& fc, const Meta& meta, std::string path, M& m) {
  auto sem = m.GetSem();
  struct {
    std::vector<MIdx> starts;
    std::vector<MIdx> sizes;
    std::vector<T> data;
    std::vector<std::vector<T>*> dataptr;
  } * ctx(sem);
  auto& t = *ctx;

  if (sem("gather")) {
    const auto bc = m.GetInBlockCells();
    t.starts.push_back(bc.GetBegin());
    t.sizes.push_back(bc.GetSize());
    t.data.resize(bc.size());
    t.dataptr.push_back(&t.data);

    m.GatherToLead(&t.starts);
    m.GatherToLead(&t.sizes);
    m.GatherToLead(&t.dataptr);
  }
  if (sem("read") && m.IsLead()) {
    const auto nblocks = t.dataptr.size();
    fassert(nblocks > 0);
    MIdx comb_start(std::numeric_limits<IntIdx>::max());
    MIdx comb_end(std::numeric_limits<IntIdx>::min());
    for (size_t b = 0; b < nblocks; ++b) {
      comb_start = comb_start.min(t.starts[b]);
      comb_end = comb_start.max(t.starts[b] + t.sizes[b]);
    }
    // size of region combined from all blocks on current rank
    const MIdx comb_size = comb_end - comb_start;

    // XXX check the blocks on current rank constitute a rectangular subarray
    fassert_equal(comb_size % t.sizes[0], MIdx(0));
    fassert_equal((comb_size / t.sizes[0]).prod(), (int)nblocks);
    // check the hyperslab has the same dimensions as the global mesh
    fassert_equal(meta.count, m.GetGlobalSize());

    typename M::IndexCells comb_indexer(comb_start, comb_size);

    fassert_equal(meta.stride, MIdx(1), ". Unsupported stride");
    fassert(
        MIdx(0) <= meta.start,
        util::Format("Negative hyperslab start={}", meta.start));
    fassert(
        meta.start + meta.count * meta.stride <= meta.dimensions,
        util::Format(
            "Hyperslab start+count*stride={} beyond array dimensions={}",
            meta.start + meta.count * meta.stride, meta.dimensions));

    const size_t elemtypesize = GetPrecision(meta.type);
    std::vector<char> buf(comb_size.prod() * elemtypesize);
#if USEFLAG(MPI)
    const MPI_Datatype elemtype = GetMpiType(meta.type);

    MPI_Datatype filetype;
    CreateSubarray(
        meta.dimensions, comb_size, meta.start + comb_start, elemtype,
        &filetype);

    MPI_Datatype memtype;
    CreateSubarray(comb_size, comb_size, MIdx(0), elemtype, &memtype);

    const MPI_Comm comm = m.GetMpiComm();
    MPI_File fh;
    try {
      MPICALL(MPI_File_open(
          comm, path.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh));
    } catch (const std::exception& e) {
      fassert(
          false,
          std::string() + e.what() + ", while opening file '" + path + "'");
    }

    MPICALL(MPI_File_set_view(
        fh, meta.seek, elemtype, filetype, "native", MPI_INFO_NULL));

    MPICALL(MPI_File_read_all(fh, buf.data(), 1, memtype, MPI_STATUS_IGNORE));
    MPI_File_close(&fh);
    MPI_Type_free(&memtype);
    MPI_Type_free(&filetype);
#else
    fassert_equal(
        meta.start, MIdx(0), ". Sequential reader only supports start=(0,0,0)");
    std::ifstream file(path, std::ios::binary);
    fassert(file.good(), "Can't open file '" + path + "' for reading");
    file.seekg(meta.seek);
    file.read(buf.data(), buf.size());
#endif

    auto from_elemtype = [type = meta.type](const void* ptr) -> T {
      switch (type) {
        case Type::UInt16:
          return *reinterpret_cast<const std::uint16_t*>(ptr);
        case Type::Float32:
          return *reinterpret_cast<const float*>(ptr);
        case Type::Float64:
          return *reinterpret_cast<const double*>(ptr);
        default:
          fassert(false);
      }
    };

    for (size_t b = 0; b < nblocks; ++b) {
      const typename M::BlockCells block(t.starts[b], t.sizes[b]);
      size_t i = 0;
      for (auto c : GRangeIn<IdxCell, dim>(comb_indexer, block)) {
        (*t.dataptr[b])[i++] = from_elemtype(&buf[c.raw() * elemtypesize]);
      }
    }
  }
  if (sem("copy")) {
    size_t i = 0;
    fc.Reinit(m);
    for (auto c : m.Cells()) {
      fc[c] = t.data[i++];
    }
  }
  if (sem()) { // XXX empty stage
  }
}

} // namespace dump
