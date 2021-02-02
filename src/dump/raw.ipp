// Created by Petr Karnakov on 01.02.2021
// Copyright 2021 ETH Zurich

#include <fstream>
#include <limits>
#include <sstream>

#include "raw.h"
#include "util/distr.h"
#include "util/format.h"
#include "util/logger.h"

namespace dump {

template <class MIdx>
void CreateSubarray(
    MIdx wsize, MIdx wsubsize, MIdx wstart, MPI_Datatype* filetype) {
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
      dim, size, subsize, start, MPI_ORDER_C, MPI_DOUBLE, filetype);
  MPI_Type_commit(filetype);
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
    const auto dim = M::dim;
    const auto nblocks = t.data.size();
    fassert(nblocks > 0);
    MIdx comb_start(std::numeric_limits<IntIdx>::max());
    MIdx comb_end(std::numeric_limits<IntIdx>::min());
    for (size_t b = 0; b < nblocks; ++b) {
      comb_start = comb_start.min(t.starts[b]);
      comb_end = comb_start.max(t.starts[b] + t.sizes[b]);
    }
    const MIdx comb_size = comb_end - comb_start;

    // XXX check the blocks on current rank constitute a rectangular subarray
    fassert_equal(comb_size % t.sizes[0], MIdx(0));
    fassert_equal((comb_size / t.sizes[0]).prod(), (int)nblocks);

    typename M::IndexCells comb_indexer(comb_start, comb_size);
    FieldCell<T> comb_field(comb_indexer);
    for (size_t b = 0; b < nblocks; ++b) {
      const typename M::BlockCells block(t.starts[b], t.sizes[b]);
      size_t i = 0;
      for (auto c : GRangeIn<IdxCell, dim>(comb_indexer, block)) {
        comb_field[c] = t.data[b][i++];
      }
    }

    MPI_Datatype filetype;
    CreateSubarray(m.GetGlobalSize(), comb_size, comb_start, &filetype);

    MPI_Datatype memtype;
    CreateSubarray(comb_size, comb_size, MIdx(0), &memtype);

    const MPI_Comm comm = m.GetMpiComm();
    MPI_File fh;
    try {
      MPICALL(MPI_File_open(
          comm, path.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL,
          &fh));
    } catch (const std::exception& e) {
      throw std::runtime_error(
          std::string() + e.what() + ", while opening file '" + path + "'");
    }
    // fassert(error == MPI_SUCCESS, "Can't open file '" + path + "'");
    MPICALL(MPI_File_set_view(
        fh, 0, MPI_DOUBLE, filetype, "native", MPI_INFO_NULL));
    MPICALL(MPI_File_write_all(
        fh, comb_field.data(), 1, memtype, MPI_STATUS_IGNORE));
    MPI_File_close(&fh);
    MPI_Type_free(&memtype);
    MPI_Type_free(&filetype);
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
    const auto dim = M::dim;
    const auto nblocks = t.dataptr.size();
    fassert(nblocks > 0);
    MIdx comb_start(std::numeric_limits<IntIdx>::max());
    MIdx comb_end(std::numeric_limits<IntIdx>::min());
    for (size_t b = 0; b < nblocks; ++b) {
      comb_start = comb_start.min(t.starts[b]);
      comb_end = comb_start.max(t.starts[b] + t.sizes[b]);
    }
    const MIdx comb_size = comb_end - comb_start;

    // XXX check the blocks on current rank constitute a rectangular subarray
    fassert_equal(comb_size % t.sizes[0], MIdx(0));
    fassert_equal((comb_size / t.sizes[0]).prod(), (int)nblocks);

    typename M::IndexCells comb_indexer(comb_start, comb_size);
    FieldCell<T> comb_field(comb_indexer);

    MPI_Datatype filetype;
    CreateSubarray(m.GetGlobalSize(), comb_size, comb_start, &filetype);

    MPI_Datatype memtype;
    CreateSubarray(comb_size, comb_size, MIdx(0), &memtype);

    const MPI_Comm comm = m.GetMpiComm();
    MPI_File fh;
    try {
      MPICALL(MPI_File_open(
          comm, path.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh));
    } catch (const std::exception& e) {
      throw std::runtime_error(
          std::string() + e.what() + ", while opening file '" + path + "'");
    }
    // fassert(error == MPI_SUCCESS, "Can't open file '" + path + "'");
    MPICALL(MPI_File_set_view(
        fh, 0, MPI_DOUBLE, filetype, "native", MPI_INFO_NULL));
    MPICALL(MPI_File_read_all(
        fh, comb_field.data(), 1, memtype, MPI_STATUS_IGNORE));
    MPI_File_close(&fh);
    MPI_Type_free(&memtype);
    MPI_Type_free(&filetype);

    for (size_t b = 0; b < nblocks; ++b) {
      const typename M::BlockCells block(t.starts[b], t.sizes[b]);
      size_t i = 0;
      for (auto c : GRangeIn<IdxCell, dim>(comb_indexer, block)) {
        (*t.dataptr[b])[i++] = comb_field[c];
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

template <class M>
void Raw<M>::WriteXmf(
    std::string xmfpath, std::string name, Type type, Vect origin, Vect spacing,
    MIdx dims, std::string rawpath) {
  std::ofstream f(xmfpath);
  f.precision(20);
  f << "<?xml version='1.0' ?>\n";
  f << "<!DOCTYPE Xdmf SYSTEM 'Xdmf.dtd' []>\n";
  f << "<Xdmf Version='2.0'>\n";

  auto getstr = [](auto v) {
    std::stringstream s;
    s.precision(20);
    bool first = true;
    for (size_t d = v.dim; d > 0;) {
      --d;
      if (!first) {
        s << ' ';
      } else {
        first = false;
      }
      s << v[d];
    }
    return s.str();
  };

  std::string numbertype;
  int precision = 0;
  switch (type) {
    case Type::UInt16:
      numbertype = "UShort";
      precision = 2;
      break;
    case Type::Float32:
      numbertype = "Float";
      precision = 4;
      break;
    case Type::Float64:
      numbertype = "Double";
      precision = 8;
      break;
    default:
      fassert(false);
  }

  f << " <Domain>\n";
  f << "   <Grid GridType='Uniform'>\n";
  f << "     <Topology TopologyType='3DCORECTMesh' Dimensions='";
  f << getstr(dims + MIdx(1)) << "'/>\n\n";

  f << "     <Geometry GeometryType='ORIGIN_DXDYDZ'>\n";
  f << "       <DataItem Name='Origin' Dimensions='" << dim
    << "' NumberType='Float' Precision='8' Format='XML'>\n";
  f << "         " << getstr(origin) << "\n";
  f << "       </DataItem>\n";
  f << "       <DataItem Name='Spacing' Dimensions='" << dim
    << "' NumberType='Float' Precision='8' Format='XML'>\n";
  f << "         " << getstr(spacing) << "\n";
  f << "       </DataItem>\n";
  f << "     </Geometry>\n\n";

  f << "     <Attribute Name='" << name
    << "' AttributeType='Scalar' Center='Cell'>\n";
  f << "       <DataItem Dimensions='" << getstr(dims);
  f << "' NumberType='" << numbertype << "' Precision='" << precision
    << "' Format='Binary'>\n";
  f << "        " << rawpath << '\n';
  f << "       </DataItem>\n";
  f << "     </Attribute>\n";
  f << "   </Grid>\n";
  f << " </Domain>\n";
  f << "</Xdmf>\n";
}

template <class M>
void Raw<M>::WriteXmf(
    std::string xmfpath, std::string name, Type type, std::string rawpath,
    const M& m) {
  const Vect origin(0);
  const Vect spacing = m.GetCellSize();
  const auto dims = m.GetGlobalSize();
  WriteXmf(xmfpath, name, type, origin, spacing, dims, rawpath);
}

} // namespace dump
