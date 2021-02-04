// Created by Petr Karnakov on 01.02.2021
// Copyright 2021 ETH Zurich

#include <fstream>
#include <limits>
#include <map>
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
std::string Raw<M>::GetXmfTemplate() {
  const std::string fmt = R"EOF(<?xml version="1.0" ?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>
<Xdmf Version="2.0">
 <Domain>
   <Grid GridType="Uniform">
     <Topology TopologyType="3DCORECTMesh" Dimensions="{nodes2} {nodes1} {nodes0}"/>
     <Geometry GeometryType="ORIGIN_DXDYDZ">
       <DataItem Name="Origin" Dimensions="3" NumberType="Float" Precision="8" Format="XML">
         {origin2} {origin1} {origin0}
       </DataItem>
       <DataItem Name="Spacing" Dimensions="3" NumberType="Float" Precision="8" Format="XML">
         {spacing2} {spacing1} {spacing0}
       </DataItem>
     </Geometry>
     <Attribute Name="{name}" AttributeType="Scalar" Center="Cell">
       <DataItem ItemType="HyperSlab" Dimensions="{countd2} {countd1} {countd0}" Type="HyperSlab">
           <DataItem Dimensions="3 3" Format="XML">
             {start2} {start1} {start0}
             {stride2} {stride1} {stride0}
             {count2} {count1} {count0}
           </DataItem>
           <DataItem Dimensions="{dim2} {dim1} {dim0}" Seek="{seek}" Precision="{precision}" NumberType="{type}" Format="Binary">
             {binpath}
           </DataItem>
       </DataItem>
     </Attribute>
   </Grid>
 </Domain>
</Xdmf>
)EOF";
  return fmt;
}

std::string Substitute(
    const std::string& fmt, const std::map<std::string, std::string>& map) {
  std::string res;
  size_t i = 0;
  std::string key;
  while (i < fmt.size()) {
    if (fmt[i] == '{') {
      key = "";
      ++i;
      for (; i < fmt.size(); ++i) {
        if (fmt[i] == '}') {
          break;
        }
        key += fmt[i];
      }
      ++i;
      fassert(map.count(key), "Key not found: " + key)
      res += map.at(key);
    } else {
      res.push_back(fmt[i]);
      ++i;
    }
  }
  return res;
}

std::map<std::string, std::string> Parse(
    const std::string& fmt, const std::string& txt) {
  std::map<std::string, std::string> res;
  size_t i = 0;
  size_t j = 0;
  std::string key;
  std::string value;
  while (i < fmt.size() && j < txt.size()) {
    if (fmt[i] == '{') {
      if (fmt[i] == txt[j]) {
        std::stringstream msg;
        msg << "Expected different characters, got '" << fmt[i]
            << "' at position " << i << " of template and '" << txt[j]
            << "' at position " << j << " of input";
        fassert(false, msg.str());
      }
      key = "";
      value = "";
      ++i;
      for (; i < fmt.size(); ++i) {
        if (fmt[i] == '}') {
          break;
        }
        key += fmt[i];
      }
      ++i;
      while (i < fmt.size() && j < txt.size() && fmt[i] != txt[j]) {
        value += txt[j];
        ++j;
      }
      res[key] = value;
    } else {
      if (fmt[i] != txt[j]) {
        std::stringstream msg;
        msg << "Expected equal characters, got '" << fmt[i] << "' at position "
            << i << " of template and '" << txt[j] << "' at position " << j
            << " of input\n";
        fassert(false, msg.str());
      }
      ++i;
      ++j;
    }
  }
  return res;
}

template <class M>
auto Raw<M>::StringToType(std::string s) -> Type {
  Type res;
  if (s == "UShort") {
    res = Type::UInt16;
  } else if (s == "Float") {
    res = Type::Float32;
  } else if (s == "Double") {
    res = Type::Float64;
  } else {
    fassert(false, "Unknown type=" + s);
  }
  return res;
}

template <class M>
std::string Raw<M>::TypeToString(Type type) {
  switch (type) {
    case Type::UInt16:
      return "UShort";
    case Type::Float32:
      return "Float";
    case Type::Float64:
      return "Double";
    default:
      fassert(false);
  }
  return "";
}

template <class M>
int Raw<M>::GetPrecision(Type type) {
  switch (type) {
    case Type::UInt16:
      return 2;
    case Type::Float32:
      return 4;
    case Type::Float64:
      return 8; 
    default:
      fassert(false);
  }
  return 8;
}

template <class M>
auto Raw<M>::ReadXmf(std::istream& buf) -> Meta {
  const std::string fmt = GetXmfTemplate();
  std::stringstream sstr;
  sstr << buf.rdbuf();
  const std::string txt = sstr.str();
  const auto map = Parse(fmt, txt);

  auto get = [&map](std::string key) {
    fassert(map.count(key), "Key not found: " + key);
    return map.at(key);
  };
  auto conv = [&map, &get](std::string key, auto& value) {
    std::stringstream(get(key)) >> value;
  };

  Meta res;
  conv("name", res.name);
  conv("binpath", res.binpath);

  conv("dim0", res.dimensions[0]);
  conv("dim1", res.dimensions[1]);
  conv("dim2", res.dimensions[2]);

  conv("start0", res.start[0]);
  conv("start1", res.start[1]);
  conv("start2", res.start[2]);

  conv("stride0", res.stride[0]);
  conv("stride1", res.stride[1]);
  conv("stride2", res.stride[2]);

  conv("count0", res.count[0]);
  conv("count1", res.count[1]);
  conv("count2", res.count[2]);

  conv("seek", res.seek);

  conv("spacing0", res.spacing[0]);
  conv("spacing1", res.spacing[1]);
  conv("spacing2", res.spacing[2]);

  conv("origin0", res.origin[0]);
  conv("origin1", res.origin[1]);
  conv("origin2", res.origin[2]);

  res.type = StringToType(get("type"));
  int precision;
  conv("precision", precision);
  fassert_equal(precision, GetPrecision(res.type));

  return res;
}

template <class M>
auto Raw<M>::ReadXmf(const std::string& xmfpath) -> Meta {
  std::ifstream f(xmfpath);
  return ReadXmf(f);
}

template <class M>
auto Raw<M>::GetMeta(MIdx start, MIdx stride, const M& m) -> Meta {
  Meta res;
  res.dimensions = m.GetGlobalSize();
  res.start = start;
  res.stride = stride;
  res.count = (res.dimensions - res.start) / stride;
  res.origin = Vect(start) * m.GetCellSize();
  res.spacing = m.GetCellSize();
  return res;
}

template <class M>
void Raw<M>::WriteXmf(std::ostream& buf, const Meta& meta) {
  std::map<std::string, std::string> map;

  auto conv = [&map](std::string key, const auto& value) {
    std::stringstream s;
    s.precision(16);
    s << value;
    map[key] = s.str();
  };

  conv("name", meta.name);
  conv("binpath", meta.binpath);

  conv("dim0", meta.dimensions[0]);
  conv("dim1", meta.dimensions[1]);
  conv("dim2", meta.dimensions[2]);

  conv("start0", meta.start[0]);
  conv("start1", meta.start[1]);
  conv("start2", meta.start[2]);

  conv("stride0", meta.stride[0]);
  conv("stride1", meta.stride[1]);
  conv("stride2", meta.stride[2]);

  conv("nodes0", meta.count[0] + 1);
  conv("nodes1", meta.count[1] + 1);
  conv("nodes2", meta.count[2] + 1);
  conv("count0", meta.count[0]);
  conv("count1", meta.count[1]);
  conv("count2", meta.count[2]);
  conv("countd0", meta.count[0]);
  conv("countd1", meta.count[1]);
  conv("countd2", meta.count[2]);

  conv("seek", meta.seek);

  conv("spacing0", meta.spacing[0]);
  conv("spacing1", meta.spacing[1]);
  conv("spacing2", meta.spacing[2]);

  conv("origin0", meta.origin[0]);
  conv("origin1", meta.origin[1]);
  conv("origin2", meta.origin[2]);

  conv("type", TypeToString(meta.type));
  conv("precision", GetPrecision(meta.type));

  buf << Substitute(GetXmfTemplate(), map);
}

template <class M>
void Raw<M>::WriteXmf(const std::string& xmfpath, const Meta& meta) {
  std::ofstream f(xmfpath);
  WriteXmf(f, meta);
}

} // namespace dump
