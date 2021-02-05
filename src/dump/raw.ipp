// Created by Petr Karnakov on 01.02.2021
// Copyright 2021 ETH Zurich

#include <cstdint>
#include <fstream>
#include <limits>
#include <map>
#include <sstream>

#include "raw.h"
#include "util/distr.h"
#include "util/format.h"
#include "util/logger.h"

namespace dump {

template <class Type>
MPI_Datatype GetMpiType(Type type) {
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
void CreateSubarray(
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

    const MPI_Datatype metatype = GetMpiType(meta.type);
    MPI_Datatype filetype;
    CreateSubarray(
        m.GetGlobalSize(), comb_size, comb_start, metatype, &filetype);

    MPI_Datatype memtype;
    CreateSubarray(comb_size, comb_size, MIdx(0), metatype, &memtype);

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
    MPICALL(
        MPI_File_set_view(fh, 0, metatype, filetype, "native", MPI_INFO_NULL));

    auto cast = [type = meta.type](T value, void* ptr) {
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

    const size_t metatypesize = GetPrecision(meta.type);
    std::vector<char> buf(comb_size.prod() * metatypesize);
    for (size_t b = 0; b < nblocks; ++b) {
      const typename M::BlockCells block(t.starts[b], t.sizes[b]);
      size_t i = 0;
      for (auto c : GRangeIn<IdxCell, dim>(comb_indexer, block)) {
        cast(t.data[b][i++], &buf[c.raw() * metatypesize]);
      }
    }

    MPICALL(MPI_File_write_all(fh, buf.data(), 1, memtype, MPI_STATUS_IGNORE));
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

    const MPI_Datatype metatype = GetMpiType(meta.type);

    MPI_Datatype filetype;
    CreateSubarray(
        meta.dimensions, comb_size, meta.start + comb_start, metatype,
        &filetype);

    MPI_Datatype memtype;
    CreateSubarray(comb_size, comb_size, MIdx(0), metatype, &memtype);

    const MPI_Comm comm = m.GetMpiComm();
    MPI_File fh;
    try {
      MPICALL(MPI_File_open(
          comm, path.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh));
    } catch (const std::exception& e) {
      throw std::runtime_error(
          std::string() + e.what() + ", while opening file '" + path + "'");
    }

    MPICALL(MPI_File_set_view(
        fh, meta.seek, metatype, filetype, "native", MPI_INFO_NULL));

    const size_t metatypesize = GetPrecision(meta.type);
    std::vector<char> buf(comb_size.prod() * metatypesize);
    MPICALL(MPI_File_read_all(fh, buf.data(), 1, memtype, MPI_STATUS_IGNORE));
    MPI_File_close(&fh);
    MPI_Type_free(&memtype);
    MPI_Type_free(&filetype);

    auto cast = [type = meta.type](const void* ptr) -> T {
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
        (*t.dataptr[b])[i++] = cast(&buf[c.raw() * metatypesize]);
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
   <Grid Name="mesh" GridType="Uniform">
     <Topology TopologyType="{dim}DCORECTMesh" Dimensions="{nodes*}"/>
     <Geometry GeometryType="{geomtype}">
       <DataItem Name="Origin" Dimensions="{dim}" NumberType="Float" Precision="8" Format="XML">
         {origin*}
       </DataItem>
       <DataItem Name="Spacing" Dimensions="{dim}" NumberType="Float" Precision="8" Format="XML">
         {spacing*}
       </DataItem>
     </Geometry>
     <Attribute Name="{name}" AttributeType="Scalar" Center="Cell">
       <DataItem ItemType="HyperSlab" Dimensions="{countd*}" Type="HyperSlab">
           <DataItem Dimensions="3 {dim}" Format="XML">
             {start*}
             {stride*}
             {count*}
           </DataItem>
           <DataItem Dimensions="{dim*}" Seek="{seek}" Precision="{precision}" NumberType="{type}" Format="Binary">
             {binpath}
           </DataItem>
       </DataItem>
     </Attribute>
   </Grid>
 </Domain>
</Xdmf>
)EOF";
  std::string res;
  size_t i = 0;
  std::string key;
  // replace '{key*}' with '{key2} {key1} {key0}'
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
      if (key.length() && key.back() == '*') {
        const auto keybase = key.substr(0, key.length() - 1);
        for (size_t d = dim; d > 0;) {
          --d;
          res += '{' + keybase + std::to_string(d) + '}';
          if (d > 0) {
            res += ' ';
          }
        }
      } else {
        res += '{' + key + '}';
      }
    } else {
      res.push_back(fmt[i]);
      ++i;
    }
  }
  return res;
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
      fassert(map.count(key), "Key not found: " + key);
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
  auto convvect = [&conv](std::string key, auto& value) {
    for (size_t d = 0; d < dim; ++d) {
      conv(key + std::to_string(d), value[d]);
    }
  };

  Meta res;
  conv("name", res.name);
  conv("binpath", res.binpath);

  int dim;
  conv("dim", dim);
  fassert_equal(dim, M::dim);

  convvect("dim", res.dimensions);
  convvect("start", res.start);
  convvect("stride", res.stride);
  convvect("count", res.count);

  conv("seek", res.seek);

  convvect("spacing", res.spacing);
  convvect("origin", res.origin);

  res.type = StringToType(get("type"));
  int precision;
  conv("precision", precision);
  fassert_equal(precision, GetPrecision(res.type));

  return res;
}

template <class M>
auto Raw<M>::ReadXmf(const std::string& xmfpath) -> Meta {
  std::ifstream f(xmfpath);
  fassert(f.good(), "Can't open file '" + xmfpath + "' for reading");
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
  auto convvect = [&conv](std::string key, const auto& value) {
    for (size_t d = 0; d < dim; ++d) {
      conv(key + std::to_string(d), value[d]);
    }
  };

  conv("name", meta.name);
  conv("binpath", meta.binpath);
  conv("dim", M::dim);

  convvect("dim", meta.dimensions);
  convvect("start", meta.start);
  convvect("stride", meta.stride);
  convvect("nodes", meta.count + MIdx(1));
  convvect("count", meta.count);
  convvect("countd", meta.count);

  conv("seek", meta.seek);

  convvect("spacing", meta.spacing);
  convvect("origin", meta.origin);

  conv("type", TypeToString(meta.type));
  conv("precision", GetPrecision(meta.type));

  const std::string geomtypes[] = {
      "ORIGIN_DX",
      "ORIGIN_DXDY",
      "ORIGIN_DXDYDZ",
      "ORIGIN_DXDYDZDW",
  };
  conv("geomtype", geomtypes[dim - 1]);

  buf << Substitute(GetXmfTemplate(), map);
}

template <class M>
void Raw<M>::WriteXmf(const std::string& xmfpath, const Meta& meta) {
  std::ofstream f(xmfpath);
  fassert(f.good(), "Can't open file '" + xmfpath + "' for writing");
  WriteXmf(f, meta);
}

} // namespace dump
