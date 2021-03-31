// Created by Petr Karnakov on 09.02.2021
// Copyright 2021 ETH Zurich

#include <fstream>
#include <sstream>

#include "parse/template.h"
#include "util/logger.h"

#include "xmf.h"

namespace dump {

template <class Vect>
std::string Xmf<Vect>::GetXmfTemplate() {
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

template <class Vect>
template <class M>
auto Xmf<Vect>::GetMeta(const M& m, MIdx start, MIdx stride) -> Meta {
  Meta res;
  res.dimensions = m.GetGlobalSize();
  res.start = start;
  res.stride = stride;
  res.count = (res.dimensions - res.start) / stride;
  res.origin = Vect(start) * m.GetCellSize();
  res.spacing = m.GetCellSize();
  return res;
}

template <class Vect>
auto Xmf<Vect>::ReadXmf(std::istream& buf) -> Meta {
  const std::string fmt = GetXmfTemplate();
  std::stringstream sstr;
  sstr << buf.rdbuf();
  const std::string txt = sstr.str();
  const auto map = parse::ParseTemplate(fmt, txt);

  auto get = [&map](std::string key) {
    fassert(map.count(key), "Key not found: " + key);
    return map.at(key);
  };
  auto conv = [&get](std::string key, auto& value) {
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

  int fdim;
  conv("dim", fdim);
  fassert_equal(fdim, Vect::dim);

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

template <class Vect>
auto Xmf<Vect>::ReadXmf(const std::string& xmfpath) -> Meta {
  std::ifstream f(xmfpath);
  fassert(f.good(), "Can't open file '" + xmfpath + "' for reading");
  return ReadXmf(f);
}

template <class Vect>
void Xmf<Vect>::WriteXmf(std::ostream& buf, const Meta& meta) {
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
  conv("dim", dim);

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

  buf << parse::SubstituteTemplate(GetXmfTemplate(), map);
}

template <class Vect>
void Xmf<Vect>::WriteXmf(const std::string& xmfpath, const Meta& meta) {
  std::ofstream f(xmfpath);
  fassert(f.good(), "Can't open file '" + xmfpath + "' for writing");
  WriteXmf(f, meta);
}

template <class Vect>
constexpr size_t Xmf<Vect>::dim;

} // namespace dump
