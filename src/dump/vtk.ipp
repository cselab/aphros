// Created by Petr Karnakov on 20.06.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <algorithm>
#include <cassert>
#include <fstream>
#include <unordered_map>
#include <unordered_set>

#include "geom/vect.h"
#include "vtk.h"

// Returns true if machine is little endian.
static inline bool IsLittleEnd() {
  int a = 1;
  return *(char*)&a == 1;
}

// Swaps endiannes.
template <class T>
T SwapEnd(T a) {
  char* p = (char*)(&a);
  size_t st = sizeof(T);
  for (size_t i = 0; i < st / 2; ++i) {
    std::swap(p[i], p[st - i - 1]);
  }
  return a;
}

// Converts to index representation.
// vv: polygons as lists of points
// Returns:
// xx: points
// pp: polygons as lists of indices
template <class Vect>
void ConvertPointsToIndices(
    const std::vector<std::vector<Vect>>& vv,
    const std::vector<std::vector<Vect>>* vvn, std::vector<Vect>& xx,
    std::vector<Vect>& nn, std::vector<std::vector<size_t>>& pp) {
  xx.resize(0);
  nn.resize(0);
  pp.resize(0);
  for (auto& v : vv) {
    pp.emplace_back();
    for (auto& x : v) {
      pp.back().push_back(xx.size());
      xx.push_back(x);
    }
  }
  if (vvn) {
    for (auto& v : *vvn) {
      for (auto& x : v) {
        nn.push_back(x);
      }
    }
  }
}

template <class V, class T = typename V::value_type>
void RemoveDuplicatesUnordered(std::vector<V>& pp) {
  // returns hash up to permutation of elements in V
  struct HashPoly {
    size_t operator()(const V& p) const noexcept {
      size_t h = 1;
      for (auto i : p) {
        h *= std::hash<T>{}(i);
      }
      return h;
    }
  };
  std::unordered_set<V, HashPoly> ppu;
  for (auto& p : pp) {
    ppu.insert(p);
  }
  pp.clear();
  for (auto p : ppu) {
    pp.push_back(p);
  }
}

// Converts to index representation merging closely located points.
// vv: polygons as lists of points
// tol: positive tolerance for matching coordinates of points
// Returns:
// xx: points
// pp: polygons as lists of indices
template <class Scal, class Vect = generic::Vect<Scal, 3>>
void ConvertPointsToIndicesMerged(
    const std::vector<std::vector<Vect>>& vv,
    const std::vector<std::vector<Vect>>* vvn, Scal tol, std::vector<Vect>& xx,
    std::vector<Vect>& nn, std::vector<std::vector<size_t>>& pp) {
  static constexpr auto dim = Vect::dim;
  struct Hash {
    size_t operator()(const Vect& x) const noexcept {
      size_t res = 0;
      for (size_t d = 0; d < dim; ++d) {
        res ^= (std::hash<Scal>{}(x[d]) << d);
      }
      return res;
    }
  };

  std::unordered_map<Vect, size_t, Hash> index; // point to index in `xx`

  auto floor = [](Vect x) {
    for (size_t i = 0; i < dim; ++i) {
      x[i] = std::floor(x[i]);
    }
    return x;
  };

  // Returns a canonical representative of one cell of size `tol`.
  auto canonical = [tol, &floor](const Vect& x) -> Vect {
    return floor(x / tol) * tol;
  };

  // Returns pointer to element found by point up to tolerance `tol`.
  auto findtol = [&canonical, tol, &index](const Vect& x) -> size_t* {
    // iterate over {-1,1}^4
    for (unsigned i = 0; i < (1 << dim); ++i) {
      Vect dx;
      for (unsigned j = 0; j < dim; ++j) {
        dx[j] = (i & (1 << j)) ? -1 : 1;
      }
      auto it = index.find(canonical(x + dx * (tol * 0.25)));
      if (it != index.end()) {
        return &it->second;
      }
    }
    return nullptr;
  };

  xx.resize(0);
  nn.resize(0);
  pp.resize(0);
  for (size_t i = 0; i < vv.size(); ++i) {
    auto& v = vv[i];
    pp.emplace_back();
    for (size_t j = 0; j < v.size(); ++j) {
      auto& x = v[j];
      const auto* p = findtol(x);
      if (p) {
        pp.back().push_back(*p);
      } else {
        const auto nextindex = xx.size();
        index[canonical(x)] = nextindex;
        pp.back().push_back(nextindex);
        xx.push_back(x);
        if (vvn) {
          nn.push_back((*vvn)[i][j]);
        }
      }
    }
  }
}

namespace dump {

template <class Vect_>
void Vtk<Vect_>::WriteVtkPoly(
    const std::string& path, const std::vector<Vect>& points,
    const std::vector<Vect>& normals,
    const std::vector<std::vector<size_t>>& poly_indices,
    const std::vector<const std::vector<typename Vect::Scal>*>& celldata,
    const std::vector<std::string>& celldata_names, const std::string& comment,
    bool poly, bool binary) {
  std::ofstream f(path.c_str(), std::ios::binary);

  auto WInt = [&f](int a) {
    if (IsLittleEnd()) {
      a = SwapEnd(a);
    }
    f.write((char*)&a, sizeof(int));
  };
  auto WFloat = [&f](float a) {
    if (IsLittleEnd()) {
      a = SwapEnd(a);
    }
    f.write((char*)&a, sizeof(float));
  };

  f.precision(16);
  f << "# vtk DataFile Version 2.0\n";
  f << comment << "\n";
  f << (binary ? "BINARY" : "ASCII") << "\n";
  f << "DATASET POLYDATA\n";

  f << "POINTS " << points.size() << " float\n";
  if (binary) {
    for (auto& x : points) {
      WFloat(x[0]);
      WFloat(x[1]);
      WFloat(Vect::dim > 2 ? x[2] : 0);
    }
    f << "\n";
  } else {
    for (auto& x : points) {
      f << x[0] << " " << x[1] << " " << (Vect::dim > 2 ? x[2] : 0) << "\n";
    }
  }

  size_t np = 0; // total number of vortices
  for (auto& p : poly_indices) {
    np += p.size();
  }
  f << (poly ? "POLYGONS" : "LINES");
  f << " " << poly_indices.size() << " " << (np + poly_indices.size()) << "\n";
  if (binary) {
    for (auto& p : poly_indices) {
      WInt(p.size());
      for (auto& k : p) {
        WInt(k);
      }
    }
    f << "\n";
  } else {
    for (auto& p : poly_indices) {
      f << p.size();
      for (auto& k : p) {
        f << " " << k;
      }
      f << "\n";
    }
  }

  // check data size
  for (auto& d : celldata) {
    (void)d;
    assert(d->size() == poly_indices.size());
  }

  // cell-based datasets
  for (size_t i = 0; i < celldata.size(); ++i) {
    auto& d = *(celldata[i]);
    auto& n = celldata_names[i];
    if (i == 0) {
      f << "CELL_DATA " << d.size() << "\n";
    }
    f << "SCALARS " << n << " float\n"
      << "LOOKUP_TABLE default\n";
    if (binary) {
      for (auto& a : d) {
        WFloat(a);
      }
    } else {
      for (auto& a : d) {
        f << a << "\n";
      }
    }
  }

  // normals
  if (!normals.empty()) {
    f << "POINT_DATA " << normals.size() << "\n";
    f << "VECTORS "
      << "normals"
      << " float\n";
    if (binary) {
      for (auto& n : normals) {
        WFloat(n[0]);
        WFloat(n[1]);
        WFloat(Vect::dim ? n[2] : 0);
      }
    } else {
      for (auto& n : normals) {
        f << n[0] << " " << n[1] << " " << (Vect::dim > 2 ? n[2] : 0) << "\n";
      }
    }
  }
}

template <class Vect_>
void Vtk<Vect_>::WriteVtkPoly(
    const std::string& path, const std::vector<std::vector<Vect>>& poly_points,
    const std::vector<std::vector<Vect>>* poly_normals,
    const std::vector<const std::vector<Scal>*>& celldata,
    const std::vector<std::string>& celldata_names, const std::string& comment,
    bool poly, bool binary, bool merge, Scal merge_tol) {
  std::vector<Vect> xx;
  std::vector<Vect> nn;
  std::vector<std::vector<size_t>> pp;
  if (merge) {
    ConvertPointsToIndicesMerged(
        poly_points, poly_normals, merge_tol, xx, nn, pp);
  } else {
    ConvertPointsToIndices(poly_points, poly_normals, xx, nn, pp);
  }
  WriteVtkPoly(
      path, xx, nn, pp, celldata, celldata_names, comment, poly, binary);
}

} // namespace dump
