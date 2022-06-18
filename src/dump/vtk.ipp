// Created by Petr Karnakov on 20.06.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <algorithm>
#include <cassert>
#include <fstream>
#include <unordered_map>
#include <unordered_set>

#include "geom/vect.h"
#include "util/logger.h"
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

// Converts polygons to lists of indices.
// poly_points: polygons as lists of points
//
// Output:
// points: points
// poly_indices: polygons as lists of indices
template <class Scal, class Vect>
void ConvertPointsToIndices(
    const std::vector<std::vector<Vect>>& poly_points,
    const std::vector<std::vector<Vect>>* poly_normals,
    const std::vector<const std::vector<Scal>*>& poly_pointdata,
    /*out*/
    std::vector<Vect>& points, std::vector<Vect>& normals,
    std::vector<std::vector<size_t>>& poly_indices,
    std::vector<std::vector<Scal>>& pointdata) {
  points.clear();
  poly_indices.clear();
  normals.clear();
  pointdata.clear();

  for (auto& v : poly_points) {
    poly_indices.emplace_back();
    for (auto& x : v) {
      poly_indices.back().push_back(points.size());
      points.push_back(x);
    }
  }
  if (poly_normals) {
    normals.reserve(points.size());
    for (auto& v : *poly_normals) {
      for (auto& n : v) {
        normals.push_back(n);
      }
    }
  }
  pointdata.reserve(poly_pointdata.size());
  for (auto ptr : poly_pointdata) {
    pointdata.push_back(*ptr);
  }
}

template <class V, class T = typename V::value_type>
void RemoveDuplicatesUnordered(std::vector<V>& poly_indices) {
  // Returns hash up to permutation of elements in V.
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
  for (auto& p : poly_indices) {
    ppu.insert(p);
  }
  poly_indices.clear();
  for (auto p : ppu) {
    poly_indices.push_back(p);
  }
}

// Converts to index representation merging closely located points.
// poly_points: polygons as lists of points
// tol: positive tolerance for matching coordinates of points
//
// Output:
// points: points
// poly_indices: polygons as lists of indices
template <class Scal, class Vect>
void ConvertPointsToIndicesMerged(
    const std::vector<std::vector<Vect>>& poly_points,
    const std::vector<std::vector<Vect>>* poly_normals,
    const std::vector<const std::vector<Scal>*>& poly_pointdata, Scal tol,
    /*out*/
    std::vector<Vect>& points, std::vector<Vect>& normals,
    std::vector<std::vector<size_t>>& poly_indices,
    std::vector<std::vector<Scal>>& pointdata) {
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

  std::unordered_map<Vect, size_t, Hash> index; // Point to index in `points`.

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

  points.clear();
  normals.clear();
  poly_indices.clear();
  pointdata.clear();
  pointdata.resize(poly_pointdata.size());
  size_t cnt = 0; // Absolute index of point.
  for (size_t i = 0; i < poly_points.size(); ++i) {
    auto& v = poly_points[i];
    poly_indices.emplace_back();
    for (size_t j = 0; j < v.size(); ++j) {
      auto& x = v[j];
      const size_t* p = findtol(x);
      if (p) { // Use index of existing point.
        poly_indices.back().push_back(*p);
      } else { // Add new point.
        const size_t nextindex = points.size();
        index[canonical(x)] = nextindex;
        poly_indices.back().push_back(nextindex);
        points.push_back(x);
        for (size_t id = 0; id < poly_pointdata.size(); ++id) {
          pointdata[id].push_back((*poly_pointdata[id])[cnt]);
        }
        if (poly_normals) {
          normals.push_back((*poly_normals)[i][j]);
        }
      }
      ++cnt;
    }
  }
  // Check size of point-datasets.
  for (auto& d : poly_pointdata) {
    fassert_equal(d->size(), cnt);
  }
}

namespace dump {

template <class Vect_>
void Vtk<Vect_>::WriteVtkPolyRaw(
    const std::string& path, const std::vector<Vect>& points,
    const std::vector<Vect>& normals,
    const std::vector<std::vector<size_t>>& poly_indices,
    const std::vector<const std::vector<typename Vect::Scal>*>& celldata,
    const std::vector<std::string>& celldata_names,
    const std::vector<const std::vector<Scal>*>& pointdata,
    const std::vector<std::string>& pointdata_names, const std::string& comment,
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

  size_t np = 0; // Total number of points.
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

  // Check data size.
  for (auto& d : celldata) {
    fassert_equal(d->size(), poly_indices.size());
  }

  // Cell-datasets.
  for (size_t i = 0; i < celldata.size(); ++i) {
    const auto& d = *(celldata[i]);
    if (i == 0) {
      f << "CELL_DATA " << d.size() << "\n";
    }
    f << "SCALARS " << celldata_names[i] << " float\n"
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

  // Normals.
  if (!normals.empty()) {
    f << "POINT_DATA " << normals.size() << "\n";
    f << "VECTORS normals float\n";
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

  // Point-datasets.
  for (size_t i = 0; i < pointdata.size(); ++i) {
    const auto& d = *(pointdata[i]);
    if (i == 0) {
      f << "POINT_DATA " << d.size() << "\n";
    }
    f << "SCALARS " << pointdata_names[i] << " float\n"
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
}

template <class Vect_>
void Vtk<Vect_>::WriteVtkPoly(
    const std::string& path, const std::vector<std::vector<Vect>>& poly_points,
    const std::vector<std::vector<Vect>>* poly_normals,
    const std::vector<const std::vector<Scal>*>& celldata,
    const std::vector<std::string>& celldata_names,
    const std::vector<const std::vector<Scal>*>& poly_pointdata,
    const std::vector<std::string>& pointdata_names, const std::string& comment,
    bool poly, bool binary, bool merge, Scal merge_tol) {
  std::vector<Vect> points;
  std::vector<Vect> normals;
  std::vector<std::vector<size_t>> poly_indices;
  std::vector<std::vector<Scal>> pointdata;
  if (merge) {
    ConvertPointsToIndicesMerged(
        poly_points, poly_normals, poly_pointdata, merge_tol, //
        /*out*/ points, normals, poly_indices, pointdata);
  } else {
    ConvertPointsToIndices(
        poly_points, poly_normals, poly_pointdata, //
        /*out*/ points, normals, poly_indices, pointdata);
  }
  std::vector<const std::vector<Scal>*> pointdata_ptr;
  for (const auto& d : pointdata) {
    pointdata_ptr.push_back(&d);
  }
  WriteVtkPolyRaw(
      path, points, normals, poly_indices, celldata, celldata_names,
      pointdata_ptr, pointdata_names, comment, poly, binary);
}

} // namespace dump
