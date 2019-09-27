#pragma once

#include <string>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>

// Returns true if machine is little endian
static inline bool IsLittleEnd() {
  int a = 1;
  return *(char*)&a == 1;
}

// Swap endiannes
template <class T>
T SwapEnd(T a) {
  char* p = (char*)(&a);
  size_t st = sizeof(T);
  for (size_t i = 0; i < st / 2; ++i) {
    std::swap(p[i], p[st - i - 1]);
  }
  return a;
}

// Writes legacy vtk polydata 
// xx: points
// nn: normals
// pp: polygons as lists of indices
// dd: list of scalar cell-based datasets of size pp.size()
// dn: list of names for dd
// fn: path
// cm: comment
// poly: true: polygons, false: lines
// bin: true: binary, false: ascii
template <class Vect, class Scal=typename Vect::value_type>
void WriteVtkPoly(const std::string& fn,
                  const std::vector<Vect>& xx,
                  const std::vector<Vect>& nn,
                  const std::vector<std::vector<size_t>>& pp,
                  const std::vector<const std::vector<Scal>*>& dd,
                  const std::vector<std::string>& dn,
                  const std::string& cm, bool poly, bool bin) {
  std::ofstream f(fn.c_str(), std::ios::binary);

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
  f << cm << "\n";
  f << (bin ? "BINARY" : "ASCII") << "\n";
  f << "DATASET POLYDATA\n";

  f << "POINTS " <<  xx.size() << " float\n";
  if (bin) {
    for (auto& x : xx) {
      WFloat(x[0]);
      WFloat(x[1]);
      WFloat(x[2]);
    }
    f << "\n";
  } else {
    for (auto& x : xx) {
      f << x[0] << " " << x[1] << " " << x[2] << "\n";
    }
  }

  size_t np = 0; // total number of vortices
  for (auto& p : pp) {
    np += p.size();
  }
  f << (poly ? "POLYGONS" : "LINES");
  f << " " << pp.size() << " " << (np + pp.size()) << "\n";
  if (bin) {
    for (auto& p : pp) {
      WInt(p.size());
      for (auto& k : p) {
        WInt(k);
      }
    }
    f << "\n";
  } else {
    for (auto& p : pp) {
      f << p.size();
      for (auto& k : p) {
        f << " " << k;
      }
      f << "\n";
    }
  }

  // check data size
  for (auto& d : dd) {
    (void) d;
    assert(d->size() == pp.size());
  }

  // cell-based datasets
  for (size_t i = 0; i < dd.size(); ++i) {
    auto& d = *(dd[i]);
    auto& n = dn[i];
    if (i == 0) {
      f << "CELL_DATA " << d.size() << "\n";
    }
    f << "SCALARS " << n << " float\n"
        << "LOOKUP_TABLE default\n";
    if (bin) {
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
  if (!nn.empty()) {
    f << "POINT_DATA " << nn.size() << "\n";
    f << "VECTORS " << "Normals" << " float\n" << "LOOKUP_TABLE default\n";
    if (bin) {
      for (auto& n : nn) {
        WFloat(n[0]);
        WFloat(n[1]);
        WFloat(n[2]);
      }
    } else {
      for (auto& n : nn) {
        f << n[0] << " " << n[1] << " " << n[2] << "\n";
      }
    }
  }
}

// Converts to index representation.
// vv: polygons as lists of points
// Returns:
// xx: points
// pp: polygons as lists of indices
template <class Vect>
void Convert(const std::vector<std::vector<Vect>>& vv, 
             const std::vector<std::vector<Vect>>* vvn, 
             std::vector<Vect>& xx, 
             std::vector<Vect>& nn, 
             std::vector<std::vector<size_t>>& pp) {
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

template <class V, class T=typename V::value_type>
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
// Returns:
// xx: points
// pp: polygons as lists of indices
template <class Vect>
void ConvertMerge(
             const std::vector<std::vector<Vect>>& vv, 
             const std::vector<std::vector<Vect>>* vvn, 
             std::vector<Vect>& xx, 
             std::vector<Vect>& nn, 
             std::vector<std::vector<size_t>>& pp) {
  struct HashPoint {
    size_t operator()(const Vect& x) const noexcept {
      const int m = 1000000;
      size_t h0 = std::hash<int>{}(int(x[0] * m));
      size_t h1 = std::hash<int>{}(int(x[1] * m));
      size_t h2 = std::hash<int>{}(int(x[2] * m));
      return h0 ^ (h1 << 1) ^ (h2 << 2);
    }
  };

  std::unordered_map<Vect, size_t, HashPoint> s;
  xx.resize(0);
  nn.resize(0);
  pp.resize(0);
  for (size_t i = 0; i < vv.size(); ++i) {
    auto& v = vv[i];
    pp.emplace_back();
    for (size_t j = 0; j < v.size(); ++j) {
      auto& x = v[j];
      if (!s.count(x)) {
        s[x] = xx.size();
        xx.push_back(x);
        if (vvn) {
          nn.push_back((*vvn)[i][j]);
        }
      }
      pp.back().push_back(s[x]);
    }
  }
}


// Writes legacy vtk polydata with cell-based dataset.
// vv: polygons as lists of points
// vvn: normals, same shape as vv
// dd: cell-based dataset of size vv.size()
// dn: name of dataset
// fn: path
// cm: comment
// poly: true: polygons, false: lines
// binary: use binary format
// merge: combine close points
template <class Vect, class Scal=typename Vect::value_type>
void WriteVtkPoly(const std::string& fn,
                  const std::vector<std::vector<Vect>>& vv,
                  const std::vector<std::vector<Vect>>* vvn,
                  const std::vector<const std::vector<Scal>*>& dd,
                  const std::vector<std::string>& dn,
                  const std::string& cm, bool poly, bool binary,
                  bool merge) {
  std::vector<Vect> xx;
  std::vector<Vect> nn;
  std::vector<std::vector<size_t>> pp;
  if (merge) {
    ConvertMerge(vv, vvn, xx, nn, pp);
  } else {
    Convert(vv, vvn, xx, nn, pp);
  }
  WriteVtkPoly(fn, xx, nn, pp, dd, dn, cm, poly, binary);
}
