#pragma once

#include <string>
#include <fstream>
#include <vector>

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
}

// Converts to index representation.
// vv: polygons as lists of points
// Returns:
// xx: points
// pp: polygons as lists of indices
template <class Vect>
void Convert(const std::vector<std::vector<Vect>>& vv, 
             std::vector<Vect>& xx, 
             std::vector<std::vector<size_t>>& pp) {
  xx.resize(0);
  pp.resize(0);
  for (auto& v : vv) {
    pp.emplace_back();
    for (auto& x : v) {
      pp.back().push_back(xx.size());
      xx.push_back(x);
    }
  }
}

// Writes legacy vtk polydata with cell-based dataset.
// vv: polygons as lists of points
// dd: cell-based dataset of size vv.size()
// dn: name of dataset
// fn: path
// cm: comment
// poly: true: polygons, false: lines
template <class Vect, class Scal=typename Vect::value_type>
void WriteVtkPoly(const std::string& fn,
                  const std::vector<std::vector<Vect>>& vv,  
                  const std::vector<const std::vector<Scal>*>& dd,
                  const std::vector<std::string>& dn,
                  const std::string& cm, bool poly, bool binary) {
  std::vector<Vect> xx;
  std::vector<std::vector<size_t>> pp;
  Convert(vv, xx, pp);
  WriteVtkPoly(fn, xx, pp, dd, dn, cm, poly, binary);
}
