#pragma once

#include <string>
#include <fstream>
#include <vector>

#include "geom/vect.h"

// Writes legacy vtk polydata 
// xx: points
// pp: polygons as lists of indices
// dd: list of scalar cell-based datasets of size pp.size()
// dn: list of names for dd
// fn: path
// cm: comment
template <class Vect, class Scal=typename Vect::value_type>
void WriteVtkPoly(const std::string& fn, 
                  const std::vector<Vect>& xx, 
                  const std::vector<std::vector<size_t>>& pp,  
                  const std::vector<const std::vector<Scal>*>& dd,
                  const std::vector<std::string>& dn,
                  const std::string& cm="") {
  std::ofstream f(fn.c_str());
  f << "# vtk DataFile Version 2.0\n";
  f << cm << "\n";
  f << "ASCII\n";
  f << "DATASET POLYDATA\n";

  f << "POINTS " <<  xx.size() << " float\n";
  for (auto& x : xx) {
    f << x[0] << " " << x[1] << " " << x[2] << "\n";
  }

  size_t np = 0; // total number of vortices
  for (auto& p : pp) {
    np += p.size();
  }
  f << "POLYGONS " << pp.size() << " " << (np + pp.size()) << "\n";
  for (auto& p : pp) {
    f << p.size();
    for (auto& k : p) {
      f << " " << k;
    }
    f << "\n";
  }

  // cell-based datasets
  for (size_t i = 0; i < dd.size(); ++i) {
    auto& d = *(dd[i]);
    auto& n = dn[i];
    f << "CELL_DATA " << d.size() << "\n"
        << "SCALARS " << n << " float\n"
        << "LOOKUP_TABLE default\n";
    for (auto& a : d) {
      f << a << "\n";
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
template <class Vect, class Scal=typename Vect::value_type>
void WriteVtkPoly(const std::string& fn,
                  const std::vector<std::vector<Vect>>& vv,  
                  const std::vector<const std::vector<Scal>*>& dd,
                  const std::vector<std::string>& dn,
                  const std::string& cm) {
  std::vector<Vect> xx;
  std::vector<std::vector<size_t>> pp;
  Convert(vv, xx, pp);
  WriteVtkPoly(fn, xx, pp, dd, dn, cm);
}
