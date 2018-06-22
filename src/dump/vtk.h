#pragma once

#include <string>
#include <fstream>
#include <vector>

#include "geom/vect.h"

// Writes legacy vtk polydata 
// xx: points
// pp: polygons as lists of indices
// fn: path
// cm: comment
template <class Vect>
void WriteVtkPoly(const std::vector<Vect>& xx, 
                  const std::vector<std::vector<size_t>>& pp,  
                  const std::string& fn,
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

// Writes legacy vtk polydata 
// vv: polygons as lists of points
// fn: path
// cm: comment
template <class Vect>
void WriteVtkPoly(const std::vector<std::vector<Vect>>& vv,  
                  const std::string& fn,
                  const std::string& cm="") {
  std::vector<Vect> xx;
  std::vector<std::vector<size_t>> pp;
  Convert(vv, xx, pp);
  WriteVtkPoly(xx, pp, fn, cm);
}
