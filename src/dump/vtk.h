// Created by Petr Karnakov on 17.06.2022
// Copyright 2022 ETH Zurich

#pragma once

#include <string>
#include <vector>

namespace dump {

template <class Vect_>
class Vtk {
 public:
  using Vect = Vect_;
  using Scal = typename Vect::Scal;

  // Writes legacy VTK polydata from polygons stored as lists of indices.
  // poly_indices: polygons as lists of indices in `points`
  // celldata: list of scalar cell-datasets of size `poly_indices.size()`
  // celldata_names: names of cell-datasets
  // pointdata: list of scalar point-datasets of size `points.size()`
  // pointdata_names: names of point-datasets
  // poly: write polygons (true) or lines (false)
  // binary: write binary file (true) or ASCII (false)
  static void WriteVtkPolyRaw(
      const std::string& path, const std::vector<Vect>& points,
      const std::vector<Vect>& normals,
      const std::vector<std::vector<size_t>>& poly_indices,
      const std::vector<const std::vector<Scal>*>& celldata,
      const std::vector<std::string>& celldata_names,
      const std::vector<const std::vector<Scal>*>& pointdata,
      const std::vector<std::string>& pointdata_names,
      const std::string& comment, bool poly, bool binary);

  // Writes legacy VTK polydata from polygons stored as lists of points.
  // poly_points: polygons as lists of points
  // poly_normals: normals, same shape as `poly_points`
  // celldata: list of scalar cell-datasets of size `poly_points.size()`
  // celldata_names: names of datasets
  // poly_pointdata: list of scalar point-datasets, the size of each array
  //                 equals the total number of points in `poly_points`; choice
  //                 of value in merged points is unspecified
  // pointdata_names: names of point-datasets
  // poly: write polygons (true) or lines (false)
  // binary: write binary file (true) or ASCII (false)
  // merge: merge close points (true)
  // merge_tol: tolerance for detection of close points
  static void WriteVtkPoly(
      const std::string& path,
      const std::vector<std::vector<Vect>>& poly_points,
      const std::vector<std::vector<Vect>>* poly_normals,
      const std::vector<const std::vector<Scal>*>& celldata,
      const std::vector<std::string>& celldata_names,
      const std::vector<const std::vector<Scal>*>& poly_pointdata,
      const std::vector<std::string>& pointdata_names,
      const std::string& comment, bool poly, bool binary, bool merge,
      Scal merge_tol = 1e-8);

  // Writes legacy VTK point data.
  // data: mapping from field name to data array.
  // Point coordinates will be read from data['x'], data['y'], data['z'].
  // Non-existing coordinate arrays will be replaces with zeros.
  static void WriteVtkPoints(
      const std::vector<std::pair<std::string, std::vector<Scal>>>& data,
      const std::string& path, const std::string& comment, bool binary);
};

} // namespace dump
