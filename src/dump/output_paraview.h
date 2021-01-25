// Created by Petr Karnakov on 25.04.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>

#include "dump/dumper.h"
#include "geom/mesh.h"
#include "output.h"

namespace output {

namespace vtk {

class SerVtk : public Ser {
 protected:
  VOut vo_;
  std::string tl_; // tl for collection
  std::string bn_; // basename
  std::ofstream oc_; // out collection (.pvd)
 public:
  SerVtk(const VOut& vo, std::string tl, std::string fn)
      : vo_(vo), tl_(tl), bn_(fn), oc_(bn_ + ".pvd") {}
};

template <class M>
class SerVtkStruct : public SerVtk {
 public:
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using MIdx = typename M::MIdx;

  // vo: instances of Out
  // tl: title
  // fn: filename
  SerVtkStruct(const VOut& vo, std::string tl, std::string fn, const M& m_)
      : SerVtk(vo, tl, fn), m(m_), fr_(0) {
    ColHead(oc_);
  }
  ~SerVtkStruct() {
    ColFoot(oc_);
  }
  void Write(double t, std::string /*tl*/) override {
    std::string fn = GetDumpName(bn_, ".vts", fr_);
    ColEntry(oc_, t, fn);
    CreateFile(fn);
    ++fr_;
  }

 private:
  void ArHead(std::ostream& out, std::string name, size_t dim) {
    out << "        <DataArray "
        << "Name=\"" << name << "\" "
        << "NumberOfComponents=\"" << dim << "\" "
        << "type=\"Float32\" format=\"ascii\">\n";
  }
  void ArFoot(std::ostream& out) {
    out << "        </DataArray>\n";
  }
  void ArData(std::ostream& out, OutFld<FieldCell<Scal>>* e) {
    auto& fc = e->GetField();

    ArHead(out, e->GetName(), 1);
    for (auto c : m.Cells()) {
      if (!std::isnan(fc[c])) {
        out << fc[c] << " ";
      } else {
        out << "0 ";
      }
    }
    out << "\n";
    ArFoot(out);
  }
  void ArData(std::ostream& out, OutFld<FieldNode<Scal>>* o) {
    auto& f = o->GetField();

    ArHead(out, o->GetName(), 1);
    for (auto n : m.Nodes()) {
      out << f[n] << " ";
    }
    out << "\n";
    ArFoot(out);
  }
  void FileHead(std::ostream& out) {
    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"StructuredGrid\" "
        << "version=\"0.1\" byte_order=\"LittleEndian\">\n";

    MIdx s = m.GetInBlockCells().GetSize();
    out << "  <StructuredGrid WholeExtent=\""
        << "0 " << s[0] << " 0 " << s[1] << " 0 " << (M::dim == 3 ? s[2] : 0)
        << "\">\n";

    out << "    <Piece Extent=\""
        << "0 " << s[0] << " 0 " << s[1] << " 0 " << (M::dim == 3 ? s[2] : 0)
        << "\">\n";
  }
  void FileFoot(std::ostream& out) {
    out << "    </Piece>\n";
    out << "  </StructuredGrid>\n";
    out << "</VTKFile>\n";
  }
  void FileData(std::ostream& out) {
    // node fields
    out << "      <PointData>\n";
    for (auto& og : vo_) {
      if (auto o = dynamic_cast<OutFld<FieldNode<Scal>>*>(og.get())) {
        o->Prepare();
        ArData(out, o);
      }
    }
    out << "      </PointData>\n";

    // cell fields
    out << "      <CellData>\n";
    for (auto& og : vo_) {
      if (auto o = dynamic_cast<OutFld<FieldCell<Scal>>*>(og.get())) {
        o->Prepare();
        ArData(out, o);
      }
    }
    out << "      </CellData>\n";

    // mesh nodes
    out << "      <Points>\n";
    ArHead(out, "m", 3);
    for (auto n : m.Nodes()) {
      Vect p = m.GetNode(n);
      for (size_t i = 0; i < 3; ++i) {
        out << (i < M::dim ? p[i] : 0.) << " ";
      }
    }
    ArFoot(out);
    out << "      </Points>\n";
  }
  void CreateFile(std::string fn) {
    std::ofstream out(fn);
    out.precision(20);
    FileHead(out);
    FileData(out);
    FileFoot(out);
  }
  void ColHead(std::ostream& out) {
    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"Collection\" "
        << "version=\"0.1\" byte_order=\"LittleEndian\">\n";
    out << "  <Collection>\n";
  }
  void ColFoot(std::ostream& out) {
    out << "  </Collection>\n";
    out << "</VTKFile>\n";
  }
  void ColEntry(std::ostream& out, double t, std::string fn) {
    out << "    <DataSet tstep=\"" << t << "\" "
        << "group=\"\" part=\"0\" "
        << "file=\"" << fn << "\"/>\n";
  }

  const M& m;
  size_t fr_; // frame index
};

} // namespace vtk

using SerVtk = vtk::SerVtk;

template <class M>
using SerVtkStruct = vtk::SerVtkStruct<M>;

} // namespace output
