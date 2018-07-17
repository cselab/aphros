#pragma once

#include <string>
#include <iostream>
#include <memory>
#include <cassert>
#include <fstream>
#include <fstream>

#include "output.h"
#include "geom/mesh.h"
#include "dump/dumper.h"

namespace output {

namespace paraview {

class SessionParaview : public Session {
 protected:
  Content content_;
  std::string title_;
  std::string filename_;
  std::ofstream collection_;
 public:
  SessionParaview(const Content& content,
                  std::string title,
                  std::string filename)
      : content_(content),
        title_(title),
        filename_(filename),
        collection_(filename_ + ".pvd") {
    collection_.sync_with_stdio(false);
  }
};

template <class Mesh>
class SessionParaviewStructured : public SessionParaview {
  using Scal = typename Mesh::Scal;
  using Vect = typename Mesh::Vect;
  using MIdx = typename Mesh::MIdx;

  const Mesh& mesh;
  size_t timestep_;
  void WriteDataArrayHeader(std::ostream& out,
                            std::string name,
                            size_t num_components) {
    out << "        <DataArray "
        << "Name=\"" << name << "\" "
        << "NumberOfComponents=\"" << num_components << "\" "
        << "type=\"Float32\" format=\"ascii\">\n";
  }
  void WriteDataArrayFooter(std::ostream& out) {
    out << "        </DataArray>\n";
  }
  void WriteField(std::ostream& out,
                  EntryField<FieldCell<Scal>>* entry) {
    auto& field = entry->GetField();

    WriteDataArrayHeader(out, entry->GetName(), 1);
    for (auto idxcell : mesh.Cells()) {
      out << field[idxcell] << " ";
    }
    out << "\n";
    WriteDataArrayFooter(out);
  }
  void WriteField(std::ostream& out,
                  EntryField<FieldNode<Scal>>* entry) {
    auto& field = entry->GetField();

    WriteDataArrayHeader(out, entry->GetName(), 1);
    for (auto idxnode : mesh.Nodes()) {
      out << field[idxnode] << " ";
    }
    out << "\n";
    WriteDataArrayFooter(out);
  }
  void WriteDataFileHeader(std::ostream& out) {
    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"StructuredGrid\" "
        << "version=\"0.1\" byte_order=\"LittleEndian\">\n";

    MIdx size = mesh.GetBlockCells().GetDimensions();
    out << "  <StructuredGrid WholeExtent=\""
        << "0 " << size[0]
        << " 0 " << size[1]
        << " 0 " << (Mesh::dim == 3 ? size[2] : 0) << "\">\n";

    out << "    <Piece Extent=\""
        << "0 " << size[0]
        << " 0 " << size[1]
        << " 0 " << (Mesh::dim == 3 ? size[2] : 0) << "\">\n";
  }
  void WriteDataFileFooter(std::ostream& out) {
    out << "    </Piece>\n";
    out << "  </StructuredGrid>\n";
    out << "</VTKFile>\n";
  }
  void WriteDataFileContent(std::ostream& out) {
    out << "      <PointData>\n";
    for (auto& entry_generic : content_) {
      if (auto entry = dynamic_cast<EntryField<FieldNode<Scal>>*>(
          entry_generic.get())) {
        entry->Prepare();
        WriteField(out, entry);
      }
    }
    out << "      </PointData>\n";

    out << "      <CellData>\n";
    for (auto& entry_generic : content_) {
      if (auto entry = dynamic_cast<EntryField<FieldCell<Scal>>*>(
          entry_generic.get())) {
        entry->Prepare();
        WriteField(out, entry);
      }
    }
    out << "      </CellData>\n";

    out << "      <Points>\n";
    WriteDataArrayHeader(out, "mesh", 3);
    for (auto idxnode : mesh.Nodes()) {
      Vect p = mesh.GetNode(idxnode);
      for (size_t i = 0; i < 3; ++i) {
        out << (i < Mesh::dim ? p[i] : 0.) << " ";
      }
    }
    WriteDataArrayFooter(out);
    out << "      </Points>\n";
  }
  void CreateDataFile(std::string datafile_name) {
    std::ofstream datafile(datafile_name);
    datafile.sync_with_stdio(false);
    WriteDataFileHeader(datafile);
    WriteDataFileContent(datafile);
    WriteDataFileFooter(datafile);
  }
  void WriteCollectionHeader(std::ostream& out) {
    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"Collection\" "
        << "version=\"0.1\" byte_order=\"LittleEndian\">\n";
    out << "  <Collection>\n";
  }
  void WriteCollectionFooter(std::ostream& out) {
    out << "  </Collection>\n";
    out << "</VTKFile>\n";
  }
  void WriteCollectionEntry(std::ostream& out,
                            double time,
                            std::string datafile_name) {
     out << "    <DataSet timestep=\"" << time << "\" "
         << "group=\"\" part=\"0\" "
         << "file=\"" << datafile_name << "\"/>\n";
  }
 public:
  SessionParaviewStructured(const Content& content,
                            std::string title,
                            std::string filename,
                            const Mesh& mesh)
      : SessionParaview(content, title, filename),
        mesh(mesh),
        timestep_(0) {
    WriteCollectionHeader(collection_);
  }
  ~SessionParaviewStructured() {
    WriteCollectionFooter(collection_);
  }
  void Write(double time, std::string /*title*/) override {
    std::string fn = GetDumpName(filename_, ".vts", timestep_);
    WriteCollectionEntry(collection_, time, fn);
    CreateDataFile(fn);
    ++timestep_;
  }
};

} // namespace paraview

using SessionParaview = paraview::SessionParaview;

template <class Mesh>
using SessionParaviewStructured = paraview::SessionParaviewStructured<Mesh>;

} // namespace output
