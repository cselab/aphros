#pragma once

#include <string>
#include <iostream>
#include <memory>
#include <cassert>
#include <fstream>
#include <functional>

#include "mesh.hpp"

namespace output {

class EntryGeneric {
  std::string name_;
 public:
  EntryGeneric(std::string name)
      : name_(name)
  {}
  virtual std::string GetName() {
    return name_;
  }
  virtual void Prepare() = 0;
  virtual ~EntryGeneric() {}
};

template <class FieldType>
class EntryField : public EntryGeneric {
 public:
  EntryField(std::string name)
      : EntryGeneric(name)
  {}
  virtual const FieldType& GetField() = 0;
};

template <class FieldType>
class EntryFieldCopy : public EntryField<FieldType> {
  const FieldType& field_;
 public:
  EntryFieldCopy(std::string name, const FieldType& field)
      : EntryField<FieldType>(name)
      , field_(field)
  {}
  void Prepare() override {}
  const FieldType& GetField() override {
    return field_;
  }
};

template <class Vect, class IdxType>
class EntryExtractScalar :
    public EntryField<GField<typename Vect::value_type, IdxType>> {
  using VectField = GField<Vect, IdxType>;
  using ScalarField = GField<typename Vect::value_type, IdxType>;
  const VectField& vect_field_;
  ScalarField scalar_field_;
  size_t component_number_;
 public:
  EntryExtractScalar(std::string name, const VectField& vect_field,
                     size_t component_number)
      : EntryField<ScalarField>(name)
      , vect_field_(vect_field)
      , scalar_field_(vect_field_.size())
      , component_number_(component_number)
  {}
  void Prepare() override {
    for (size_t i = 0; i < vect_field_.size(); ++i) {
      IdxType idx(i);
      scalar_field_[idx] = vect_field_[idx][component_number_];
    }
  }
  const ScalarField& GetField() override {
    return scalar_field_;
  }
};

template <class Value, class Idx, class Mesh>
class EntryFunction :
    public EntryField<GField<Value, Idx>> {
  using ScalarField = GField<Value, Idx>;
  using Function = std::function<Value (Idx)>;
  const Mesh& mesh_;
  ScalarField scalar_field_;
  Function function_;
 public:
  EntryFunction(std::string name, const Mesh& mesh, const Function& function)
      : EntryField<ScalarField>(name)
      , mesh_(mesh)
      , scalar_field_(mesh_)
      , function_(function)
  {}
  void Prepare() override {
    for (auto idx : mesh_.template Get<Idx>()) {
      scalar_field_[idx] = function_(idx);
    }
  }
  const ScalarField& GetField() override {
    return scalar_field_;
  }
};

template <class Value>
class EntryScalar : public EntryGeneric {
 public:
  EntryScalar(std::string name)
      : EntryGeneric(name)
  {}
  virtual Value GetValue() = 0;
};

template <class Value>
class EntryScalarFunction :
    public EntryScalar<Value> {
  using Function = std::function<Value ()>;
  Function function_;
 public:
  EntryScalarFunction(std::string name, const Function& function)
      : EntryScalar<Value>(name)
      , function_(function)
  {}
  void Prepare() override {}
  Value GetValue() override {
    return function_();
  }
};

using Content = std::vector<std::shared_ptr<EntryGeneric>>;

class Session {
 public:
  virtual void Write(double time = 0., std::string title = "") = 0;
  virtual ~Session() {}
};

namespace plain {

// Requires structured mesh
template <class Mesh>
class SessionPlain : public Session {
  using Scal = typename Mesh::Scal;
  using MIdx = typename Mesh::MIdx;
  const Mesh& mesh;
  Content content_;
  std::ostream& out_;
  std::ofstream output_file_;
  void WriteHeader() {
    out_ << "Cell: ";
    for (auto& entry_generic : content_)
    if (auto entry = dynamic_cast<EntryField<FieldCell<Scal>>*>(
        entry_generic.get())) {
      out_ << entry->GetName() << " ";
    }
    out_ << "\n";
    out_ << "Face: ";
    for (auto& entry_generic : content_)
    if (auto entry = dynamic_cast<EntryField<FieldFace<Scal>>*>(
        entry_generic.get())) {
      out_ << entry->GetName() << " ";
    }
    out_ << "\n";
    out_ << "Node: ";
    for (auto& entry_generic : content_)
    if (auto entry = dynamic_cast<EntryField<FieldNode<Scal>>*>(
        entry_generic.get())) {
      out_ << entry->GetName() << " ";
    }
    out_ << "\n";

    out_ << "Dim: ";
    MIdx dim = mesh.GetBlockCells().GetDimensions();
    for (size_t i = 0; i < dim.size(); ++i) {
      out_ << dim[i] << " ";
    }
    out_ << std::endl;
  }
  void WriteFooter() {
    out_.flush();
  }
  template <class FieldType>
  bool TryWriteField(EntryGeneric* entry_generic) {
    if (auto entry = dynamic_cast<EntryField<FieldType>*>(
        entry_generic)) {
      auto& field = entry->GetField();
      for (size_t i = 0; i < field.size(); ++i) {
        out_ << field[typename FieldType::IdxType(i)] << " ";
      }
      out_ << "\n";
      return true;
    }
    return false;
  }
 public:
  SessionPlain(const Content& content, std::string filename, const Mesh& mesh)
      : mesh(mesh)
      , content_(content)
      , out_(output_file_)
  {
    output_file_.open(filename);
    WriteHeader();
  }
  void Write(double time = 0., std::string title = "") override {
    out_ << "Zone: time " << time
        << ", title " << title << std::endl;
    for (auto& entry_generic : content_) {
      entry_generic->Prepare();
    }
    for (auto& entry_generic : content_) {
      TryWriteField<FieldCell<Scal>>(entry_generic.get());
    }
    for (auto& entry_generic : content_) {
      TryWriteField<FieldFace<Scal>>(entry_generic.get());
    }
    for (auto& entry_generic : content_) {
      TryWriteField<FieldNode<Scal>>(entry_generic.get());
    }
    out_ << std::endl;
  }
  ~SessionPlain() {
    WriteFooter();
  }
};

template <class Scal>
class SessionPlainScalar : public Session {
 public:
  SessionPlainScalar(const Content& content, std::string filename)
      : content_(content)
      , out_(output_file_)
  {
    output_file_.open(filename);
    for (auto& entry_generic : content_) {
      out_ << entry_generic->GetName() << " ";
    }
    out_ << std::endl;
  }
  void Write(double /*time = 0.*/, std::string /*title = ""*/) override {
    out_.precision(16);
    for (auto& eg : content_) { // entry generic
      eg->Prepare();
      if (auto e = dynamic_cast<EntryScalar<Scal>*>(eg.get())) { 
        out_ << e->GetValue() << " ";
      } else if (auto e = dynamic_cast<EntryScalar<int>*>(eg.get())) {
        out_ << e->GetValue() << " ";
      } else {
        throw std::runtime_error(
            "SessionPlainScalar: Unknown entry type");
      }
    }

    out_ << std::endl;
  }
 private:
  Content content_;
  std::ostream& out_;
  std::ofstream output_file_;
};


} // namespace plain

template <class Scal>
using SessionPlainScalar = plain::SessionPlainScalar<Scal>;

template <class Mesh>
using SessionPlain = plain::SessionPlain<Mesh>;

} // namespace output
