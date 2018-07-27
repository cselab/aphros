#pragma once

#include <string>
#include <iostream>
#include <memory>
#include <cassert>
#include <fstream>
#include <functional>

#include "geom/mesh.h"

namespace output {

class Out {
 public:
  // n: name
  Out(std::string n) : n_(n) {}

  virtual std::string GetName() { return n_; }

  virtual void Prepare() = 0;

  virtual ~Out() {}

 private:
  std::string n_;
};

// F: field
template <class F>
class OutFld : public Out {
 public:
  OutFld(std::string n) : Out(n) {}

  virtual const F& GetField() = 0;
};

// V: value type
// I: index type
// M: m
template <class V, class I, class M>
class OutFldFunc : public OutFld<GField<V, I>> {
 public:
  using F = GField<V, I>;

  // Constructor.
  // n: name
  // m: mesh
  // f: function returning output value
  OutFldFunc(std::string n, const M& m, const std::function<V(I)>& u)
      : OutFld<F>(n), m(m), f_(m), u_(u) {}

  // Updates field from function.
  void Prepare() override {
    for (auto i : m.template GetIn<I>()) {
      f_[i] = u_(i);
    }
  }
  const F& GetField() override {
    return f_;
  }

 private:
  const M& m; // mesh
  F f_; // field
  std::function<V(I)> u_; // function returning output value
};

template <class V>
class OutScal : public Out {
 public:
  OutScal(std::string n) : Out(n) {}
  virtual V GetValue() = 0;
};

template <class V>
class OutScalFunc : public OutScal<V> {
 public:
  // Constructor.
  // n: name
  // u: function returning output value
  OutScalFunc(std::string n, const std::function<V()>& u)
      : OutScal<V>(n), u_(u) {}

  void Prepare() override {}

  V GetValue() override { return u_(); }

 private:
  std::function<V()> u_; // function for single value
};

using Content = std::vector<std::shared_ptr<Out>>;

class Session {
 public:
  virtual void Write(double time = 0., std::string title = "") = 0;
  virtual ~Session() {}
};

namespace plain {

// Requires structured m
template <class M>
class SessionPlain : public Session {
  using Scal = typename M::Scal;
  using MIdx = typename M::MIdx;
  const M& m;
  Content content_;
  std::ostream& out_;
  std::ofstream output_file_;
  void WriteHeader() {
    out_ << "Cell: ";
    for (auto& entry_generic : content_)
    if (auto entry = dynamic_cast<OutFld<FieldCell<Scal>>*>(
        entry_generic.get())) {
      out_ << entry->GetName() << " ";
    }
    out_ << "\n";
    out_ << "Face: ";
    for (auto& entry_generic : content_)
    if (auto entry = dynamic_cast<OutFld<FieldFace<Scal>>*>(
        entry_generic.get())) {
      out_ << entry->GetName() << " ";
    }
    out_ << "\n";
    out_ << "Node: ";
    for (auto& entry_generic : content_)
    if (auto entry = dynamic_cast<OutFld<FieldNode<Scal>>*>(
        entry_generic.get())) {
      out_ << entry->GetName() << " ";
    }
    out_ << "\n";

    out_ << "Dim: ";
    MIdx dim = m.GetInBlockCells().GetSize();
    for (size_t i = 0; i < dim.size(); ++i) {
      out_ << dim[i] << " ";
    }
    out_ << std::endl;
  }
  void WriteFooter() {
    out_.flush();
  }
  template <class FieldType>
  bool TryWriteField(Out* entry_generic) {
    if (auto entry = dynamic_cast<OutFld<FieldType>*>(
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
  SessionPlain(const Content& content, std::string filename, const M& m)
      : m(m)
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
      if (auto e = dynamic_cast<OutScal<Scal>*>(eg.get())) { 
        out_ << e->GetValue() << " ";
      } else if (auto e = dynamic_cast<OutScal<int>*>(eg.get())) {
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

template <class M>
using SessionPlain = plain::SessionPlain<M>;

} // namespace output
