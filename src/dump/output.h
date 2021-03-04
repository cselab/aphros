// Created by Petr Karnakov on 25.04.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <cassert>
#include <fstream>
#include <functional>
#include <iostream>
#include <memory>
#include <string>

#include "geom/mesh.h"

namespace output {

class Out {
 public:
  // n: name
  Out(std::string n) : n_(n) {}

  virtual std::string GetName() {
    return n_;
  }

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
  OutFldFunc(std::string n, const M& m_, const std::function<V(I)>& u)
      : OutFld<F>(n), m(m_), f_(m), u_(u) {}

  // Updates field from function.
  void Prepare() override {
    for (auto i : m.template GetRangeIn<I>()) {
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
  virtual V second() = 0;
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

  V second() override {
    return u_();
  }

 private:
  std::function<V()> u_; // function for single value
};

using VOut = std::vector<std::shared_ptr<Out>>;

// Series
class Ser {
 public:
  virtual void Write(double time, std::string title) = 0;
  virtual ~Ser() {}
};

namespace plain {

template <class Scal>
class SerScalPlain : public Ser {
 public:
  // vo: instances of Out
  // fn: output filename
  SerScalPlain(const VOut& vo, std::string fn) : vo_(vo) {
    out_.open(fn);
    out_.precision(16);
    for (auto& o : vo_) {
      out_ << o->GetName() << " ";
    }
    out_ << std::endl;
  }
  void Write(double /*time*/, std::string /*title*/) override {
    for (auto& og : vo_) { // out generic
      og->Prepare();
      if (auto oscal = dynamic_cast<OutScal<Scal>*>(og.get())) {
        out_ << oscal->second() << " ";
      } else if (auto oint = dynamic_cast<OutScal<int>*>(og.get())) {
        out_ << oint->second() << " ";
      } else {
        fassert(false, "SerScalPlain: Unknown entry type");
      }
    }

    out_ << std::endl;
  }

 private:
  VOut vo_;
  std::ofstream out_;
};

} // namespace plain

template <class Scal>
using SerScalPlain = plain::SerScalPlain<Scal>;

} // namespace output
