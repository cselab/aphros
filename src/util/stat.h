#pragma once

#include <algorithm>
#include <functional>
#include <iosfwd>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include "geom/mesh.h"
#include "solver/embed.h"
#include "util/logger.h"

template <class M>
class Stat {
 public:
  enum class Reduction { sum, min, max };
  using Scal = typename M::Scal;

  struct Entry {
    std::string name;
    std::string desc;
    Reduction op;
    std::function<Scal()> func;
  };

  Stat(M& m, const Embed<M>* eb = nullptr) : m(m), eb_(eb) {}
  Stat(const Stat&) = delete;

  void Add(
      std::string name, std::string desc, Reduction op,
      std::function<Scal(IdxCell c, const M&)> func) {
    AddName(name);
    AddMeshLoop(name, desc, op, func, m);
  };
  void Add(
      std::string name, std::string desc, Reduction op,
      std::function<Scal(IdxCell c, const Embed<M>&)> func) {
    AddName(name);
    fassert(eb_, "Add: Can't add the entry to Stat created without Embed.");
    AddMeshLoop(name, desc, op, func, *eb_);
  };
  void Add(
      std::string name, std::string desc, Reduction op,
      std::function<Scal()> func) {
    AddName(name);
    entries_.emplace(name, Entry{name, desc, op, func});
  };
  void Update() {
    auto sem = m.GetSem();
    if (sem()) {
      for (auto& p : entries_) {
        values_[p.first] = p.second.func();
      }
      for (auto& p : values_) {
        m.Reduce(&p.second, "sum");
      }
    }
    if (sem()) {
    }
  }
  void WriteHeader(std::ostream& out) const {
    bool first = true;
    for (auto n : names_) {
      first || out << ' ', first = false;
      out << n;
    }
    out << std::endl;
  }
  void WriteValues(std::ostream& out) const {
    bool first = true;
    for (auto n : names_) {
      first || out << ' ', first = false;
      out << values_.at(n);
    }
    out << std::endl;
  }
  void WriteSummary(std::ostream& out) const {
    for (auto n : names_) {
      out << n << ": " << entries_.at(n).desc << '\n';
    }
  }
  void SortNames() {
    std::sort(names_.begin(), names_.end());
  }
  const Entry& GetEntry(std::string name) const {
    return entries_[name];
  }
  const std::vector<std::string>& GetNames() const {
    return names_;
  }
  std::string GetDesc(std::string name) const {
    return entries_[name].desc;
  }

 private:
  void CheckExisting(std::string name) const {
    fassert(!entries_.count(name), "Add: Entry '" + name + "' already exists.");
  }
  void AddName(std::string name) {
    CheckExisting(name);
    names_.push_back(name);
  }
  template <class MEB>
  void AddMeshLoop(
      std::string name, std::string desc, Reduction op,
      std::function<Scal(IdxCell c, const MEB&)> func, const MEB& meb) {
    switch (op) {
      case Reduction::sum:
        entries_.emplace(
            name, //
            Entry{name, desc, op, [func, &meb = meb]() {
                    Scal sum = 0;
                    for (auto c : meb.Cells()) {
                      sum += func(c, meb);
                    }
                    return sum;
                  }});
        break;
      case Reduction::max:
        entries_.emplace(
            name, //
            Entry{name, desc, op, [func, &meb = meb]() {
                    Scal max = -std::numeric_limits<Scal>::max();
                    for (auto c : meb.Cells()) {
                      max = std::max(max, func(c, meb));
                    }
                    return max;
                  }});
        break;
      case Reduction::min:
        entries_.emplace(
            name, //
            Entry{name, desc, op, [func, &meb = meb]() {
                    Scal min = std::numeric_limits<Scal>::max();
                    for (auto c : meb.Cells()) {
                      min = std::min(min, func(c, meb));
                    }
                    return min;
                  }});
        break;
    }
  }

  M& m;
  const Embed<M>* eb_;
  std::vector<std::string> names_;
  std::map<std::string, Entry> entries_;
  std::map<std::string, Scal> values_;
};
