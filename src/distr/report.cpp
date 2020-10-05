// Created by Petr Karnakov on 21.07.2018
// Copyright 2018 ETH Zurich

#include <iomanip>
#include <sstream>
#include <vector>

#include "report.h"

// Splits string by whitespace, skips every second.
std::vector<std::string> SplitSkip(const std::string& str) {
  std::vector<std::string> names; // result
  std::stringstream buf(str);
  buf >> std::skipws;
  while (buf) {
    std::string s;
    buf >> s;
    names.push_back(s);
    buf >> s; // skip
  }
  return names;
}

struct Node {
  const double kNone = -1.; // empty time
  std::vector<Node> nodes; // children
  std::string stage; // stage name
  double time; // time
  Node(const std::string& s_, double t_) : stage(s_), time(t_) {}
  Node(const std::string& s_) : stage(s_), time(kNone) {}
  Node* Find(const std::string& s) {
    for (auto& n : nodes) {
      if (n.stage == s) {
        return &n;
      }
    }
    return nullptr;
  }
  Node* Add(const std::string& s) {
    nodes.emplace_back(s, kNone);
    return &nodes.back();
  }
  void Print(std::ostream& out, std::string pre, double ta) const {
    auto fl = out.flags();
    out << pre << stage << " [" << std::setprecision(5) << time << " s, "
        << std::setprecision(3) << 100. * time / ta << "%]" << std::endl;
    for (auto& n : nodes) {
      n.Print(out, pre + "|     ", ta);
    }
    out.flags(fl);
  }
  void FillTime() {
    double ts = 0.; // sum
    if (time == kNone) {
      for (auto& n : nodes) {
        if (n.time == kNone) {
          n.FillTime();
        }
        ts += n.time;
      }
      time = ts;
    }
  }
};

void ParseReport(
    const std::map<std::string, double>& timings, std::ostream& out) {
  std::vector<std::string> names_prev; // list of strings previous
  Node root("all");
  for (auto& it : timings) {
    std::vector<std::string> names = SplitSkip(it.first);
    if (names.empty() || (names.size() == 1 && names[0] == "")) {
      names = {"other"};
    }

    // collect names from root to last common node with names_prev
    size_t i = 0;
    Node* pos = &root;
    while (i < names.size() && i < names_prev.size() &&
           names[i] == names_prev[i]) {
      pos = pos->Find(names[i]);
      ++i;
    }
    while (i < names.size()) {
      pos = pos->Add(names[i]);
      ++i;
    }
    pos->time = it.second;
    names_prev.swap(names);
  }

  root.FillTime();
  root.Print(out, "", root.time);
}
