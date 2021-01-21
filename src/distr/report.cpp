// Created by Petr Karnakov on 21.07.2018
// Copyright 2018 ETH Zurich

#include <iomanip>
#include <sstream>
#include <vector>

#include "report.h"
#include "util/format.h"

// Splits string by whitespace, skips every second.
std::vector<std::string> SplitBySeparator(
    const std::string& str, const std::string& sep) {
  std::vector<std::string> res;
  size_t i = 0;
  while (true) {
    const size_t imatch = str.find(sep, i);
    if (imatch != std::string::npos) {
      res.push_back(str.substr(i, imatch - i));
      i = imatch + sep.length();
    } else {
      res.push_back(str.substr(i));
      break;
    }
  }
  return res;
}

namespace suspender_tree {

struct Node {
  Node(const std::string& name_, double time_) : name(name_), time(time_) {}
  Node* FindChild(const std::string& name_) {
    for (auto& node : children) {
      if (node.name == name_) {
        return &node;
      }
    }
    return nullptr;
  }
  Node* AppendChild(const std::string& name_) {
    children.emplace_back(name_, 0);
    return &children.back();
  }
  std::vector<Node> children;
  std::string name;
  double time;
};

void PrintTree(
    const Node& root, std::ostream& out, std::string prefix, double timeall) {
  out << util::Format(
      "{}{} [{:.3f} s, {:.2f}%]\n", prefix, root.name, root.time,
      root.time * 100 / timeall);
  for (const auto& node : root.children) {
    PrintTree(node, out, prefix + "|     ", timeall);
  }
}
void PrintTree(const Node& root, std::ostream& out) {
  PrintTree(root, out, "", root.time);
}
void AccumulateTimeFromLeaves(Node& root) {
  if (root.children.empty()) {
    return;
  }
  double timesum = 0;
  for (auto& node : root.children) {
    AccumulateTimeFromLeaves(node);
    timesum += node.time;
  }
  root.time = timesum;
}

} // namespace suspender_tree

void ParseReport(
    const std::map<std::string, double>& timings, std::ostream& out) {
  using namespace suspender_tree;

  std::vector<std::string> names_prev; // list of strings previous
  Node root("all", 0);
  for (auto& pair : timings) {
    std::vector<std::string> names = SplitBySeparator(pair.first, " --> ");
    if (names.empty() || (names.size() == 1 && names[0] == "")) {
      names = {"other"};
    }

    // collect names from root to last common node with names_prev
    size_t i = 0;
    Node* pos = &root;
    while (i < names.size() && i < names_prev.size() &&
           names[i] == names_prev[i]) {
      pos = pos->FindChild(names[i]);
      ++i;
    }
    while (i < names.size()) {
      pos = pos->AppendChild(names[i]);
      ++i;
    }
    pos->time = pair.second;
    names_prev.swap(names);
  }

  AccumulateTimeFromLeaves(root);
  PrintTree(root, out);
}
