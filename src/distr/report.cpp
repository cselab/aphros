#include <sstream>
#include <vector>
#include <iomanip>

#include "report.h"

// Splits string by whitespace, skips every second.
std::vector<std::string> Split(const std::string& str) {
  std::vector<std::string> ss; // result
  std::stringstream st(str);
  st >> std::skipws;
  while (st) {
    std::string s;
    st >> s;
    ss.push_back(s);
    st >> s;  // skip 
  }
  return ss;
};

// Node
struct N { 
  const double kNone = -1.; // empty time
  std::vector<N> nn; // children
  std::string s; // stage name
  double t; // time
  N(const std::string& s, double t) : s(s), t(t) {}
  N(const std::string& s) : s(s), t(kNone) {}
  N* Find(const std::string& s) {
    for (auto& n : nn) {
      if (n.s == s) {
        return &n;
      }
    }
    return nullptr;
  }
  N* Add(const std::string& s) {
    nn.emplace_back(s, kNone);
    return &nn.back();
  }
  void Print(std::ostream& out, std::string pre, double ta) const {
    auto fl = out.flags();
    out << pre << s << " [" 
        << std::setprecision(5) << t << " s, "
        << std::setprecision(3) << 100. * t / ta << "%]" << std::endl;
    for (auto& n : nn) {
      n.Print(out, pre + "|     ", ta);
    }
    out.flags(fl);
  }
  void FillTime() {
    double ts = 0.; // sum
    if (t == kNone) {
      for (auto& n : nn) {
        if (n.t == kNone) {
          n.FillTime();
        } 
        ts += n.t;
      }
      t = ts;
    }
  }
};

void ParseReport(const std::map<std::string, double>& mp, std::ostream& out) {
  std::vector<std::string> ss0; // list of strings previous
  N r("all"); // root
  // split strings, mp stores graph traversal in depth-first order
  for (auto it : mp) {
    std::vector<std::string> ss = Split(it.first);
    if (ss.empty() || (ss.size() == 1 && ss[0] == "")) {
      ss = {"other"};
    }

    // from root to last common node with ss0
    size_t i = 0;
    N* p = &r; // position
    while (i < ss.size() && i < ss0.size() && ss[i] == ss0[i]) {
      p = p->Find(ss[i]);
      ++i;
    }
    while (i < ss.size()) {
      p = p->Add(ss[i]);
      ++i;
    }
    p->t = it.second;
    ss0 = ss;
  }

  r.FillTime();
  r.Print(out, "", r.t);
}
