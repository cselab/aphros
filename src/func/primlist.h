#pragma once

#include <cmath>
#include <stdexcept>
#include <functional>
#include <limits>
#include <sstream>
#include <vector>

template <class Scal>
struct UPrimList {
  static constexpr size_t dim = 3;
  using Vect = GVect<Scal, dim>;

  struct Primitive {
    Vect c;  // center
    Vect r;  // axes in coordinate directions
  };

  static std::vector<Primitive> Parse(std::string fn, bool verb) {
    std::vector<Primitive> pp;
    std::ifstream f(fn);
    if (!f.good() && verb) {
      throw std::runtime_error("Can't open particle list '" + fn + "'");
    }

    f >> std::skipws;
    // Read until eof
    while (true) {
      Primitive p;
      // Read single particle: x y z r
      // first character (to skip empty strings)
      char c;
      f >> c;
      if (f.good()) {
        std::string s;
        std::getline(f, s);
        if (c == '#') {
          continue;
        }
        s = c + s; // append first character
        std::stringstream st(s);
        st >> p.c[0] >> p.c[1] >> p.c[2];
        st >> p.r[0];
        if (st.fail()) {
          throw std::runtime_error("list: missing rx in '" + s + "'");
        }
        st >> p.r[1];
        if (st.fail()) {
          p.r[1] = p.r[0];
          p.r[2] = p.r[0];
        } else {
          st >> p.r[2];
          if (st.fail()) {
            p.r[2] = p.r[0];
          }
        }
        pp.push_back(p);
      } else {
        break;
      }
    }

    if (verb) {
      std::cout << "Read " << pp.size() << " particles from " 
          << "'" << fn << "'" << std::endl;
    }
    return pp;
  }
};
