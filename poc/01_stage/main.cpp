#include <iostream>
#include <omp.h>
#include <cmath>
#include "unistd.h"
#include <iomanip>
#include <omp.h>
#include <array>
#include <cassert>
#include <sstream>


void test() {
  #pragma omp parallel 
  {
    #pragma omp critical
    std::cout << "hello: " << omp_get_thread_num() << std::endl;
    #pragma omp single
    std::cout << "hello: " << omp_get_thread_num() << std::endl;
    #pragma omp master
    std::cout << "hello: " << omp_get_thread_num() << std::endl;
  }
}



#define N 1000
float F(float); 
float x[N+2], u[N+2], p;

void First() {
  for (int n = 1; n <= N; ++n) {
    float xn = x[n];
    float un = u[n];
    float xnp = xn + un;
    float f = F(xnp);
    x[n+1] = xnp;
    u[n+1] = un + f;
    p += u[n-1] * f;
  }
}

void Second() {
  for (int n = 1; n <= N; ++n) {
    float xn = x[n];
    float un = u[n];
    float xnp = xn + un;
    float f = F(xnp);
    x[n+1] = xnp;
    u[n+1] = un + f;
    p += (un - F(xn)) * f;
  }
}

float F(float x) {
  return x;
}

#include <list>
#include <stdexcept>
class Mesh {
 public:
  struct SC { // stage counter
    int c; // current
    int t; // target
    SC(int c, int t) : c(c), t(t) {}
  };
  class Stage {
    const Mesh& m;
    std::string name_;
   public:
    Stage(const Mesh& m, std::string name="") 
    : m(m), name_(name)
    {
      //fprintf(stderr, "[%s] Stage()\n", name_.c_str());
      auto& s = m.lsc_;
      auto& i = m.lsci_;
      if (std::next(i) == s.end()) {
        //fprintf(stderr, "[%s] emplace_back(0, 0)\n", name_.c_str());
        s.emplace_back(0, 0);
      }
      ++i;
      i->c = 0;
    }
    ~Stage() {
      //fprintf(stderr, "[%s] ~Stage()\n", name_.c_str());
      auto& s = m.lsc_;
      auto& i = m.lsci_;
      assert(!s.empty());
      assert(i != s.end());

      assert(i != s.begin());
      auto ip = std::prev(i);

      if (std::next(i) == s.end()) {
        // all lower calls done, next stage
        //fprintf(stderr, "[%s] all lower calls done\n", name_.c_str());
        ++i->t;
        if (i->c == i->t) { 
          // all stages done, remove current call
          // i->c is number of stages
          //fprintf(stderr, "[%s] all stages done, remove current call\n",
          //    name_.c_str());
          s.pop_back();
        }
      } else {
        //fprintf(stderr, "[%s] keep stage\n", name_.c_str());
      }
      i = ip;
      //fprintf(stderr, "[%s] %s\n", name_.c_str(), m.LscStr().c_str());
    }
    bool operator()() {
      auto& i = m.lsci_;
      //fprintf(stderr, "[%s] operator() c=%d, t=%d\n", name_.c_str(), i->c, i->t);
      return i->c++ == i->t;
    }
  };
  friend Stage;
  Mesh() 
  : lsc_(1, SC(-1,-1)), lsci_(lsc_.begin())
  {}
  Stage GetStage(std::string name="") const {
    return Stage(*this, name);
  }
  std::string LscStr() const {
    std::stringstream b;
    for (auto sc : lsc_) {
      b << "(" << sc.c << " " << sc.t << ") ";
    }
    return b.str();
  }
  bool CallPending() const {
    return lsc_.size() != 1;
  }

 private:
  using LSC = std::list<SC>;
  mutable LSC lsc_;
  mutable LSC::iterator lsci_;
};

void Add(const Mesh& m) {
  auto st = m.GetStage("add");
  if (st()) {
    std::cerr << "Add stage0" << std::endl;
  } 
  if (st()) {
    std::cerr << "Add stage1" << std::endl;
  } 
}



void Grad(const Mesh& m) {
  auto st = m.GetStage("grad");
  if (st()) {
    std::cerr << "Grid stage0" << std::endl;
    Add(m);
  } 
  if (st()) {
    std::cerr << "Grid stage1" << std::endl;
    Add(m);
  } 
}

int main() {
  Mesh m;

  int i = 0;
  do {
    std::cerr << "Grad(m) i=" << i << std::endl;
    Grad(m);
    ++i;
  } while (m.CallPending());

  return 0;
}
