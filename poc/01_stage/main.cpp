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
#include <vector>

using Scal = double;

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
      if (i->c == i->t) {
        fprintf(stderr, "[%s] operator() %s\n", name_.c_str(), m.LscStr().c_str());
      }
      return i->c++ == i->t;
    }
  };
  friend Stage;
  explicit Mesh(size_t n) 
  : lsc_(1, SC(-1,-1)), lsci_(lsc_.begin()), u(n)
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
  void SetName(std::string name) {
    name_ = name;
  }
  std::string GetName() const {
    return name_;
  }
  void Comm(std::vector<Scal>& u) {
    cm_.push_back(&u);
  }
  std::vector<Scal> u;
  std::vector<std::vector<Scal>*> cm_;

 private:
  using LSC = std::list<SC>;
  mutable LSC lsc_;
  mutable LSC::iterator lsci_;
  std::string name_;
};

void CommCommit(std::vector<Mesh>& mm) {
  const size_t n = mm.size();
  for (size_t i = 0; i < n; ++i) {
    for (size_t k = 0; k < mm[i].cm_.size(); ++k) {
      const size_t il = (i + n - 1) % n;
      const size_t ir = (i + 1) % n;
      assert(mm[il].cm_[k]->size() == mm[i].cm_[k]->size());
      assert(mm[ir].cm_[k]->size() == mm[i].cm_[k]->size());

      fprintf(stderr, "comm=%d src=%d dstl=%d dstr=%d\n",
          k, i, il, ir);
      (*mm[il].cm_[k]).back() = (*mm[i].cm_[k]).front();
      (*mm[ir].cm_[k]).front() = (*mm[i].cm_[k]).back();
    }
  }
}

void Add(const Mesh& m) {
  auto st = m.GetStage(m.GetName() + "_add");
  if (st()) {
  } 
  if (st()) {
  } 
}

void Grad(Mesh& m) {
  auto st = m.GetStage(m.GetName() + "_grad");
  if (st()) {
    //Add(m);
    m.Comm(m.u);
  } 
  if (st()) {
    //Add(m);
  } 
}

int main() {
  const int n = 2;
  const int nu = 10;
  // init mesh
  std::vector<Mesh> mm;
  for (size_t i = 0; i < n; ++i) {
    mm.emplace_back(nu);
    mm.back().SetName("mesh" + std::to_string(i));
  }

  // comp and comm
  do {
    for (Mesh& m : mm) {
      Grad(m);
    }
    for (Mesh& m : mm) {
      assert(mm[0].CallPending() == m.CallPending());
    }
    CommCommit(mm);
  } while (mm[0].CallPending());

  return 0;
}
