#include <iostream>
#include <omp.h>
#include <cmath>
#include "unistd.h"
#include <iomanip>
#include <omp.h>
#include <array>
#include <cassert>
#include <sstream>
#include <fstream>


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

const size_t wh = 1; // width halo

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
        //fprintf(stderr, "[%s] operator() %s\n", name_.c_str(), m.LscStr().c_str());
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

bool CallPending(const std::vector<Mesh>& mm) {
  if (mm.empty()) {
    return false;
  }
  for (const Mesh& m : mm) {
    assert(mm[0].CallPending() == m.CallPending());
  }
  return mm[0].CallPending();
}

void CommCommit(std::vector<Mesh>& mm) {
  const size_t n = mm.size();
  for (size_t i = 0; i < n; ++i) {
    const size_t il = (i + n - 1) % n;
    const size_t ir = (i + 1) % n;
    assert(mm[il].cm_.size() == mm[i].cm_.size());
    assert(mm[ir].cm_.size() == mm[i].cm_.size());
    //fprintf(stderr, "src=%d dstl=%d dstr=%d\n", i, il, ir);
    for (size_t k = 0; k < mm[i].cm_.size(); ++k) {
      auto& u = *mm[i].cm_[k];
      auto& ul = *mm[il].cm_[k];
      auto& ur = *mm[ir].cm_[k];
      for (size_t h = 0; h < wh; ++h) {
        ul[ul.size() - wh + h] = u[wh + h];
        ur[h] = u[u.size() - wh + h - wh];
      }
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

void Adv(Mesh& m, bool rec=true) {
  auto st = m.GetStage(m.GetName() + "_adv");
  auto& u = m.u;
  const Scal cfl = 0.3;
  const size_t n = u.size();
  if (st()) {
    for (size_t i = wh; i < n-wh; ++i) {
      const size_t il = i - 1;
      const size_t ir = i + 1;
      u[i] += -cfl * (u[i] - u[il]);
    }
    m.Comm(m.u);
  } 
  if (1 && st()) {
    for (size_t i = wh; i < n-wh; ++i) {
      const size_t il = i - 1;
      const size_t ir = i + 1;
      u[i] += -cfl * (u[i] - u[il]);
    }
    m.Comm(m.u);
  } 
  if (rec && st()) {
    Adv(m, false);
  } 
}

//TODO no_nested_stages flag
//     or allow_nested_stages flag to Stage

void Steps(Mesh& m, size_t nt, bool root) {
  auto st = m.GetStage(m.GetName() + "_steps");
  for (size_t t = 0; t < nt; ++t) {
    if (st()) {
      root && fprintf(stderr, "t=%03d\n", t);
    }
    if (st()) {
      Adv(m);
    }
  }
}

size_t t;
void StepsGlbt(Mesh& m, size_t nt, bool root) {
  auto st = m.GetStage(m.GetName() + "_steps");
  for (; t < nt; ) {
    if (st()) {
      root && fprintf(stderr, "t=%03d\n", t);
      root && ++t;
      break;
    }
    if (st()) {
      Adv(m);
      break;
    }

  }
}

void Init(Mesh& m, size_t mi, size_t mn) {
  auto st = m.GetStage(m.GetName() + "_init");
  if (st()) {
    const size_t nu = m.u.size();
    for (size_t j = 0; j < nu-2*wh; ++j) {
      const size_t gnu = (nu - 2 * wh) * mn;
      const size_t gj = mi * (nu - 2*wh) + j;
      const Scal x = double(gj) / (gnu - 1);
      m.u[j+wh] = sin(x * 10) * x * (1-x) * sin(x * 5);
    }
    m.Comm(m.u);
  }
}

void Write(Mesh& m, std::ofstream& f) {
  const size_t nu = m.u.size();
  for (size_t j = wh; j < nu-wh; ++j) {
    f << m.u[j] << "\n";
  }
}

void WriteAll(std::vector<Mesh>& mm, std::string fn) {
  std::ofstream f(fn);
  for (Mesh& m : mm) {
    Write(m, f);
  }
}

int main() {
  const size_t n = 10;
  const size_t nu = 50 / n + 2 * wh;
  // init mesh
  std::vector<Mesh> mm;
  for (size_t i = 0; i < n; ++i) {
    mm.emplace_back(nu);
    auto& m = mm.back();
    m.SetName("mesh" + std::to_string(i));
  }

  // init
  do {
    for (size_t i = 0; i < n; ++i) {
      Init(mm[i], i, n);
    }
    CommCommit(mm);
  } while (CallPending(mm));

  WriteAll(mm, "a.dat");

  // comp and comm
  const size_t nt = 10;
  if (0) {
    for (size_t t = 0; t < nt; ++t) {
      fprintf(stderr, "t=%03d\n", t);
      do {
        for (Mesh& m : mm) {
          Adv(m);
        }
        CommCommit(mm);
      } while (CallPending(mm));
    }
  }
  if (1) {
    do {
      for (size_t i = 0; i < n; ++i) {
        Steps(mm[i], nt, i == 0);
      }
      CommCommit(mm);
    } while (CallPending(mm));
  }
  if (0) {
    t = 0;
    do {
      for (size_t i = 0; i < n; ++i) {
        StepsGlbt(mm[i], nt, i == 0);
      }
      CommCommit(mm);
    } while (CallPending(mm));
  }

  WriteAll(mm, "b.dat");

  return 0;
}
