#undef NDEBUG
#include <sstream>
#include <iostream>
#include <cassert>
#include <functional>
#include <cmath>

#include "geom/mesh.h"
#include "solver/solver.h"

const int dim = 3;
using MIdx = GMIdx<dim>;
using IdxCell = IdxCell;
using IdxFace = IdxFace;
using Dir = GDir<dim>;
using Scal = double;
using Vect = GVect<Scal, dim>;
using Mesh = MeshStructured<Scal, dim>;

bool Cmp(Scal a, Scal b) {
  return std::abs(a - b) < 1e-12;
}

template <class T>
bool Cmp(T* a, T* b) {
  return a == b;
}

bool Cmp(size_t a, size_t b) {
  return a == b;
}

template <class V>
typename V::value_type Prod(const V& v) {
  auto r = v[0];
  for (size_t i = 1; i < v.size(); ++i) {
    r *= v[i];
  }
  return r;
}

template <class Idx, class M>
typename Mesh::Scal DiffMax(
    const GField<typename Mesh::Scal, Idx>& u,
    const GField<typename Mesh::Scal, Idx>& v,
    const M& m) {
  using Scal = typename Mesh::Scal;
  Scal r = 0;
  for (auto i : m.template Get<Idx>()) {
    r = std::max(r, std::abs(u[i] - v[i]));
  }
  return r;
}

template <class Idx, class M>
typename Mesh::Scal DiffMax(
    const GField<typename M::Vect, Idx>& u,
    const GField<typename M::Vect, Idx>& v,
    const M& m) {
  using Scal = typename Mesh::Scal;
  Scal r = 0;
  for (auto i : m.template Get<Idx>()) {
    r = std::max(r, (u[i] - v[i]).norminf());
  }
  return r;
}


#define CMP(a, b) \
  assert(Cmp(a, b)); 

// Print CMP
#define PCMP(a, b) \
  std::cerr << #a << "=" << a << ", " << #b << "=" << b << std::endl; \
  CMP(a, b); 

// Echo Execute
#define EE(...); std::cerr << "\n" << #__VA_ARGS__ << std::endl; __VA_ARGS__;


Mesh GetMesh(MIdx s /*size in cells*/) {
  Rect<Vect> dom(Vect(0.1, 0.2, 0.1), Vect(1.1, 1.2, 1.3));
  MIdx b(-2, -3, -4); // lower index
  int hl = 1;         // halos 
  Vect doms = dom.GetDimensions();
  Vect h = dom.GetDimensions() / Vect(s);
  return InitUniformMesh<Mesh>(dom, b, s, hl, true);
}

template <class T, class Idx, class M>
void Eval(std::function<T(Vect)> f,   
          GField<T, Idx>& r,
          const M& m) {
  r.Reinit(m);
  for (auto i : m.template GetAll<Idx>()) {
    r[i] = f(m.GetCenter(i));
  }
}

template <class M>
Scal RunInterp(std::function<Scal(Vect)> uf, const M& m) {
  // Init field on all cells (including halos)
  FieldCell<Scal> cf;
  Eval(uf, cf, m);

  // Empty bc
  MapFace<std::shared_ptr<solver::CondFace>> bc; 

  // Interpolate to faces
  FieldFace<Scal> ff = solver::Interpolate(cf, bc, m);

  // Init reference on faces
  FieldFace<Scal> fr;  
  Eval(uf, fr, m);

  return DiffMax(ff, fr, m);
}

Scal RunGrad(std::function<Scal(Vect)> uf, 
             std::function<Vect(Vect)> ugr,  // reference
             const Mesh& m) {
  // Init field on all cells (including halos)
  FieldCell<Scal> f;
  Eval(uf, f, m);

  // Empty bc
  MapFace<std::shared_ptr<solver::CondFace>> bc; 

  // Interpolate to faces
  FieldFace<Scal> ff = solver::Interpolate(f, bc, m);

  // Gradient on cells
  FieldCell<Vect> g = solver::Gradient(ff, m);

  // Init reference on cells
  FieldCell<Vect> gr;  
  Eval(ugr, gr, m);

  return DiffMax(g, gr, m);
}

void Single(std::function<Scal(Vect)> f) {
  auto m = GetMesh(MIdx(5, 4, 3));
  auto e = RunInterp(f, m);
  std::cerr << std::scientific
      << "err: " << e
      << std::endl;
  CMP(e, 0.);
}

void SingleG(std::function<Scal(Vect)> f, std::function<Vect(Vect)> gr) {
  auto m = GetMesh(MIdx(5, 4, 3));
  auto e = RunGrad(f, gr, m);
  std::cerr << std::scientific
      << "err: " << e
      << std::endl;
  CMP(e, 0.);
}

void Conv(std::function<Scal(Vect)> f) {
  MIdx s(2, 3, 3); // initial mesh size
  Scal rh = 2; // refinement factor (for number of cells in one direction)

  Scal ep = -1;
  Scal ord = 0.;
  size_t n = 0;
  while (s.prod() < 1e4) {
    auto m = GetMesh(s);
    Scal e = RunInterp(f, m); 
    // e = C * h ^ ord
    // ep / e = rh ^ ord
    if (ep > 0) {
      ord = std::log(ep / e) / std::log(rh);
    }
    std::cerr 
        << "s=" << s 
        << "\terr=" << e 
        << "\tord=" << ord
        << std::endl;
    s = MIdx(Vect(s) * rh);
    ep = e;
    ++n;
  }
  assert(n > 1);
  assert(ord > 1.8 || ep < 1e-12);
}

void ConvG(std::function<Scal(Vect)> f, std::function<Vect(Vect)> r) {
  MIdx s(2, 3, 3); // initial mesh size
  Scal rh = 2; // refinement factor (for number of cells in one direction)

  Scal ep = -1;
  Scal ord = 0.;
  size_t n = 0;
  while (s.prod() < 1e4) {
    auto m = GetMesh(s);
    Scal e = RunGrad(f, r, m); 
    // e = C * h ^ ord
    // ep / e = rh ^ ord
    if (ep > 0) {
      ord = std::log(ep / e) / std::log(rh);
    }
    std::cerr 
        << "s=" << s 
        << "\terr=" << e 
        << "\tord=" << ord
        << std::endl;
    s = MIdx(Vect(s) * rh);
    ep = e;
    ++n;
  }
  assert(n > 1);
  assert(ord > 1.8 || ep < 1e-12);
}

Scal sqr(Scal a) {
  return a * a;
}

Scal cube(Scal a) {
  return a * a * a;
}

template <class T>
std::ostream& operator<<(std::ostream& o, const std::vector<T>& v) {
  char c = '(';
  for (auto& a : v) {
    o << c << a;
    c = ',';
  }
  o << ")";
  return o;
}

void TestGrad() {
  EE(Single([](Vect x) { return 0.; }));
  EE(Single([](Vect x) { return 1.; }));
  EE(Single([](Vect x) { return x[0]; }));
  EE(Single([](Vect x) { return x[0] * x[1] * x[2]; }));
  EE(Conv([](Vect x) { return std::sin(x[0] + sqr(x[1]) + cube(x[2])); }));

  EE(SingleG([](Vect x) { return 1.; },
             [](Vect x) { return Vect(0., 0., 0.); }));
  EE(SingleG([](Vect x) { return x[0]; },
             [](Vect x) { return Vect(1., 0., 0.); }));
  EE(SingleG([](Vect x) { return sqr(x[0]); },
             [](Vect x) { return Vect(2 * x[0], 0., 0.); }));
  EE(SingleG([](Vect x) { return sqr(x[0]); },
             [](Vect x) { return Vect(2 * x[0], 0., 0.); }));
  EE(SingleG([](Vect x) { return sqr(x[0]) * sqr(x[1]) * sqr(x[2]); },
             [](Vect x) { return Vect(
                 2 * sqr(x[1]) * sqr(x[2]) * x[0], 
                 2 * sqr(x[2]) * sqr(x[0]) * x[1], 
                 2 * sqr(x[0]) * sqr(x[1]) * x[2]
                 ); }));

  EE(ConvG([](Vect x) { return std::sin(x[0] + sqr(x[1]) + cube(x[2])); },
             [](Vect x) { return Vect(
                 std::cos(x[0] + sqr(x[1]) + cube(x[2])),
                 std::cos(x[0] + sqr(x[1]) + cube(x[2])) * 2. * x[1],
                 std::cos(x[0] + sqr(x[1]) + cube(x[2])) * 3. * sqr(x[2])
                 ); }));
}

void TestCoeff() {
  auto p = [](Scal x, const std::vector<Scal>& z, 
              const std::vector<Scal>& ke, size_t b=0) {
    std::vector<Scal> k;
    if (b == 0) {
      k = solver::GetGradCoeffs(x, z);
    } else {
      k = solver::GetGradCoeffs(x, z, b);
    }
    std::cerr
        << "x=" << x
        << " z=" << z
        << " k=" << k
        << " b=" << b
        << " ke=" << ke 
        << std::endl;
    for (size_t i = 0; i < k.size(); ++i) {
      assert(std::abs(k[i] - ke[i]) < 1e-12);
    }
  };

  p(0, {-1, 0}, {-1, 1});
  p(0, {-2, -1, 0}, {0.5, -2., 1.5});
  p(0, {-2, -1, 0}, {0, -1, 1}, 1);
  p(0, {-1, 0, 1}, {-0.5, 0, 0.5});
}

int main() {
  TestGrad();
  TestCoeff();
}
