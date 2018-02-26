#include <sstream>
#include <iostream>
#include <cassert>
#include <functional>
#include <cmath>

#include "mesh.hpp"
#include "mesh3d.hpp"
#include "solver.hpp"

const int dim = 3;
using MIdx = geom::GMIdx<dim>;
using IdxCell = geom::IdxCell;
using IdxFace = geom::IdxFace;
using Dir = geom::GDir<dim>;
using Scal = double;
using Vect = geom::GVect<Scal, dim>;
using Mesh = geom::MeshStructured<Scal, dim>;

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

template <class Idx, class Mesh>
typename Mesh::Scal Max(
    const geom::GField<typename Mesh::Scal, Idx>& u,
    const Mesh& m) {
  using Scal = typename Mesh::Scal;
  Scal r = 0;
  for (auto i : m.template Get<Idx>()) {
    r = std::max(r, u[i]);
  }
  return r;
}

template <class Idx, class Mesh>
typename Mesh::Scal DiffMax(
    const geom::GField<typename Mesh::Scal, Idx>& u,
    const geom::GField<typename Mesh::Scal, Idx>& v,
    const Mesh& m) {
  using Scal = typename Mesh::Scal;
  Scal r = 0;
  for (auto i : m.template Get<Idx>()) {
    r = std::max(r, std::abs(u[i] - v[i]));
  }
  return r;
}

#define CMP(a, b) \
  assert(Cmp(a, b)); 

// Print CMP
#define PCMP(a, b) \
  std::cerr << #a << "=" << a << ", " << #b << "=" << b << std::endl; \
  CMP(a, b); 

// EE - Echo Execute
#define EE(...); std::cerr << #__VA_ARGS__ << std::endl; __VA_ARGS__;


Mesh GetMesh(MIdx s /*size in cells*/) {
  geom::Rect<Vect> dom(Vect(0.1, 0.2, 0.1), Vect(1.1, 1.2, 1.3));
  MIdx b(-2, -3, -4); // lower index
  int hl = 1;         // halos 
  Vect doms = dom.GetDimensions();
  Vect h = dom.GetDimensions() / Vect(s);
  return geom::InitUniformMesh<Mesh>(dom, b, s, hl);
}

Scal TestInterp0(std::function<Scal(Vect)> f, const Mesh& m) {
  // Init field on all cells (including halos)
  geom::FieldCell<Scal> fc(m);
  for (auto i : m.AllCells()) {
    Vect x = m.GetCenter(i);
    fc[i] = f(x);
  }

  // Empty bc
  geom::MapFace<std::shared_ptr<solver::ConditionFace>> mfc; 

  // Inteatpolate to faces
  geom::FieldFace<Scal> ff = solver::Interpolate(fc, mfc, m);

  // Init reference on faces
  geom::FieldFace<Scal> ffr(m);  
  for (auto i : m.AllFaces()) {
    Vect x = m.GetCenter(i);
    ffr[i] = f(x);
  }
  return DiffMax(ff, ffr, m);
}

void Single(std::function<Scal(Vect)> f) {
  auto m = GetMesh(MIdx(5, 4, 3));
  auto e = TestInterp0(f, m);
  std::cerr << std::scientific
      << "err: " << e
      << std::endl;
  CMP(e, 0.);
}

void Conv(std::function<Scal(Vect)> f) {
  MIdx s(2, 3, 4); // initial mesh size
  Scal rh = 2; // refinement factor (for number of cells in one direction)

  Scal ep = -1;
  Scal ord = 0.;
  size_t n = 0;
  while (Prod(s) < 1e5) {
    auto m = GetMesh(s);
    Scal e = TestInterp0(f, m); 
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

int main() {
  EE(Single([](Vect x) { return 0.; }));
  EE(Single([](Vect x) { return 1.; }));
  EE(Single([](Vect x) { return x[0]; }));
  EE(Single([](Vect x) { return x[0] * x[1] * x[2]; }));
  EE(Conv([](Vect x) { return std::sin(x[0] + sqr(x[1]) + cube(x[2])); }));
}
