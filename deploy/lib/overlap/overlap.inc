#include <Eigen/Core>
#include <Eigen/Geometry>

#include <algorithm>
#include <array>
#include <bitset>
#include <cassert>
#include <cmath>
#include <iterator>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <type_traits>
#include <utility>

// typedefs
typedef double scalar_t;
typedef Eigen::Matrix<scalar_t, 3, 1, Eigen::DontAlign> vector_t;

// constants
const scalar_t pi = scalar_t(4) * std::atan(scalar_t(1.0));

namespace detail {

static const scalar_t tinyEpsilon(2 * std::numeric_limits<scalar_t>::epsilon());

static const scalar_t mediumEpsilon(1e2 * tinyEpsilon);
static const scalar_t largeEpsilon(1e-10);

// Robust calculation of the normal vector of a polygon using Newell's method
// and a pre-calculated center.
// Ref: Christer Ericson - Real-Time Collision Detection (2005)
template <typename Iterator>
inline vector_t normalNewell(
    Iterator begin, Iterator end, const vector_t& center) {
  const size_t count = end - begin;
  vector_t n(vector_t::Zero());

  for (size_t i = 0; i < count; ++i)
    n += (*(begin + i) - center).cross(*(begin + ((i + 1) % count)) - center);

  scalar_t length = n.stableNorm();

  if (length)
    return n / length;
  else
    return n;
}

// This implementation of double_prec is based on:
// T.J. Dekker, A floating-point technique for extending the available
// precision, http://dx.doi.org/10.1007/BF01397083

template <typename T>
struct double_prec_constant;

template <>
struct double_prec_constant<float> {
  // Constant used to split double precision values:
  // 2^(24 - 24/2) + 1 = 2^12 + 1 = 4097
  static const uint32_t value = 4097;
};

template <>
struct double_prec_constant<double> {
  // Constant used to split double precision values:
  // 2^(53 - int(53/2)) + 1 = 2^27 + 1 = 134217729
  static const uint32_t value = 134217729;
};

// For GCC and Clang an attribute can be used to control the FP precision...
#if defined(__GNUC__) && !defined(__clang__) && !defined(__ICC) && \
    !defined(__INTEL_COMPILER)
#define ENFORCE_EXACT_FPMATH_ATTR __attribute__((__target__("ieee-fp")))
#else
#define ENFORCE_EXACT_FPMATH_ATTR
#endif

// ... whereas ICC requires a pragma.
#if defined(__ICC) || defined(__INTEL_COMPILER)
#define ENFORCE_EXACT_FPMATH_ATTR
#define USE_EXACT_FPMATH_PRAGMA 1
#endif

template <typename T>
class double_prec;

template <typename T>
inline double_prec<T> operator+(
    const double_prec<T>& lhs,
    const double_prec<T>& rhs) ENFORCE_EXACT_FPMATH_ATTR;

template <typename T>
inline double_prec<T> operator-(
    const double_prec<T>& lhs,
    const double_prec<T>& rhs) ENFORCE_EXACT_FPMATH_ATTR;

template <typename T>
inline double_prec<T> operator*(
    const double_prec<T>& lhs,
    const double_prec<T>& rhs)ENFORCE_EXACT_FPMATH_ATTR;

template <typename T>
class double_prec {
 private:
  static const uint32_t c = detail::double_prec_constant<T>::value;

  template <typename TF>
  friend double_prec<TF> operator+(
      const double_prec<TF>&, const double_prec<TF>&);

  template <typename TF>
  friend double_prec<TF> operator-(
      const double_prec<TF>&, const double_prec<TF>&);

  template <typename TF>
  friend double_prec<TF> operator*(
      const double_prec<TF>&, const double_prec<TF>&);

 public:
  inline double_prec() : h_(0), l_(0) {}

  // This constructor requires floating point operations in accordance
  // with IEEE754 to perform the proper splitting. To allow full
  // optimization of all other parts of the code, precise floating point
  // ops are only requested here. Unfortunately the way to do this is
  // extremely compiler dependent.
  inline double_prec(const T& val) ENFORCE_EXACT_FPMATH_ATTR : h_(0), l_(0) {
#ifdef USE_EXACT_FPMATH_PRAGMA
#pragma float_control(precise, on)
#endif

    T p = val * T(c);
    h_ = (val - p) + p;
    l_ = val - h_;
  }

 private:
  inline explicit double_prec(const T& h, const T& l) : h_(h), l_(l) {}

 public:
  inline const T& high() const {
    return h_;
  }

  inline const T& low() const {
    return l_;
  }

  inline T value() const {
    return h_ + l_;
  }

  template <typename TOther>
  inline TOther convert() const {
    return TOther(h_) + TOther(l_);
  }

 private:
  T h_;
  T l_;
};

template <typename T>
inline double_prec<T> operator+(
    const double_prec<T>& lhs, const double_prec<T>& rhs) {
#ifdef USE_EXACT_FPMATH_PRAGMA
#pragma float_control(precise, on)
#endif

  T h = lhs.h_ + rhs.h_;
  T l = std::abs(lhs.h_) >= std::abs(rhs.h_)
            ? ((((lhs.h_ - h) + rhs.h_) + lhs.l_) + rhs.l_)
            : ((((rhs.h_ - h) + lhs.h_) + rhs.l_) + lhs.l_);

  T c = h + l;

  return double_prec<T>(c, (h - c) + l);
}

template <typename T>
inline double_prec<T> operator-(
    const double_prec<T>& lhs, const double_prec<T>& rhs) {
#ifdef USE_EXACT_FPMATH_PRAGMA
#pragma float_control(precise, on)
#endif

  T h = lhs.h_ - rhs.h_;
  T l = std::abs(lhs.h_) >= std::abs(rhs.h_)
            ? ((((lhs.h_ - h) - rhs.h_) - rhs.l_) + lhs.l_)
            : ((((-rhs.h_ - h) + lhs.h_) + lhs.l_) - rhs.l_);

  T c = h + l;

  return double_prec<T>(c, (h - c) + l);
}

template <typename T>
inline double_prec<T> operator*(
    const double_prec<T>& lhs, const double_prec<T>& rhs) {
#ifdef USE_EXACT_FPMATH_PRAGMA
#pragma float_control(precise, on)
#endif

  double_prec<T> l(lhs.h_);
  double_prec<T> r(rhs.h_);

  T p = l.h_ * r.h_;
  T q = l.h_ * r.l_ + l.l_ * r.h_;
  T v = p + q;

  double_prec<T> c(v, ((p - v) + q) + l.l_ * r.l_);
  c.l_ = ((lhs.h_ + lhs.l_) * rhs.l_ + lhs.l_ * rhs.h_) + c.l_;
  T z = c.value();

  return double_prec<T>(z, (c.h_ - z) + c.l_);
}

// Ref: J.R. Shewchuk - Lecture Notes on Geometric Robustness
//      http://www.cs.berkeley.edu/~jrs/meshpapers/robnotes.pdf
inline scalar_t orient2D(
       scalar_t A0, scalar_t A1,
       scalar_t B0, scalar_t B1,
       scalar_t C0, scalar_t C1) {
  typedef double_prec<scalar_t> real_t;

  real_t a0(A0);
  real_t a1(A1);
  real_t b0(B0);
  real_t b1(B1);
  real_t c0(C0);
  real_t c1(C1);

  real_t result = (a0 - c0) * (b1 - c1) - (a1 - c1) * (b0 - c0);

  return result.convert<scalar_t>();
}

// Numerically robust calculation of the normal of the triangle defined by
// the points a, b, and c.
// Ref: J.R. Shewchuk - Lecture Notes on Geometric Robustness
//      http://www.cs.berkeley.edu/~jrs/meshpapers/robnotes.pdf
inline vector_t triangleNormal(
    const vector_t& a, const vector_t& b, const vector_t& c) {
  scalar_t xy = orient2D(a[0], a[1], b[0], b[1], c[0], c[1]);
  scalar_t yz = orient2D(a[1], a[2], b[1], b[2], c[1], c[2]);
  scalar_t zx = orient2D(a[2], a[0], b[2], b[0], c[2], c[0]);
  return vector_t(yz, zx, xy).normalized();
}

// Numerically robust routine to calculate the angle between normalized
// vectors.
// Ref: http://www.plunk.org/~hatch/rightway.php
inline scalar_t angle(const vector_t& v0, const vector_t& v1) {
  if (v0.dot(v1) < scalar_t(0))
    return pi - scalar_t(2) * std::asin(scalar_t(0.5) * (v0 + v1).stableNorm());
  else
    return scalar_t(2) * std::asin(scalar_t(0.5) * (v0 - v1).stableNorm());
}

template <typename Derived0, typename Derived1>
inline std::array<vector_t, 2> gramSchmidt(
    const Eigen::MatrixBase<Derived0>& arg0,
    const Eigen::MatrixBase<Derived1>& arg1) {
  vector_t v0(arg0.normalized());
  vector_t v1(arg1);

  std::array<vector_t, 2> result;
  result[0] = v0;
  result[1] = (v1 - v1.dot(v0) * v0).normalized();

  return result;
}

inline scalar_t clamp(
    scalar_t value, scalar_t min, scalar_t max, scalar_t limit) {
  assert(min <= max && limit >= scalar_t(0));

  value = (value < min && value > (min - limit)) ? min : value;
  value = (value > max && value < (max + limit)) ? max : value;

  return value;
}

} // namespace detail

class Transformation {
 public:
  Transformation(const vector_t& t, const scalar_t& s)
      : translation(t), scaling(s) {}

  vector_t translation;
  scalar_t scaling;
};

template <size_t VertexCount>
class Polygon {
 private:
  static_assert(
      VertexCount >= 3 && VertexCount <= 4,
      "Only triangles and quadrilateral are supported.");

 public:
  static const size_t vertex_count = VertexCount;

 protected:
  Polygon() : vertices(), center(), normal(), area() {}

  template <typename... Types>
  Polygon(const vector_t& v0, Types... verts)
      : vertices{{v0, verts...}}, center(), normal(), area() {
    center = scalar_t(1.0 / vertex_count) *
             std::accumulate(
                 vertices.begin(), vertices.end(), vector_t::Zero().eval());

    // For a quadrilateral, Newell's method can be simplified
    // significantly.
    // Ref: Christer Ericson - Real-Time Collision Detection (2005)
    if (VertexCount == 4) {
      normal = ((vertices[2] - vertices[0]).cross(vertices[3] - vertices[1]))
                   .normalized();
    } else {
      normal = detail::normalNewell(vertices.begin(), vertices.end(), center);
    }
  }

  void apply(const Transformation& t) {
    for (auto& v : vertices)
      v = t.scaling * (v + t.translation);

    center = t.scaling * (center + t.translation);
  }

 public:
  bool isPlanar(const scalar_t epsilon = detail::largeEpsilon) const {
    if (VertexCount == 3) return true;

    for (auto& v : vertices)
      if (std::abs(normal.dot(v - center)) > epsilon) return false;

    return true;
  }

 public:
  std::array<vector_t, vertex_count> vertices;
  vector_t center;
  vector_t normal;
  scalar_t area;
};

class Triangle : public Polygon<3> {
 public:
  Triangle() : Polygon<3>() {}

  template <typename... Types>
  Triangle(const vector_t& v0, Types... verts) : Polygon<3>(v0, verts...) {
    init();
  }

  void apply(const Transformation& t) {
    Polygon<3>::apply(t);
    init();
  }

 private:
  void init() {
    area = scalar_t(0.5) *
           ((vertices[1] - vertices[0]).cross(vertices[2] - vertices[0]))
               .stableNorm();
  }
};

class Quadrilateral : public Polygon<4> {
 public:
  Quadrilateral() : Polygon<4>() {}

  template <typename... Types>
  Quadrilateral(const vector_t& v0, Types... verts) : Polygon<4>(v0, verts...) {
    init();
  }

  void apply(const Transformation& t) {
    Polygon<4>::apply(t);
    init();
  }

 private:
  void init() {
    area = scalar_t(0.5) *
           (((vertices[1] - vertices[0]).cross(vertices[2] - vertices[0]))
                .stableNorm() +
            ((vertices[2] - vertices[0]).cross(vertices[3] - vertices[0]))
                .stableNorm());
  }
};

// Forward declarations of the mesh elements.
class Tetrahedron;
class Wedge;
class Hexahedron;

namespace detail {

// Some tricks are required to keep this code header-only.
template <typename T, typename Nil>
struct mappings;

template <typename Nil>
struct mappings<Tetrahedron, Nil> {
  // Map edges of a tetrahedron to vertices and faces.
  static const uint32_t edge_mapping[6][2][2];

  // Map vertices of a tetrahedron to edges and faces.
  // 0: local IDs of the edges intersecting at this vertex
  // 1: 0 if the edge is pointing away from the vertex, 1 otherwise
  // 2: faces joining at the vertex
  static const uint32_t vertex_mapping[4][3][3];

  // This mapping contains the three sets of the two edges for each of the
  // faces joining at a vertex. The indices are mapped to the local edge IDs
  // using the first value field of the 'vertex_mapping' table.
  static const uint32_t face_mapping[3][2];
};

template <typename Nil>
const uint32_t mappings<Tetrahedron, Nil>::edge_mapping[6][2][2] = {
    {{0, 1}, {0, 1}}, {{1, 2}, {0, 2}}, {{2, 0}, {0, 3}},
    {{0, 3}, {1, 3}}, {{1, 3}, {1, 2}}, {{2, 3}, {2, 3}}};

template <typename Nil>
const uint32_t mappings<Tetrahedron, Nil>::vertex_mapping[4][3][3] = {
    {{0, 2, 3}, {0, 1, 0}, {0, 1, 3}},
    {{0, 1, 4}, {1, 0, 0}, {0, 1, 2}},
    {{1, 2, 5}, {1, 0, 0}, {0, 2, 3}},
    {{3, 4, 5}, {1, 1, 1}, {1, 3, 2}}};

template <typename Nil>
const uint32_t mappings<Tetrahedron, Nil>::face_mapping[3][2] = {
    {0, 1}, {0, 2}, {1, 2}};

typedef mappings<Tetrahedron, void> tet_mappings;

} // namespace detail

class Tetrahedron : public detail::tet_mappings {
 public:
  template <typename... Types>
  Tetrahedron(const vector_t& v0, Types... verts)
      : vertices{{v0, verts...}}, faces(), center(), volume() {
#ifndef NDEBUG
    // Make sure the ordering of the vertices is correct.
    assert(
        (vertices[1] - vertices[0])
            .cross(vertices[2] - vertices[0])
            .dot(vertices[3] - vertices[0]) >= scalar_t(0));
#endif // NDEBUG

    init();
  }

  Tetrahedron(const std::array<vector_t, 4>& verts)
      : vertices(verts), faces(), center(), volume() {
    init();
  }

  Tetrahedron()
      : vertices{{vector_t::Zero(), vector_t::Zero(), vector_t::Zero(),
                  vector_t::Zero()}}
      , faces()
      , center()
      , volume() {}

  void apply(const Transformation& t) {
    for (auto& v : vertices)
      v = t.scaling * (v + t.translation);

    for (auto& f : faces)
      f.apply(t);

    center = scalar_t(0.25) *
             std::accumulate(
                 vertices.begin(), vertices.end(), vector_t::Zero().eval());

    volume = calcVolume();
  }

  scalar_t surfaceArea() const {
    scalar_t area(0);
    for (const auto& f : faces)
      area += f.area;

    return area;
  }

 private:
  void init() {
    // 0: v2, v1, v0
    faces[0] = Triangle(vertices[2], vertices[1], vertices[0]);

    // 1: v0, v1, v3
    faces[1] = Triangle(vertices[0], vertices[1], vertices[3]);

    // 2: v1, v2, v3
    faces[2] = Triangle(vertices[1], vertices[2], vertices[3]);

    // 3: v2, v0, v3
    faces[3] = Triangle(vertices[2], vertices[0], vertices[3]);

    center = scalar_t(0.25) *
             std::accumulate(
                 vertices.begin(), vertices.end(), vector_t::Zero().eval());

    volume = calcVolume();
  }

  scalar_t calcVolume() const {
    return scalar_t(1.0 / 6.0) *
           std::abs((vertices[0] - vertices[3])
                        .dot((vertices[1] - vertices[3])
                                 .cross(vertices[2] - vertices[3])));
  }

 public:
  std::array<vector_t, 4> vertices;
  std::array<Triangle, 4> faces;
  vector_t center;
  scalar_t volume;
};

namespace detail {

template <typename Nil>
struct mappings<Wedge, Nil> {
  // Map edges of a wedge to vertices and faces.
  static const uint32_t edge_mapping[9][2][2];

  // Map vertices of a wedge to edges and faces.
  // 0: local IDs of the edges intersecting at this vertex
  // 1: 0 if the edge is pointing away from the vertex, 1 otherwise
  // 2: faces joining at the vertex
  static const uint32_t vertex_mapping[6][3][3];

  // This mapping contains the three sets of the two edges for each of the
  // faces joining at a vertex. The indices are mapped to the local edge IDs
  // using the first value field of the 'vertex_mapping' table.
  static const uint32_t face_mapping[3][2];
};

template <typename Nil>
const uint32_t mappings<Wedge, Nil>::edge_mapping[9][2][2] = {
    {{0, 1}, {0, 1}}, {{1, 2}, {0, 2}}, {{2, 0}, {0, 3}},
    {{0, 3}, {1, 3}}, {{1, 4}, {1, 2}}, {{2, 5}, {2, 3}},
    {{3, 4}, {1, 4}}, {{4, 5}, {2, 4}}, {{5, 3}, {3, 4}}};

template <typename Nil>
const uint32_t mappings<Wedge, Nil>::vertex_mapping[6][3][3] = {
    {{0, 2, 3}, {0, 1, 0}, {0, 1, 3}}, {{0, 1, 4}, {1, 0, 0}, {0, 1, 2}},
    {{1, 2, 5}, {1, 0, 0}, {0, 2, 3}},

    {{3, 6, 8}, {1, 0, 1}, {1, 3, 4}}, {{4, 6, 7}, {1, 1, 0}, {1, 2, 4}},
    {{5, 7, 8}, {1, 1, 0}, {2, 3, 4}}};

template <typename Nil>
const uint32_t mappings<Wedge, Nil>::face_mapping[3][2] = {
    {0, 1}, {0, 2}, {1, 2}};

typedef mappings<Wedge, void> wedge_mappings;

} // namespace detail

class Wedge : public detail::wedge_mappings {
 public:
  template <typename... Types>
  Wedge(const vector_t& v0, Types... verts)
      : vertices{{v0, verts...}}, faces(), center(), volume() {
    init();
  }

  Wedge(const std::array<vector_t, 6>& verts)
      : vertices(verts), faces(), center(), volume() {
    init();
  }

  Wedge()
      : vertices{{vector_t::Zero(), vector_t::Zero(), vector_t::Zero(),
                  vector_t::Zero(), vector_t::Zero(), vector_t::Zero()}}
      , faces()
      , center()
      , volume() {}

  void apply(const Transformation& t) {
    for (auto& v : vertices)
      v = t.scaling * (v + t.translation);

    for (auto& f : faces)
      f.apply(t);

    center = scalar_t(1.0 / 6.0) *
             std::accumulate(
                 vertices.begin(), vertices.end(), vector_t::Zero().eval());

    volume = calcVolume();
  }

  scalar_t surfaceArea() const {
    scalar_t area(0);
    for (const auto& f : faces)
      area += f.area;

    return area;
  }

 private:
  void init() {
    // All faces of the wedge are stored as quadrilaterals, so an
    // additional point is inserted between v0 and v1.
    // 0: v2, v1, v0, v02
    faces[0] = Quadrilateral(
        vertices[2], vertices[1], vertices[0],
        scalar_t(0.5) * (vertices[0] + vertices[2]));

    // 1: v0, v1, v4, v3
    faces[1] =
        Quadrilateral(vertices[0], vertices[1], vertices[4], vertices[3]);

    // 2: v1, v2, v5, v4
    faces[2] =
        Quadrilateral(vertices[1], vertices[2], vertices[5], vertices[4]);

    // 3: v2, v0, v3, v5
    faces[3] =
        Quadrilateral(vertices[2], vertices[0], vertices[3], vertices[5]);

    // All faces of the wedge are stored as quadrilaterals, so an
    // additional point is inserted between v3 and v5.
    // 4: v3, v4, v5, v53
    faces[4] = Quadrilateral(
        vertices[3], vertices[4], vertices[5],
        scalar_t(0.5) * (vertices[5] + vertices[3]));

    center = scalar_t(1.0 / 6.0) *
             std::accumulate(
                 vertices.begin(), vertices.end(), vector_t::Zero().eval());

    volume = calcVolume();
  }

  scalar_t calcVolume() const {
    // The wedge is treated as a degenerate hexahedron here by adding
    // two fake vertices v02 and v35.
    vector_t diagonal(vertices[5] - vertices[0]);

    return scalar_t(1.0 / 6.0) *
           (diagonal.dot(
               ((vertices[1] - vertices[0]).cross(vertices[2] - vertices[4])) +
               ((vertices[3] - vertices[0])
                    .cross(
                        vertices[4] -
                        scalar_t(0.5) * (vertices[3] + vertices[5]))) +
               ((scalar_t(0.5) * (vertices[0] + vertices[2]) - vertices[0])
                    .cross(
                        scalar_t(0.5) * (vertices[3] + vertices[5]) -
                        vertices[2]))));
  }

 public:
  std::array<vector_t, 6> vertices;
  std::array<Quadrilateral, 5> faces;
  vector_t center;
  scalar_t volume;
};

namespace detail {

template <typename Nil>
struct mappings<Hexahedron, Nil> {
  // Map edges of a hexahedron to vertices and faces.
  static const uint32_t edge_mapping[12][2][2];

  // Map vertices of a hexahedron to edges and faces.
  // 0: local IDs of the edges intersecting at this vertex
  // 1: 0 if the edge is pointing away from the vertex, 1 otherwise
  // 2: faces joining at the vertex
  static const uint32_t vertex_mapping[8][3][3];

  // This mapping contains the three sets of the two edges for each of the
  // faces joining at a vertex. The indices are mapped to the local edge IDs
  // using the first value field of the 'vertex_mapping' table.
  static const uint32_t face_mapping[3][2];
};

template <typename Nil>
const uint32_t mappings<Hexahedron, Nil>::edge_mapping[12][2][2] = {
    {{0, 1}, {0, 1}}, {{1, 2}, {0, 2}}, {{2, 3}, {0, 3}}, {{3, 0}, {0, 4}},

    {{0, 4}, {1, 4}}, {{1, 5}, {1, 2}}, {{2, 6}, {2, 3}}, {{3, 7}, {3, 4}},

    {{4, 5}, {1, 5}}, {{5, 6}, {2, 5}}, {{6, 7}, {3, 5}}, {{7, 4}, {4, 5}}};

template <typename Nil>
const uint32_t mappings<Hexahedron, Nil>::vertex_mapping[8][3][3] = {
    {{0, 3, 4}, {0, 1, 0}, {0, 1, 4}},  {{0, 1, 5}, {1, 0, 0}, {0, 1, 2}},
    {{1, 2, 6}, {1, 0, 0}, {0, 2, 3}},  {{2, 3, 7}, {1, 0, 0}, {0, 3, 4}},

    {{4, 8, 11}, {1, 0, 1}, {1, 4, 5}}, {{5, 8, 9}, {1, 1, 0}, {1, 2, 5}},
    {{6, 9, 10}, {1, 1, 0}, {2, 3, 5}}, {{7, 10, 11}, {1, 1, 0}, {3, 4, 5}}};

template <typename Nil>
const uint32_t mappings<Hexahedron, Nil>::face_mapping[3][2] = {
    {0, 1}, {0, 2}, {1, 2}};

typedef mappings<Hexahedron, void> hex_mappings;

} // namespace detail

class Hexahedron : public detail::hex_mappings {
 public:
  template <typename... Types>
  Hexahedron(const vector_t& v0, Types... verts)
      : vertices{{v0, verts...}}, faces(), center(), volume() {
    init();
  }

  Hexahedron(const std::array<vector_t, 8>& verts)
      : vertices(verts), faces(), center(), volume() {
    init();
  }

  void apply(const Transformation& t) {
    for (auto& v : vertices)
      v = t.scaling * (v + t.translation);

    for (auto& f : faces)
      f.apply(t);

    center = scalar_t(1.0 / 8.0) *
             std::accumulate(
                 vertices.begin(), vertices.end(), vector_t::Zero().eval());

    volume = calcVolume();
  }

  scalar_t surfaceArea() const {
    scalar_t area(0);
    for (const auto& f : faces)
      area += f.area;

    return area;
  }

 private:
  void init() {
    // 0: v3, v2, v1, v0
    faces[0] =
        Quadrilateral(vertices[3], vertices[2], vertices[1], vertices[0]);

    // 1: v0, v1, v5, v4
    faces[1] =
        Quadrilateral(vertices[0], vertices[1], vertices[5], vertices[4]);

    // 2: v1, v2, v6, v5
    faces[2] =
        Quadrilateral(vertices[1], vertices[2], vertices[6], vertices[5]);

    // 3: v2, v3, v7, v6
    faces[3] =
        Quadrilateral(vertices[2], vertices[3], vertices[7], vertices[6]);

    // 4: v3, v0, v4, v7
    faces[4] =
        Quadrilateral(vertices[3], vertices[0], vertices[4], vertices[7]);

    // 5: v4, v5, v6, v7
    faces[5] =
        Quadrilateral(vertices[4], vertices[5], vertices[6], vertices[7]);

    center = scalar_t(1.0 / 8.0) *
             std::accumulate(
                 vertices.begin(), vertices.end(), vector_t::Zero().eval());

    volume = calcVolume();
  }

  scalar_t calcVolume() const {
    vector_t diagonal(vertices[6] - vertices[0]);

    return scalar_t(1.0 / 6.0) *
           diagonal.dot(
               ((vertices[1] - vertices[0]).cross(vertices[2] - vertices[5])) +
               ((vertices[4] - vertices[0]).cross(vertices[5] - vertices[7])) +
               ((vertices[3] - vertices[0]).cross(vertices[7] - vertices[2])));
  }

 public:
  std::array<vector_t, 8> vertices;
  std::array<Quadrilateral, 6> faces;
  vector_t center;
  scalar_t volume;
};

class Sphere {
 public:
  Sphere(const vector_t& c, scalar_t r)
      : center(c), radius(r), volume(scalar_t(4.0 / 3.0 * pi) * r * r * r) {}

  scalar_t capVolume(scalar_t h) const {
    if (h <= scalar_t(0))
      return scalar_t(0);
    else if (h >= scalar_t(2) * radius)
      return volume;
    else
      return scalar_t(pi / 3.0) * h * h * (scalar_t(3) * radius - h);
  }

  scalar_t capSurfaceArea(scalar_t h) const {
    if (h <= scalar_t(0))
      return scalar_t(0);
    else if (h >= scalar_t(2) * radius)
      return surfaceArea();
    else
      return scalar_t(2 * pi) * radius * h;
  }

  scalar_t diskArea(scalar_t h) const {
    if (h <= scalar_t(0) || h >= scalar_t(2) * radius)
      return scalar_t(0);
    else
      return pi * h * (scalar_t(2) * radius - h);
  }

  scalar_t surfaceArea() const {
    return (scalar_t(4) * pi) * (radius * radius);
  }

 public:
  vector_t center;
  scalar_t radius;
  scalar_t volume;
};

class Plane {
 public:
  Plane(const vector_t& c, const vector_t& n) : center(c), normal(n) {}

 public:
  vector_t center;
  vector_t normal;
};

class AABB {
 public:
  AABB(
      const vector_t& minimum =
          vector_t::Constant(std::numeric_limits<scalar_t>::infinity()),
      const vector_t& maximum =
          vector_t::Constant(-std::numeric_limits<scalar_t>::infinity()))
      : min(minimum), max(maximum) {}

  bool intersects(const AABB& aabb) const {
    if ((min.array() > aabb.max.array()).any() ||
        (max.array() < aabb.min.array()).any())
      return false;

    return true;
  }

  AABB overlap(const AABB& aabb) const {
    return AABB(min.cwiseMax(aabb.min), max.cwiseMin(aabb.max));
  }

  bool contains(const vector_t& p) const {
    if ((p.array() < min.array()).any() || (p.array() > max.array()).any())
      return false;

    return true;
  }

  void include(const vector_t& point) {
    min = min.cwiseMin(point);
    max = max.cwiseMax(point);
  }

  template <size_t N>
  void include(const std::array<vector_t, N>& points) {
    for (const auto& p : points)
      include(p);
  }

  scalar_t volume() const {
    vector_t size(max - min);

    return size[0] * size[1] * size[2];
  }

 public:
  vector_t min, max;
};

// Decomposition of a tetrahedron into 4 tetrahedra.
inline void decompose(
    const Tetrahedron& tet, std::array<Tetrahedron, 4>& tets) {
  tets[0] = Tetrahedron(
      tet.vertices[0], tet.vertices[1], tet.vertices[2], tet.center);

  tets[1] = Tetrahedron(
      tet.vertices[0], tet.vertices[1], tet.center, tet.vertices[3]);

  tets[2] = Tetrahedron(
      tet.vertices[1], tet.vertices[2], tet.center, tet.vertices[3]);

  tets[3] = Tetrahedron(
      tet.vertices[2], tet.vertices[0], tet.center, tet.vertices[3]);
}

// Decomposition of a hexahedron into 2 wedges.
inline void decompose(const Hexahedron& hex, std::array<Wedge, 2>& wedges) {
  wedges[0] = Wedge(
      hex.vertices[0], hex.vertices[1], hex.vertices[2], hex.vertices[4],
      hex.vertices[5], hex.vertices[6]);

  wedges[1] = Wedge(
      hex.vertices[0], hex.vertices[2], hex.vertices[3], hex.vertices[4],
      hex.vertices[6], hex.vertices[7]);
}

// Decomposition of a hexahedron into 5 tetrahedra.
inline void decompose(const Hexahedron& hex, std::array<Tetrahedron, 5>& tets) {
  tets[0] = Tetrahedron(
      hex.vertices[0], hex.vertices[1], hex.vertices[2], hex.vertices[5]);

  tets[1] = Tetrahedron(
      hex.vertices[0], hex.vertices[2], hex.vertices[7], hex.vertices[5]);

  tets[2] = Tetrahedron(
      hex.vertices[0], hex.vertices[2], hex.vertices[3], hex.vertices[7]);

  tets[3] = Tetrahedron(
      hex.vertices[0], hex.vertices[5], hex.vertices[7], hex.vertices[4]);

  tets[4] = Tetrahedron(
      hex.vertices[2], hex.vertices[7], hex.vertices[5], hex.vertices[6]);
}

// Decomposition of a hexahedron into 6 tetrahedra.
inline void decompose(const Hexahedron& hex, std::array<Tetrahedron, 6>& tets) {
  tets[0] = Tetrahedron(
      hex.vertices[0], hex.vertices[5], hex.vertices[7], hex.vertices[4]);

  tets[1] = Tetrahedron(
      hex.vertices[0], hex.vertices[1], hex.vertices[7], hex.vertices[5]);

  tets[2] = Tetrahedron(
      hex.vertices[1], hex.vertices[6], hex.vertices[7], hex.vertices[5]);

  tets[3] = Tetrahedron(
      hex.vertices[0], hex.vertices[7], hex.vertices[2], hex.vertices[3]);

  tets[4] = Tetrahedron(
      hex.vertices[0], hex.vertices[7], hex.vertices[1], hex.vertices[2]);

  tets[5] = Tetrahedron(
      hex.vertices[1], hex.vertices[7], hex.vertices[6], hex.vertices[2]);
}

inline bool contains(const Sphere& s, const vector_t& p) {
  return (s.center - p).squaredNorm() <= s.radius * s.radius;
}

// The (convex!) polygon is assumed to be planar, making this a 2D problem.
// Check the projection of the point onto the plane of the polygon for
// containment within the polygon.
template <size_t VertexCount>
bool contains(const Polygon<VertexCount>& poly, const vector_t& point) {
  const vector_t proj(
      point - poly.normal.dot(point - poly.center) * poly.normal);

  for (size_t n = 0; n < poly.vertices.size(); ++n) {
    const auto& v0 = poly.vertices[n];
    const auto& v1 = poly.vertices[(n + 1) % poly.vertices.size()];
    vector_t base(scalar_t(0.5) * (v0 + v1));
    vector_t edge(v1 - v0);

    // Note: Only the sign of the projection is of interest, so this vector
    // does not have to be normalized.
    vector_t dir(edge.cross(poly.normal));

    // Check whether the projection of the point lies inside of the
    // polygon.
    if (dir.dot(proj - base) > scalar_t(0)) return false;
  }

  return true;
}

inline bool contains(const Tetrahedron& tet, const vector_t& p) {
  for (const auto& f : tet.faces)
    if (f.normal.dot(p - f.center) > scalar_t(0)) return false;

  return true;
}

inline bool contains(const Wedge& wedge, const vector_t& p) {
  for (const auto& f : wedge.faces)
    if (f.normal.dot(p - f.center) > scalar_t(0)) return false;

  return true;
}

inline bool contains(const Hexahedron& hex, const vector_t& p) {
  for (const auto& f : hex.faces)
    if (f.normal.dot(p - f.center) > scalar_t(0)) return false;

  return true;
}

inline bool intersect(const Sphere& s, const Plane& p) {
  scalar_t proj = p.normal.dot(s.center - p.center);

  return proj * proj - s.radius * s.radius < scalar_t(0);
}

template <size_t VertexCount>
inline bool intersect(const Sphere& s, const Polygon<VertexCount>& poly) {
  return intersect(s, Plane(poly.center, poly.normal)) &&
         contains(poly, s.center);
}

inline std::pair<std::array<scalar_t, 2>, size_t> lineSphereIntersection(
    const vector_t& origin, const vector_t& direction, const Sphere& s) {
  std::array<scalar_t, 2> solutions = {
      {std::numeric_limits<scalar_t>::infinity(),
       std::numeric_limits<scalar_t>::infinity()}};

  vector_t originRel(origin - s.center);
  scalar_t a = direction.squaredNorm();

  if (a == scalar_t(0)) return std::make_pair(solutions, 0);

  scalar_t b = scalar_t(2) * direction.dot(originRel);
  scalar_t c = originRel.squaredNorm() - s.radius * s.radius;

  scalar_t discriminant = b * b - scalar_t(4) * a * c;
  if (discriminant > scalar_t(0)) {
    // Two real roots.
    scalar_t q =
        scalar_t(-0.5) * (b + std::copysign(std::sqrt(discriminant), b));

    solutions[0] = q / a;
    solutions[1] = c / q;

    if (solutions[0] > solutions[1]) std::swap(solutions[0], solutions[1]);

    return std::make_pair(solutions, 2);
  } else if (std::abs(discriminant) == scalar_t(0)) {
    // Double real root.
    solutions[0] = (scalar_t(-0.5) * b) / a;
    solutions[1] = solutions[0];

    return std::make_pair(solutions, 1);
  } else {
    // No real roots.
    return std::make_pair(solutions, 0);
  }
}

namespace detail {

// Calculate the volume of a regularized spherical wedge defined by the radius,
// the distance of the intersection point from the center of the sphere and the
// angle.
inline scalar_t regularizedWedge(scalar_t r, scalar_t d, scalar_t alpha) {
#ifndef NDEBUG
  // Clamp slight deviations of the angle to valid range.
  if (alpha < scalar_t(0) && alpha > -detail::tinyEpsilon) alpha = scalar_t(0);

  if (alpha > scalar_t(0.5 * pi) && alpha < scalar_t(0.5 * pi) + tinyEpsilon)
    alpha = scalar_t(0.5 * pi);
#endif

  // Check the parameters for validity (debug version only).
  assert(r > scalar_t(0));
  assert(d >= scalar_t(0) && d <= r);
  assert(alpha >= scalar_t(0) && alpha <= scalar_t(0.5 * pi));

  const scalar_t sinAlpha = std::sin(alpha);
  const scalar_t cosAlpha = std::cos(alpha);

  const scalar_t a = d * sinAlpha;
  const scalar_t b = std::sqrt(std::abs(r * r - d * d));
  const scalar_t c = d * cosAlpha;

  return scalar_t(1.0 / 3.0) * a * b * c +
         a * (scalar_t(1.0 / 3.0) * a * a - r * r) * std::atan2(b, c) +
         scalar_t(2.0 / 3.0) * r * r * r *
             std::atan2(sinAlpha * b, cosAlpha * r);
}

// Wrapper around the above function handling correctly handling the case of
// alpha > pi/2 and negative z.
inline scalar_t regularizedWedge(
    scalar_t r, scalar_t d, scalar_t alpha, scalar_t z) {
  if (z >= scalar_t(0)) {
    if (alpha > scalar_t(0.5 * pi)) {
      scalar_t h = r - z;

      return scalar_t(pi / 3.0) * h * h * (scalar_t(3) * r - h) -
             regularizedWedge(r, d, pi - alpha);
    } else {
      return regularizedWedge(r, d, alpha);
    }
  } else {
    scalar_t vHem = scalar_t(2.0 / 3.0 * pi) * r * r * r;

    if (alpha > scalar_t(0.5 * pi)) {
      return vHem - regularizedWedge(r, d, pi - alpha);
    } else {
      scalar_t h = r + z;
      scalar_t vCap = scalar_t(pi / 3.0) * h * h * (scalar_t(3) * r - h);

      return vHem - (vCap - regularizedWedge(r, d, alpha));
    }
  }
}

// Calculate the surface area of a regularized spherical wedge defined by the
// radius, the distance of the intersection point from the center of the sphere
// and the angle.
// Ref: Gibson, K. D. & Scheraga, H. A.: Exact calculation of the volume and
//      surface area of fused hard-sphere molecules with unequal atomic radii,
//      Molecular Physics, 1987, 62, 1247-1265
inline scalar_t regularizedWedgeArea(scalar_t r, scalar_t z, scalar_t alpha) {
#ifndef NDEBUG
  // Clamp slight deviations of the angle to valid range.
  if (alpha < scalar_t(0) && alpha > -detail::tinyEpsilon) alpha = scalar_t(0);

  if (alpha > pi && alpha < pi + tinyEpsilon) alpha = pi;
#endif

  // Check the parameters for validity (debug version only).
  assert(r > scalar_t(0));
  assert(z >= -r && z <= r);
  assert(alpha >= scalar_t(0) && alpha <= pi);

  if (alpha < tinyEpsilon || std::abs(r * r - z * z) <= tinyEpsilon)
    return scalar_t(0);

  const scalar_t sinAlpha = std::sin(alpha);
  const scalar_t cosAlpha = std::cos(alpha);
  const scalar_t factor = scalar_t(1) / std::sqrt(std::abs(r * r - z * z));

  // Clamp slight deviations of the argument to acos() to valid range.
  const scalar_t arg0 = clamp(
      r * cosAlpha * factor, scalar_t(-1), scalar_t(1), detail::tinyEpsilon);

  const scalar_t arg1 = clamp(
      (z * cosAlpha * factor) / sinAlpha, scalar_t(-1), scalar_t(1),
      detail::tinyEpsilon);

  // Check the argument to acos() for validity (debug version only).
  assert(scalar_t(-1) <= arg0 && arg0 <= scalar_t(1));
  assert(scalar_t(-1) <= arg1 && arg1 <= scalar_t(1));

  return scalar_t(2) * r * r * std::acos(arg0) -
         scalar_t(2) * r * z * std::acos(arg1);
}

} // namespace detail

// Depending on the dimensionality, either the volume or external surface area
// of the general wedge is computed.
template <size_t Dim>
inline scalar_t generalWedge(
    const Sphere& s, const Plane& p0, const Plane& p1, const vector_t& d) {
  static_assert(
      Dim == 2 || Dim == 3, "Invalid dimensionality, must be 2 or 3.");

  scalar_t dist(d.stableNorm());

  if (dist < detail::tinyEpsilon) {
    // The wedge (almost) touches the center, the volume depends only on
    // the angle.
    scalar_t angle = pi - detail::angle(p0.normal, p1.normal);

    if (Dim == 2) {
      return scalar_t(2) * s.radius * s.radius * angle;
    } else {
      return scalar_t(2.0 / 3.0) * s.radius * s.radius * s.radius * angle;
    }
  }

  scalar_t s0 = d.dot(p0.normal);
  scalar_t s1 = d.dot(p1.normal);

  // Detect degenerated general spherical wedge that can be treated as
  // a regularized spherical wedge.
  if (std::abs(s0) < detail::tinyEpsilon ||
      std::abs(s1) < detail::tinyEpsilon) {
    scalar_t angle = pi - detail::angle(p0.normal, p1.normal);

    if (Dim == 2) {
      return detail::regularizedWedgeArea(
          s.radius, std::abs(s0) > std::abs(s1) ? s0 : s1, angle);
    } else {
      return detail::regularizedWedge(
          s.radius, dist, angle, std::abs(s0) > std::abs(s1) ? s0 : s1);
    }
  }

  vector_t dUnit(d * (scalar_t(1) / dist));
  if (dist < detail::largeEpsilon)
    dUnit = detail::gramSchmidt(p0.normal.cross(p1.normal), dUnit)[1];

  // Check the planes specify a valid setup (debug version only).
  assert(p0.normal.dot(p1.center - p0.center) <= scalar_t(0));
  assert(p1.normal.dot(p0.center - p1.center) <= scalar_t(0));

  // Calculate the angles between the vector from the sphere center
  // to the intersection line and the normal vectors of the two planes.
  scalar_t alpha0 = detail::angle(p0.normal, dUnit);
  scalar_t alpha1 = detail::angle(p1.normal, dUnit);

  scalar_t dir0 = dUnit.dot((s.center + d) - p0.center);
  scalar_t dir1 = dUnit.dot((s.center + d) - p1.center);

  if (s0 >= scalar_t(0) && s1 >= scalar_t(0)) {
    alpha0 = scalar_t(0.5 * pi) - std::copysign(alpha0, dir0);
    alpha1 = scalar_t(0.5 * pi) - std::copysign(alpha1, dir1);

    if (Dim == 2) {
      return detail::regularizedWedgeArea(s.radius, s0, alpha0) +
             detail::regularizedWedgeArea(s.radius, s1, alpha1);
    } else {
      return detail::regularizedWedge(s.radius, dist, alpha0, s0) +
             detail::regularizedWedge(s.radius, dist, alpha1, s1);
    }
  } else if (s0 < scalar_t(0) && s1 < scalar_t(0)) {
    alpha0 =
        scalar_t(0.5 * pi) + std::copysign(scalar_t(1), dir0) * (alpha0 - pi);

    alpha1 =
        scalar_t(0.5 * pi) + std::copysign(scalar_t(1), dir1) * (alpha1 - pi);

    if (Dim == 2) {
      return s.surfaceArea() -
             (detail::regularizedWedgeArea(s.radius, -s0, alpha0) +
              detail::regularizedWedgeArea(s.radius, -s1, alpha1));
    } else {
      return s.volume - (detail::regularizedWedge(s.radius, dist, alpha0, -s0) +
                         detail::regularizedWedge(s.radius, dist, alpha1, -s1));
    }
  } else {
    alpha0 = scalar_t(0.5 * pi) -
             std::copysign(scalar_t(1), dir0 * s0) *
                 (alpha0 - (s0 < scalar_t(0) ? pi : scalar_t(0)));

    alpha1 = scalar_t(0.5 * pi) -
             std::copysign(scalar_t(1), dir1 * s1) *
                 (alpha1 - (s1 < scalar_t(0) ? pi : scalar_t(0)));

    if (Dim == 2) {
      scalar_t area0 =
          detail::regularizedWedgeArea(s.radius, std::abs(s0), alpha0);

      scalar_t area1 =
          detail::regularizedWedgeArea(s.radius, std::abs(s1), alpha1);

      return std::max(area0, area1) - std::min(area0, area1);
    } else {
      scalar_t volume0 =
          detail::regularizedWedge(s.radius, dist, alpha0, std::abs(s0));

      scalar_t volume1 =
          detail::regularizedWedge(s.radius, dist, alpha1, std::abs(s1));

      return std::max(volume0, volume1) - std::min(volume0, volume1);
    }
  }
}

template <typename T>
struct array_size;

template <typename T, size_t N>
struct array_size<std::array<T, N>> {
  static constexpr size_t value() {
    return N;
  }
};

template <typename Element>
constexpr size_t nrEdges() {
  return std::is_same<Element, Hexahedron>::value
             ? 12
             : (std::is_same<Element, Wedge>::value
                    ? 9
                    : (std::is_same<Element, Tetrahedron>::value ? 6 : -1));
}

// Workaround for the Intel compiler, as it does not yet support constexpr for
// template arguments.
template <typename Element>
struct element_trait {
  static const size_t nrVertices =
      array_size<decltype(Element::vertices)>::value();

  static const size_t nrFaces = array_size<decltype(Element::faces)>::value();
};

// Depending on the dimensionality, either the volume or external surface area
// of the general wedge is computed.
template <size_t Dim, typename Element>
scalar_t generalWedge(
    const Sphere& sphere, const Element& element, size_t edge,
    const std::array<std::array<vector_t, 2>, nrEdges<Element>()>&
        intersections) {
  static_assert(
      Dim == 2 || Dim == 3, "Invalid dimensionality, must be 2 or 3.");

  const auto& f0 = element.faces[Element::edge_mapping[edge][1][0]];
  const auto& f1 = element.faces[Element::edge_mapping[edge][1][1]];

  vector_t edgeCenter(
      scalar_t(0.5) * ((intersections[edge][0] +
                        element.vertices[Element::edge_mapping[edge][0][0]]) +
                       (intersections[edge][1] +
                        element.vertices[Element::edge_mapping[edge][0][1]])));

  Plane p0(f0.center, f0.normal);
  Plane p1(f1.center, f1.normal);

  return generalWedge<Dim>(sphere, p0, p1, edgeCenter - sphere.center);
}

template <typename Element>
scalar_t overlap(const Sphere& sOrig, const Element& elementOrig) {
  static_assert(
      std::is_same<Element, Tetrahedron>::value ||
          std::is_same<Element, Wedge>::value ||
          std::is_same<Element, Hexahedron>::value,
      "Invalid element type detected.");

  // Construct AABBs and perform a coarse overlap detection.
  AABB sAABB(
      sOrig.center - vector_t::Constant(sOrig.radius),
      sOrig.center + vector_t::Constant(sOrig.radius));

  AABB eAABB;
  eAABB.include(elementOrig.vertices);

  if (!sAABB.intersects(eAABB)) return scalar_t(0);

  // Use scaled and shifted versions of the sphere and the element.
  Transformation transformation(-sOrig.center, scalar_t(1) / sOrig.radius);

  Sphere s(vector_t::Zero(), scalar_t(1));

  Element element(elementOrig);
  element.apply(transformation);

  // Constants: Number of vertices and faces.
  static const size_t nrVertices = element_trait<Element>::nrVertices;
  static const size_t nrFaces = element_trait<Element>::nrFaces;

  size_t vOverlap = 0;
  // Check whether the vertices lie on or outside of the sphere.
  for (const auto& vertex : element.vertices)
    if ((s.center - vertex).squaredNorm() <= s.radius * s.radius) ++vOverlap;

  // Check for trivial case: All vertices inside of the sphere, resulting in
  // a full overlap.
  if (vOverlap == nrVertices) return elementOrig.volume;

  // Sanity check: All faces of the mesh element have to be planar.
  for (const auto& face : element.faces)
    if (!face.isPlanar())
      throw std::runtime_error("Non-planer face detected in element!");

  // Sets of overlapping primitives.
  std::bitset<nrVertices> vMarked;
  std::bitset<nrEdges<Element>()> eMarked;
  std::bitset<nrFaces> fMarked;

  // Initial value: Volume of the full sphere.
  scalar_t result = s.volume;

  // The intersection points between the single edges and the sphere, this
  // is needed later on.
  std::array<std::array<vector_t, 2>, nrEdges<Element>()> eIntersections;

  // Process all edges of the element.
  for (size_t n = 0; n < nrEdges<Element>(); ++n) {
    vector_t start(element.vertices[Element::edge_mapping[n][0][0]]);
    vector_t direction(
        element.vertices[Element::edge_mapping[n][0][1]] - start);

    auto solutions = lineSphereIntersection(start, direction, s);

    // No intersection between the edge and the sphere, where intersection
    // points close to the surface of the sphere are ignored.
    // Or:
    // The sphere cuts the edge twice, no vertex is inside of the
    // sphere, but the case of the edge only touching the sphere has to
    // be avoided.
    if (!solutions.second ||
        (solutions.first[0] >= scalar_t(1) - detail::mediumEpsilon) ||
        solutions.first[1] <= detail::mediumEpsilon ||
        (solutions.first[0] > scalar_t(0) && solutions.first[1] < scalar_t(1) &&
         (solutions.first[1] - solutions.first[0] < detail::largeEpsilon))) {
      continue;
    } else {
      vMarked[Element::edge_mapping[n][0][0]] =
          solutions.first[0] < scalar_t(0);

      vMarked[Element::edge_mapping[n][0][1]] =
          solutions.first[1] > scalar_t(1);
    }

    // Store the two intersection points of the edge with the sphere for
    // later usage.
    eIntersections[n][0] =
        solutions.first[0] * direction +
        (start - element.vertices[Element::edge_mapping[n][0][0]]);

    eIntersections[n][1] =
        solutions.first[1] * direction +
        (start - element.vertices[Element::edge_mapping[n][0][1]]);

    eMarked[n] = true;

    // If the edge is marked as having an overlap, the two faces forming it
    // have to be marked as well.
    fMarked[Element::edge_mapping[n][1][0]] = true;
    fMarked[Element::edge_mapping[n][1][1]] = true;
  }

  // Check whether the dependencies for a vertex intersection are fulfilled.
  for (size_t n = 0; n < nrVertices; ++n) {
    if (!vMarked[n]) continue;

    bool edgesValid = true;
    for (size_t eN = 0; eN < 3; ++eN) {
      size_t edgeId = Element::vertex_mapping[n][0][eN];
      edgesValid &= eMarked[edgeId];
    }

    // If not all three edges intersecting at this vertex where marked, the
    // sphere is only touching.
    if (!edgesValid) vMarked[n] = false;
  }

  // Process all faces of the element, ignoring the edges as those where
  // already checked above.
  for (size_t n = 0; n < nrFaces; ++n)
    if (intersect(s, element.faces[n])) fMarked[n] = true;

  // Trivial case: The center of the sphere overlaps the element, but the
  // sphere does not intersect any of the faces of the element, meaning the
  // sphere is completely contained within the element.
  if (!fMarked.count() && contains(element, s.center)) return sOrig.volume;

  // Spurious intersection: The initial intersection test was positive, but
  // the detailed checks revealed no overlap.
  if (!vMarked.count() && !eMarked.count() && !fMarked.count())
    return scalar_t(0);

  // Iterate over all the marked faces and subtract the volume of the cap cut
  // off by the plane.
  for (size_t n = 0; n < nrFaces; ++n) {
    if (!fMarked[n]) continue;

    const auto& f = element.faces[n];
    scalar_t dist = f.normal.dot(s.center - f.center);
    scalar_t vCap = s.capVolume(s.radius + dist);

    result -= vCap;
  }

  // Handle the edges and add back the volume subtracted twice above in the
  // processing of the faces.
  for (size_t n = 0; n < nrEdges<Element>(); ++n) {
    if (!eMarked[n]) continue;

    scalar_t edgeCorrection =
        generalWedge<3, Element>(s, element, n, eIntersections);

    result += edgeCorrection;
  }

  // Handle the vertices and subtract the volume added twice above in the
  // processing of the edges.
  for (size_t n = 0; n < nrVertices; ++n) {
    if (!vMarked[n]) continue;

    // Collect the points where the three edges intersecting at this
    // vertex intersect the sphere.
    // Both the relative and the absolute positions are required.
    std::array<vector_t, 3> intersectionPointsRelative;
    std::array<vector_t, 3> intersectionPoints;
    for (size_t e = 0; e < 3; ++e) {
      auto edgeIdx = Element::vertex_mapping[n][0][e];
      intersectionPointsRelative[e] =
          eIntersections[edgeIdx][Element::vertex_mapping[n][1][e]];

      intersectionPoints[e] =
          intersectionPointsRelative[e] + element.vertices[n];
    }

    // This triangle is constructed by hand to have more freedom of how
    // the normal vector is calculated.
    Triangle coneTria;
    coneTria.vertices = {
        {intersectionPoints[0], intersectionPoints[1], intersectionPoints[2]}};

    coneTria.center = scalar_t(1.0 / 3.0) * std::accumulate(
                                                intersectionPoints.begin(),
                                                intersectionPoints.end(),
                                                vector_t::Zero().eval());

    // Calculate the normal of the triangle defined by the intersection
    // points in relative coordinates to improve accuracy.
    // Also use double the normal precision to calculate this normal.
    coneTria.normal = detail::triangleNormal(
        intersectionPointsRelative[0], intersectionPointsRelative[1],
        intersectionPointsRelative[2]);

    // The area of this triangle is never needed, so it is set to an
    // invalid value.
    coneTria.area = std::numeric_limits<scalar_t>::infinity();

    std::array<std::pair<size_t, scalar_t>, 3> distances;
    for (size_t i = 0; i < 3; ++i)
      distances[i] =
          std::make_pair(i, intersectionPointsRelative[i].squaredNorm());

    std::sort(
        distances.begin(), distances.end(),
        [](const std::pair<size_t, scalar_t>& a,
           const std::pair<size_t, scalar_t>& b) -> bool {
          return a.second < b.second;
        });

    if (distances[1].second < distances[2].second * detail::largeEpsilon) {
      // Use the general spherical wedge defined by the edge with the
      // non-degenerated intersection point and the normals of the
      // two faces forming it.
      scalar_t correction = generalWedge<3, Element>(
          s, element, Element::vertex_mapping[n][0][distances[2].first],
          eIntersections);

      result -= correction;

      continue;
    }

    scalar_t tipTetVolume =
        scalar_t(1.0 / 6.0) *
        std::abs(-intersectionPointsRelative[2].dot(
            (intersectionPointsRelative[0] - intersectionPointsRelative[2])
                .cross(
                    intersectionPointsRelative[1] -
                    intersectionPointsRelative[2])));

    // Make sure the normal points in the right direction i.e. away from
    // the center of the element.
    if (coneTria.normal.dot(element.center - coneTria.center) > scalar_t(0)) {
      coneTria.normal = -coneTria.normal;
    }

    Plane plane(coneTria.center, coneTria.normal);

    scalar_t dist = coneTria.normal.dot(s.center - coneTria.center);
    scalar_t capVolume = s.capVolume(s.radius + dist);

    // The cap volume is tiny, so the corrections will be even smaller.
    // There is no way to actually calculate them with reasonable
    // precision, so just the volume of the tetrahedron at the tip is
    // used.
    if (capVolume < detail::tinyEpsilon) {
      result -= tipTetVolume;
      continue;
    }

    // Calculate the volume of the three spherical segments between
    // the faces joining at the vertex and the plane through the
    // intersection points.
    scalar_t segmentVolume = 0;

    for (size_t e = 0; e < 3; ++e) {
      const auto& f = element.faces[Element::vertex_mapping[n][2][e]];
      uint32_t e0 = Element::face_mapping[e][0];
      uint32_t e1 = Element::face_mapping[e][1];

      vector_t center(
          scalar_t(0.5) * (intersectionPoints[e0] + intersectionPoints[e1]));

      scalar_t wedgeVolume = generalWedge<3>(
          s, plane, Plane(f.center, -f.normal), center - s.center);

      segmentVolume += wedgeVolume;
    }

    // Calculate the volume of the cone and clamp it to zero.
    scalar_t coneVolume =
        std::max(tipTetVolume + capVolume - segmentVolume, scalar_t(0));

    // Sanity check: detect negative cone volume.
    assert(coneVolume > -std::sqrt(detail::tinyEpsilon));

    result -= coneVolume;

    // Sanity check: detect negative intermediate result.
    assert(result > -std::sqrt(detail::tinyEpsilon));
  }

  // In case of different sized objects the error can become quite large,
  // so a relative limit is used.
  scalar_t maxOverlap = std::min(s.volume, element.volume);
  const scalar_t limit(
      std::sqrt(std::numeric_limits<scalar_t>::epsilon()) * maxOverlap);

  // Clamp tiny negative volumes to zero.
  if (result < scalar_t(0) && result > -limit) return scalar_t(0);

  // Clamp results slightly too large.
  if (result > maxOverlap && result - maxOverlap < limit)
    return std::min(sOrig.volume, elementOrig.volume);

  // Perform a sanity check on the final result (debug version only).
  assert(result >= scalar_t(0) && result <= maxOverlap);

  // Scale the overlap volume back for the original objects.
  result = (result / s.volume) * sOrig.volume;

  return result;
}

template <typename Iterator>
scalar_t overlap(const Sphere& s, Iterator eBegin, Iterator eEnd) {
  scalar_t sum(0);

  for (Iterator it = eBegin; it != eEnd; ++it)
    sum += overlap(s, *it);

  return sum;
}

// Calculate the surface area of the sphere and the element that are contained
// within the common or intersecting part of the geometries, respectively.
// The returned array of size (N + 2), with N being the number of vertices,
// holds (in this order):
//   - surface area of the region of the sphere intersecting the element
//   - for each face of the element: area contained within the sphere
//   - total surface area of the element intersecting the sphere
template <
    typename Element, size_t NrFaces = element_trait<Element>::nrFaces + 2>
auto overlapArea(const Sphere& sOrig, const Element& elementOrig)
    -> std::array<scalar_t, NrFaces> {
  static_assert(
      NrFaces == element_trait<Element>::nrFaces + 2,
      "Invalid number of faces for the element provided.");

  static_assert(
      std::is_same<Element, Tetrahedron>::value ||
          std::is_same<Element, Wedge>::value ||
          std::is_same<Element, Hexahedron>::value,
      "Invalid element type detected.");

  // Constants: Number of vertices and faces.
  static const size_t nrVertices = element_trait<Element>::nrVertices;
  static const size_t nrFaces = element_trait<Element>::nrFaces;

  // Initial value: Zero overlap.
  std::array<scalar_t, nrFaces + 2> result;
  result.fill(scalar_t(0));

  // Construct AABBs and perform a coarse overlap detection.
  AABB sAABB(
      sOrig.center - vector_t::Constant(sOrig.radius),
      sOrig.center + vector_t::Constant(sOrig.radius));

  AABB eAABB;
  eAABB.include(elementOrig.vertices);

  if (!sAABB.intersects(eAABB)) return result;

  // Use scaled and shifted versions of the sphere and the element.
  Transformation transformation(-sOrig.center, scalar_t(1) / sOrig.radius);

  Sphere s(vector_t::Zero(), scalar_t(1));

  Element element(elementOrig);
  element.apply(transformation);

  size_t vOverlap = 0;
  // Check whether the vertices lie on or outside of the sphere.
  for (const auto& vertex : element.vertices)
    if ((s.center - vertex).squaredNorm() <= s.radius * s.radius) ++vOverlap;

  // Check for trivial case: All vertices inside of the sphere, resulting in
  // a full coverage of all faces.
  if (vOverlap == nrVertices) {
    for (size_t n = 0; n < nrFaces; ++n) {
      result[n + 1] = elementOrig.faces[n].area;
      result[nrFaces + 1] += elementOrig.faces[n].area;
    }

    return result;
  }

  // Sanity check: All faces of the mesh element have to be planar.
  for (const auto& face : element.faces)
    if (!face.isPlanar())
      throw std::runtime_error("Non-planer face detected in element!");

  // Sets of overlapping primitives.
  std::bitset<nrVertices> vMarked;
  std::bitset<nrEdges<Element>()> eMarked;
  std::bitset<nrFaces> fMarked;

  // The intersection points between the single edges and the sphere, this
  // is needed later on.
  std::array<std::array<vector_t, 2>, nrEdges<Element>()> eIntersections;

  // Cache the squared radius of the disk formed by the intersection between
  // the planes defined by each face and the sphere.
  std::array<scalar_t, nrFaces> intersectionRadiusSq;

  // Process all edges of the element.
  for (size_t n = 0; n < nrEdges<Element>(); ++n) {
    vector_t start(element.vertices[Element::edge_mapping[n][0][0]]);
    vector_t direction(
        element.vertices[Element::edge_mapping[n][0][1]] - start);

    auto solutions = lineSphereIntersection(start, direction, s);

    // No intersection between the edge and the sphere, where intersection
    // points close to the surface of the sphere are ignored.
    // Or:
    // The sphere cuts the edge twice, no vertex is inside of the
    // sphere, but the case of the edge only touching the sphere has to
    // be avoided.
    if (!solutions.second ||
        solutions.first[0] >= scalar_t(1) - detail::mediumEpsilon ||
        solutions.first[1] <= detail::mediumEpsilon ||
        (solutions.first[0] > scalar_t(0) && solutions.first[1] < scalar_t(1) &&
         solutions.first[1] - solutions.first[0] < detail::largeEpsilon)) {
      continue;
    } else {
      vMarked[Element::edge_mapping[n][0][0]] =
          solutions.first[0] < scalar_t(0);

      vMarked[Element::edge_mapping[n][0][1]] =
          solutions.first[1] > scalar_t(1);
    }

    // Store the two intersection points of the edge with the sphere for
    // later usage.
    eIntersections[n][0] =
        solutions.first[0] * direction +
        (start - element.vertices[Element::edge_mapping[n][0][0]]);

    eIntersections[n][1] =
        solutions.first[1] * direction +
        (start - element.vertices[Element::edge_mapping[n][0][1]]);

    eMarked[n] = true;

    // If the edge is marked as having an overlap, the two faces forming it
    // have to be marked as well.
    fMarked[Element::edge_mapping[n][1][0]] = true;
    fMarked[Element::edge_mapping[n][1][1]] = true;
  }

  // Check whether the dependencies for a vertex intersection are fulfilled.
  for (size_t n = 0; n < nrVertices; ++n) {
    if (!vMarked[n]) continue;

    bool edgesValid = true;
    for (size_t eN = 0; eN < 3; ++eN) {
      size_t edgeId = Element::vertex_mapping[n][0][eN];
      edgesValid &= eMarked[edgeId];
    }

    // If not all three edges intersecting at this vertex where marked, the
    // sphere is only touching.
    if (!edgesValid) vMarked[n] = false;
  }

  // Process all faces of the element, ignoring the edges as those where
  // already checked above.
  for (size_t n = 0; n < nrFaces; ++n)
    if (intersect(s, element.faces[n])) fMarked[n] = true;

  // Trivial case: The center of the sphere overlaps the element, but the
  // sphere does not intersect any of the faces of the element, meaning the
  // sphere is completely contained within the element.
  if (!fMarked.count() && contains(element, s.center)) {
    result[0] = sOrig.surfaceArea();

    return result;
  }

  // Spurious intersection: The initial intersection test was positive, but
  // the detailed checks revealed no overlap.
  if (!vMarked.count() && !eMarked.count() && !fMarked.count()) return result;

  // Initial value for the surface of the sphere: Surface area of the full
  // sphere.
  result[0] = s.surfaceArea();

  // Iterate over all the marked faces and calculate the area of the disk
  // defined by the plane as well as the cap surfaces.
  for (size_t n = 0; n < nrFaces; ++n) {
    if (!fMarked[n]) continue;

    const auto& f = element.faces[n];
    scalar_t dist = f.normal.dot(s.center - f.center);
    result[0] -= s.capSurfaceArea(s.radius + dist);
    result[n + 1] = s.diskArea(s.radius + dist);
  }

  // Handle the edges and subtract the area of the respective disk cut off by
  // the edge and add back the surface area of the spherical wedge defined
  // by the edge.
  for (size_t n = 0; n < nrEdges<Element>(); ++n) {
    if (!eMarked[n]) continue;

    result[0] += generalWedge<2, Element>(s, element, n, eIntersections);

    // The intersection points are relative to the vertices forming the
    // edge.
    const vector_t chord =
        ((element.vertices[Element::edge_mapping[n][0][0]] +
          eIntersections[n][0]) -
         (element.vertices[Element::edge_mapping[n][0][1]] +
          eIntersections[n][1]));

    const scalar_t chordLength = chord.stableNorm();

    // Each edge belongs to two faces, indexed via
    // Element::edge_mapping[n][1][{0,1}].
    for (size_t e = 0; e < 2; ++e) {
      const auto faceIdx = Element::edge_mapping[n][1][e];
      const auto& f = element.faces[faceIdx];

      // Height of the spherical cap cut off by the plane containing the
      // face.
      const scalar_t dist = f.normal.dot(s.center - f.center) + s.radius;
      intersectionRadiusSq[faceIdx] = dist * (scalar_t(2) * s.radius - dist);

      // Calculate the height of the triangular segment in the plane of
      // the base.
      const scalar_t factor = std::sqrt(std::max(
          scalar_t(0), intersectionRadiusSq[faceIdx] -
                           scalar_t(0.25) * chordLength * chordLength));

      const scalar_t theta =
          scalar_t(2) * std::atan2(chordLength, scalar_t(2) * factor);

      scalar_t area = scalar_t(0.5) * intersectionRadiusSq[faceIdx] *
                      (theta - std::sin(theta));

      // FIXME: Might not be necessary to use the center of the chord.
      const vector_t chordCenter =
          scalar_t(0.5) * ((element.vertices[Element::edge_mapping[n][0][0]] +
                            eIntersections[n][0]) +
                           (element.vertices[Element::edge_mapping[n][0][1]] +
                            eIntersections[n][1]));

      const vector_t proj(
          s.center - f.normal.dot(s.center - f.center) * f.normal);

      // If the projected sphere center and the face center fall on
      // opposite sides of the edge, the area has to be inverted.
      if (chord.cross(proj - chordCenter)
              .dot(chord.cross(f.center - chordCenter)) < scalar_t(0)) {
        area = intersectionRadiusSq[faceIdx] * pi - area;
      }

      result[faceIdx + 1] -= area;
    }
  }

  // Handle the vertices and add the area subtracted twice above in the
  // processing of the edges.

  // First, handle the spherical surface area of the intersection.
  // This is to a large part code duplicated from the volume calculation.
  // TODO: Unify the area and volume calculation to remove duplicate code.
  for (size_t n = 0; n < nrVertices; ++n) {
    if (!vMarked[n]) continue;

    // Collect the points where the three edges intersecting at this
    // vertex intersect the sphere.
    // Both the relative and the absolute positions are required.
    std::array<vector_t, 3> intersectionPointsRelative;
    std::array<vector_t, 3> intersectionPoints;
    for (size_t e = 0; e < 3; ++e) {
      auto edgeIdx = Element::vertex_mapping[n][0][e];
      intersectionPointsRelative[e] =
          eIntersections[edgeIdx][Element::vertex_mapping[n][1][e]];

      intersectionPoints[e] =
          intersectionPointsRelative[e] + element.vertices[n];
    }

    // This triangle is constructed by hand to have more freedom of how
    // the normal vector is calculated.
    Triangle coneTria;
    coneTria.vertices = {
        {intersectionPoints[0], intersectionPoints[1], intersectionPoints[2]}};

    coneTria.center = scalar_t(1.0 / 3.0) * std::accumulate(
                                                intersectionPoints.begin(),
                                                intersectionPoints.end(),
                                                vector_t::Zero().eval());

    // Calculate the normal of the triangle defined by the intersection
    // points in relative coordinates to improve accuracy.
    // Also use double the normal precision to calculate this normal.
    coneTria.normal = detail::triangleNormal(
        intersectionPointsRelative[0], intersectionPointsRelative[1],
        intersectionPointsRelative[2]);

    // The area of this triangle is never needed, so it is set to an
    // invalid value.
    coneTria.area = std::numeric_limits<scalar_t>::infinity();

    std::array<std::pair<size_t, scalar_t>, 3> distances;
    for (size_t i = 0; i < 3; ++i)
      distances[i] =
          std::make_pair(i, intersectionPointsRelative[i].squaredNorm());

    std::sort(
        distances.begin(), distances.end(),
        [](const std::pair<size_t, scalar_t>& a,
           const std::pair<size_t, scalar_t>& b) -> bool {
          return a.second < b.second;
        });

    if (distances[1].second < distances[2].second * detail::largeEpsilon) {
      // Use the general spherical wedge defined by the edge with the
      // non-degenerated intersection point and the normals of the
      // two faces forming it.
      scalar_t correction = generalWedge<2, Element>(
          s, element, Element::vertex_mapping[n][0][distances[2].first],
          eIntersections);

      result[0] -= correction;

      continue;
    }

    // Make sure the normal points in the right direction, i.e., away from
    // the center of the element.
    if (coneTria.normal.dot(element.center - coneTria.center) > scalar_t(0)) {
      coneTria.normal = -coneTria.normal;
    }

    Plane plane(coneTria.center, coneTria.normal);

    scalar_t dist = coneTria.normal.dot(s.center - coneTria.center);
    scalar_t capSurface = s.capSurfaceArea(s.radius + dist);

    // If cap surface area is small, the corrections will be even smaller.
    // There is no way to actually calculate them with reasonable
    // precision, so they are just ignored.
    if (capSurface < detail::largeEpsilon) continue;

    // Calculate the surface area of the three spherical segments between
    // the faces joining at the vertex and the plane through the
    // intersection points.
    scalar_t segmentSurface = 0;
    for (size_t e = 0; e < 3; ++e) {
      const auto& f = element.faces[Element::vertex_mapping[n][2][e]];
      uint32_t e0 = Element::face_mapping[e][0];
      uint32_t e1 = Element::face_mapping[e][1];

      vector_t center(
          scalar_t(0.5) * (intersectionPoints[e0] + intersectionPoints[e1]));

      segmentSurface += generalWedge<2>(
          s, plane, Plane(f.center, -f.normal), center - s.center);
    }

    // Calculate the surface area of the cone and clamp it to zero.
    scalar_t coneSurface = std::max(capSurface - segmentSurface, scalar_t(0));

    result[0] -= coneSurface;

    // Sanity checks: detect negative/excessively large intermediate
    // result.
    assert(result[0] > -std::sqrt(detail::tinyEpsilon));
    assert(result[0] < s.surfaceArea() + detail::tinyEpsilon);
  }

  // Second, correct the intersection area of the facets.
  for (size_t n = 0; n < nrVertices; ++n) {
    if (!vMarked[n]) continue;

    // Iterate over all the faces joining at this vertex.
    for (size_t f = 0; f < 3; ++f) {
      // Determine the two edges of this face intersecting at the
      // vertex.
      uint32_t e0 = Element::face_mapping[f][0];
      uint32_t e1 = Element::face_mapping[f][1];
      std::array<uint32_t, 2> edgeIndices = {
          {Element::vertex_mapping[n][0][e0],
           Element::vertex_mapping[n][0][e1]}};

      // Extract the (relative) intersection points of these edges with
      // the sphere furthest from the vertex.
      std::array<vector_t, 2> intersectionPoints = {
          {eIntersections[edgeIndices[0]][Element::vertex_mapping[n][1][e0]],

           eIntersections[edgeIndices[1]][Element::vertex_mapping[n][1][e1]]}};

      // Together with the vertex, this determines the triangle
      // representing one part of the correction.
      const scalar_t triaArea =
          scalar_t(0.5) *
          (intersectionPoints[0].cross(intersectionPoints[1])).stableNorm();

      // The second component is the segment defined by the face and the
      // intersection points.
      const scalar_t chordLength =
          (intersectionPoints[0] - intersectionPoints[1]).stableNorm();

      const auto faceIdx = Element::vertex_mapping[n][2][f];

      // TODO: Cache theta for each edge.
      const scalar_t theta =
          scalar_t(2) *
          std::atan2(
              chordLength,
              scalar_t(2) * std::sqrt(std::max(
                                scalar_t(0), intersectionRadiusSq[faceIdx] -
                                                 scalar_t(0.25) * chordLength *
                                                     chordLength)));

      scalar_t segmentArea = scalar_t(0.5) * intersectionRadiusSq[faceIdx] *
                             (theta - std::sin(theta));

      // Determine if the (projected) center of the sphere lies within
      // the triangle or not. If not, the segment area has to be
      // corrected.
      const vector_t d(
          scalar_t(0.5) * (intersectionPoints[0] + intersectionPoints[1]));

      const auto& face = element.faces[faceIdx];
      const vector_t proj(
          s.center - face.normal.dot(s.center - face.center) * face.normal);

      if (d.dot((proj - element.vertices[n]) - d) > scalar_t(0)) {
        segmentArea = intersectionRadiusSq[faceIdx] * pi - segmentArea;
      }

      result[faceIdx + 1] += triaArea + segmentArea;

      // Sanity checks: detect excessively large intermediate result.
      assert(
          result[faceIdx + 1] <
          element.faces[faceIdx].area + std::sqrt(detail::largeEpsilon));
    }
  }

  // Scale the surface areas back for the original objects and clamp
  // values within reasonable limits.
  const scalar_t scaling = sOrig.radius / s.radius;
  const scalar_t sLimit(
      std::sqrt(std::numeric_limits<scalar_t>::epsilon()) * s.surfaceArea());

  // As the precision of the area calculation deteriorates quickly with a
  // increasing size ratio between the element and the sphere, the precision
  // limit applied to the sphere is used as the lower limit for the facets.
  const scalar_t fLimit(std::max(
      sLimit, std::sqrt(std::numeric_limits<scalar_t>::epsilon()) *
                  element.surfaceArea()));

  // Sanity checks: detect negative/excessively large results for the
  // surface area of the facets.
#ifndef NDEBUG
  for (size_t n = 0; n < nrFaces; ++n) {
    assert(result[n + 1] > -fLimit);
    assert(result[n + 1] <= element.faces[n].area + fLimit);
  }
#endif // NDEBUG

  // Surface of the sphere.
  result[0] = detail::clamp(result[0], scalar_t(0), s.surfaceArea(), sLimit);
  result[0] *= (scaling * scaling);

  // Surface of the mesh element.
  for (size_t f = 0; f < nrFaces; ++f) {
    auto& value = result[f + 1];
    value = detail::clamp(value, scalar_t(0), element.faces[f].area, fLimit);

    value = value * (scaling * scaling);
  }

  result.back() =
      std::accumulate(result.begin() + 1, result.end() - 1, scalar_t(0));

  // Perform some more sanity checks on the final result (debug version
  // only).
  assert(scalar_t(0) <= result[0] && result[0] <= sOrig.surfaceArea());

  assert(
      scalar_t(0) <= result.back() &&
      result.back() <= elementOrig.surfaceArea());

  return result;
}
