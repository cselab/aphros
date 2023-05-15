#include <cstdint>

using uint32_t = std::uint32_t;

static Scal util_max(Scal x, Scal y) {
  return x > y ? x : y;
}

static Scal util_stable_norm(Scal x, Scal y, Scal z) {
  Scal i, j;
  x = std::abs(x);
  y = std::abs(y);
  z = std::abs(z);
  i = util_max(x, util_max(y, z));
  if (i == 0.0)
    return 0.0;
  else {
    j = 1 / i;
    x *= j;
    y *= j;
    z *= j;
    return i * std::sqrt(x * x + y * y + z * z);
  }
}

static Scal util_norm(Scal x, Scal y, Scal z) {
  return std::sqrt(x * x + y * y + z * z);
}

static Scal util_squared_norm(Scal x, Scal y, Scal z) {
  return x * x + y * y + z * z;
}

/* s . ((a - b) ~ (x - y)) */
static Scal util_diff_dot_cross(
    const Scal* s, const Scal* a, const Scal* b, const Scal* x, const Scal* y) {
  Scal u, v, w;
  Scal h, j, k;

  u = a[0] - b[0];
  v = a[1] - b[1];
  w = a[2] - b[2];
  h = x[0] - y[0];
  j = x[1] - y[1];
  k = x[2] - y[2];
  return -s[0] * (j * w - k * v) - s[1] * (k * u - h * w) -
         s[2] * (h * v - j * u);
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

class double_prec;
inline double_prec operator+(const double_prec& lhs, const double_prec& rhs)
    ENFORCE_EXACT_FPMATH_ATTR;

inline double_prec operator-(const double_prec& lhs, const double_prec& rhs)
    ENFORCE_EXACT_FPMATH_ATTR;

inline double_prec operator*(const double_prec& lhs, const double_prec& rhs)
    ENFORCE_EXACT_FPMATH_ATTR;

class double_prec {
 private:
  static const uint32_t c = double_prec_constant<Scal>::value;

  friend double_prec operator+(const double_prec&, const double_prec&);

  friend double_prec operator-(const double_prec&, const double_prec&);

  friend double_prec operator*(const double_prec&, const double_prec&);

 public:
  double_prec() : h_(0), l_(0) {}

  // This constructor requires floating point operations in accordance
  // with IEEE754 to perform the proper splitting. To allow full
  // optimization of all other parts of the code, precise floating point
  // ops are only requested here. Unfortunately the way to do this is
  // extremely compiler dependent.
  double_prec(const Scal& val) ENFORCE_EXACT_FPMATH_ATTR : h_(0), l_(0) {
#ifdef USE_EXACT_FPMATH_PRAGMA
#pragma float_control(precise, on)
#endif

    Scal p = val * Scal(c);
    h_ = (val - p) + p;
    l_ = val - h_;
  }

 private:
  explicit double_prec(const Scal& h, const Scal& l) : h_(h), l_(l) {}

 public:
  const Scal& high() const {
    return h_;
  }

  const Scal& low() const {
    return l_;
  }

  Scal value() const {
    return h_ + l_;
  }

  Scal convert() const {
    return Scal(h_) + Scal(l_);
  }

 private:
  Scal h_;
  Scal l_;
};

inline double_prec operator+(const double_prec& lhs, const double_prec& rhs) {
#ifdef USE_EXACT_FPMATH_PRAGMA
#pragma float_control(precise, on)
#endif

  Scal h = lhs.h_ + rhs.h_;
  Scal l = std::abs(lhs.h_) >= std::abs(rhs.h_)
               ? ((((lhs.h_ - h) + rhs.h_) + lhs.l_) + rhs.l_)
               : ((((rhs.h_ - h) + lhs.h_) + rhs.l_) + lhs.l_);

  Scal c = h + l;

  return double_prec(c, (h - c) + l);
}

inline double_prec operator-(const double_prec& lhs, const double_prec& rhs) {
#ifdef USE_EXACT_FPMATH_PRAGMA
#pragma float_control(precise, on)
#endif

  Scal h = lhs.h_ - rhs.h_;
  Scal l = std::abs(lhs.h_) >= std::abs(rhs.h_)
               ? ((((lhs.h_ - h) - rhs.h_) - rhs.l_) + lhs.l_)
               : ((((-rhs.h_ - h) + lhs.h_) + lhs.l_) - rhs.l_);

  Scal c = h + l;

  return double_prec(c, (h - c) + l);
}

inline double_prec operator*(const double_prec& lhs, const double_prec& rhs) {
#ifdef USE_EXACT_FPMATH_PRAGMA
#pragma float_control(precise, on)
#endif

  double_prec l(lhs.h_);
  double_prec r(rhs.h_);

  Scal p = l.h_ * r.h_;
  Scal q = l.h_ * r.l_ + l.l_ * r.h_;
  Scal v = p + q;

  double_prec c(v, ((p - v) + q) + l.l_ * r.l_);
  c.l_ = ((lhs.h_ + lhs.l_) * rhs.l_ + lhs.l_ * rhs.h_) + c.l_;
  Scal z = c.value();

  return double_prec(z, (c.h_ - z) + c.l_);
}
static const Scal pi = Scal(4) * std::atan(Scal(1));
static const Scal tinyEpsilon = 2 * std::numeric_limits<Scal>::epsilon();
static const Scal mediumEpsilon = 1e2 * tinyEpsilon;
static const Scal largeEpsilon = 1e-10;

// Ref: J.R. Shewchuk - Lecture Notes on Geometric Robustness
//      http://www.cs.berkeley.edu/~jrs/meshpapers/robnotes.pdf
static Scal orient2D(
    double_prec a0, double_prec a1, double_prec b0, double_prec b1,
    double_prec c0, double_prec c1) {
  double_prec result = (a0 - c0) * (b1 - c1) - (a1 - c1) * (b0 - c0);
  return result.convert();
}

// Numerically robust calculation of the normal of the triangle defined by
// the points a, b, and c.
// Ref: J.R. Shewchuk - Lecture Notes on Geometric Robustness
//      http://www.cs.berkeley.edu/~jrs/meshpapers/robnotes.pdf
static void triangleNormal(
    const Scal* a, const Scal* b, const Scal* c, Scal* ans) {
  ans[0] = orient2D(a[1], a[2], b[1], b[2], c[1], c[2]);
  ans[1] = orient2D(a[2], a[0], b[2], b[0], c[2], c[0]);
  ans[2] = orient2D(a[0], a[1], b[0], b[1], c[0], c[1]);
  Scal n = util_norm(ans[0], ans[1], ans[2]);
  if (n != 0.0) {
    ans[0] /= Scal(n);
    ans[1] /= Scal(n);
    ans[2] /= Scal(n);
  }
}

// Numerically robust routine to calculate the angle between normalized
// vectors.
// Ref: http://www.plunk.org/~hatch/rightway.php
static Scal vangle(const Scal* a, const Scal* b) {
  if (a[0] * b[0] + a[1] * b[1] + a[2] * b[2] < 0)
    return pi -
           2 * asin(
                   util_stable_norm(a[0] + b[0], a[1] + b[1], a[2] + b[2]) / 2);
  else
    return 2 *
           asin(util_stable_norm(a[0] - b[0], a[1] - b[1], a[2] - b[2]) / 2);
}

static void gramSchmidt(const Scal* x, const Scal* y, Scal* c) {
  Scal s;
  s = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
  s = (s == 0.0) ? 0.0 : (x[0] * y[0] + x[1] * y[1] + x[2] * y[2]) / s;
  c[0] = y[0] - s * x[0];
  c[1] = y[1] - s * x[1];
  c[2] = y[2] - s * x[2];
  s = util_norm(c[0], c[1], c[2]);
  if (s != 0.0) {
    c[0] /= s;
    c[1] /= s;
    c[2] /= s;
  }
}

class Polygon {
 public:
  Polygon(){};
  Polygon(const Scal* v0, const Scal* v1, const Scal* v2, const Scal* v3) {
    Scal n, a, b, c, x, y, z;
    int i;
    for (i = 0; i < 3; i++) {
      vertices[0][i] = v0[i];
      vertices[1][i] = v1[i];
      vertices[2][i] = v2[i];
      vertices[3][i] = v3[i];
    }
    center[0] = (v0[0] + v1[0] + v2[0] + v3[0]) / 4;
    center[1] = (v0[1] + v1[1] + v2[1] + v3[1]) / 4;
    center[2] = (v0[2] + v1[2] + v2[2] + v3[2]) / 4;
    a = v2[0] - v0[0];
    b = v2[1] - v0[1];
    c = v2[2] - v0[2];
    x = v3[0] - v1[0];
    y = v3[1] - v1[1];
    z = v3[2] - v1[2];
    normal[0] = b * z - c * y;
    normal[1] = c * x - a * z;
    normal[2] = a * y - b * x;
    n = util_norm(normal[0], normal[1], normal[2]);
    if (n != 0.0) {
      normal[0] /= n;
      normal[1] /= n;
      normal[2] /= n;
    }
  }

  void apply(const Scal* translation, Scal scaling) {
    size_t i;
    for (i = 0; i < sizeof vertices / sizeof *vertices; i++) {
      vertices[i][0] = scaling * (vertices[i][0] + translation[0]);
      vertices[i][1] = scaling * (vertices[i][1] + translation[1]);
      vertices[i][2] = scaling * (vertices[i][2] + translation[2]);
    }
    center[0] = scaling * (center[0] + translation[0]);
    center[1] = scaling * (center[1] + translation[1]);
    center[2] = scaling * (center[2] + translation[2]);
  }

 public:
  Scal vertices[4][3];
  Scal center[3];
  Scal normal[3];
};

enum { nrEdges = 12, nrVertices = 8, nrFaces = 6 };
static const uint32_t edge_mapping[nrEdges][2][2] = {
    {{0, 1}, {0, 1}}, {{1, 2}, {0, 2}}, {{2, 3}, {0, 3}}, {{3, 0}, {0, 4}},

    {{0, 4}, {1, 4}}, {{1, 5}, {1, 2}}, {{2, 6}, {2, 3}}, {{3, 7}, {3, 4}},

    {{4, 5}, {1, 5}}, {{5, 6}, {2, 5}}, {{6, 7}, {3, 5}}, {{7, 4}, {4, 5}}};

static const uint32_t vertex_mapping[nrVertices][3][3] = {
    {{0, 3, 4}, {0, 1, 0}, {0, 1, 4}},  {{0, 1, 5}, {1, 0, 0}, {0, 1, 2}},
    {{1, 2, 6}, {1, 0, 0}, {0, 2, 3}},  {{2, 3, 7}, {1, 0, 0}, {0, 3, 4}},

    {{4, 8, 11}, {1, 0, 1}, {1, 4, 5}}, {{5, 8, 9}, {1, 1, 0}, {1, 2, 5}},
    {{6, 9, 10}, {1, 1, 0}, {2, 3, 5}}, {{7, 10, 11}, {1, 1, 0}, {3, 4, 5}}};

static const uint32_t face_mapping[3][2] = {{0, 1}, {0, 2}, {1, 2}};

class Hexahedron {
 public:
  Hexahedron(
      const Scal* v0, const Scal* v1, const Scal* v2, const Scal* v3,
      const Scal* v4, const Scal* v5, const Scal* v6, const Scal* v7)
      : faces{Polygon(v3, v2, v1, v0), Polygon(v0, v1, v5, v4),
              Polygon(v1, v2, v6, v5), Polygon(v2, v3, v7, v6),
              Polygon(v3, v0, v4, v7), Polygon(v4, v5, v6, v7)} {
    int i;
    for (i = 0; i < 3; i++) {
      vertices[0][i] = v0[i];
      vertices[1][i] = v1[i];
      vertices[2][i] = v2[i];
      vertices[3][i] = v3[i];
      vertices[4][i] = v4[i];
      vertices[5][i] = v5[i];
      vertices[6][i] = v6[i];
      vertices[7][i] = v7[i];
    }
    center[0] =
        (v0[0] + v1[0] + v2[0] + v3[0] + v4[0] + v5[0] + v6[0] + v7[0]) / 8;
    center[1] =
        (v0[1] + v1[1] + v2[1] + v3[1] + v4[1] + v5[1] + v6[1] + v7[1]) / 8;
    center[2] =
        (v0[2] + v1[2] + v2[2] + v3[2] + v4[2] + v5[2] + v6[2] + v7[2]) / 8;
    volume = calcVolume();
  }

  void apply(const Scal* translation, Scal scaling) {
    int i;
    for (i = 0; i < nrVertices; i++) {
      vertices[i][0] = scaling * (vertices[i][0] + translation[0]);
      vertices[i][1] = scaling * (vertices[i][1] + translation[1]);
      vertices[i][2] = scaling * (vertices[i][2] + translation[2]);
    }
    for (i = 0; i < nrFaces; i++)
      faces[i].apply(translation, scaling);
    for (i = 0; i < 3; i++)
      center[i] =
          (vertices[0][i] + vertices[1][i] + vertices[2][i] + vertices[3][i] +
           vertices[4][i] + vertices[5][i] + vertices[6][i] + vertices[7][i]) /
          8;
    volume = calcVolume();
  }

  Scal calcVolume() const {
    Scal d[3];
    d[0] = vertices[6][0] - vertices[0][0];
    d[1] = vertices[6][1] - vertices[0][1];
    d[2] = vertices[6][2] - vertices[0][2];
    return (util_diff_dot_cross(
                d, vertices[1], vertices[0], vertices[2], vertices[5]) +
            util_diff_dot_cross(
                d, vertices[4], vertices[0], vertices[5], vertices[7]) +
            util_diff_dot_cross(
                d, vertices[3], vertices[0], vertices[7], vertices[2])) /
           6;
  }

 public:
  Scal vertices[nrVertices][3];
  Polygon faces[nrFaces];
  Scal center[3];
  Scal volume;
};

class Sphere {
 public:
  Sphere(Scal x, Scal y, Scal z, Scal r)
      : radius(r), volume(Scal(4.0 / 3.0 * pi) * r * r * r) {
    center[0] = x;
    center[1] = y;
    center[2] = z;
  }

  Scal capVolume(Scal h) const {
    if (h <= Scal(0))
      return Scal(0);
    else if (h >= Scal(2) * radius)
      return volume;
    else
      return Scal(pi / 3.0) * h * h * (Scal(3) * radius - h);
  }

 public:
  Scal center[3];
  Scal radius;
  Scal volume;
};

class Plane {
 public:
  Plane(const Scal* c, const Scal* n) {
    center[0] = c[0];
    center[1] = c[1];
    center[2] = c[2];
    normal[0] = n[0];
    normal[1] = n[1];
    normal[2] = n[2];
  }

 public:
  Scal center[3];
  Scal normal[3];
};

class AABB {
 public:
  AABB() {
    min[0] = min[1] = min[2] = std::numeric_limits<Scal>::infinity();
    max[0] = max[1] = max[2] = -std::numeric_limits<Scal>::infinity();
  }
  AABB(const Scal* mi, const Scal* ma) {
    min[0] = mi[0];
    min[1] = mi[1];
    min[2] = mi[2];
    max[0] = ma[0];
    max[1] = ma[1];
    max[2] = ma[2];
  }

  bool intersects(const AABB* aabb) const {
    return !(
        min[0] > aabb->max[0] || min[1] > aabb->max[1] ||
        min[2] > aabb->max[2] || max[0] < aabb->min[0] ||
        max[1] < aabb->min[1] || max[2] < aabb->min[2]);
  }

  void include(int n, const Scal (*p)[3]) {
    int i;
    for (i = 0; i < n; i++) {
      if (p[i][0] < min[0]) min[0] = p[i][0];
      if (p[i][1] < min[1]) min[1] = p[i][1];
      if (p[i][2] < min[2]) min[2] = p[i][2];
      if (p[i][0] > max[0]) max[0] = p[i][0];
      if (p[i][1] > max[1]) max[1] = p[i][1];
      if (p[i][2] > max[2]) max[2] = p[i][2];
    }
  }

 public:
  Scal min[3];
  Scal max[3];
};

// The (convex!) polygon is assumed to be planar, making this a 2D problem.
// Check the projection of the point onto the plane of the polygon for
// containment within the polygon.
static bool contains(const Polygon* poly, const Scal* point) {
  Scal proj[3];
  Scal dot = poly->normal[0] * (point[0] - poly->center[0]) +
             poly->normal[1] * (point[1] - poly->center[1]) +
             poly->normal[2] * (point[2] - poly->center[2]);
  proj[0] = point[0] - dot * poly->normal[0];
  proj[1] = point[1] - dot * poly->normal[1];
  proj[2] = point[2] - dot * poly->normal[2];

  const int nsize = sizeof poly->vertices / sizeof poly->vertices[0];
  for (size_t n = 0; n < nsize; ++n) {
    const Scal* v0 = poly->vertices[n];
    const Scal* v1 = poly->vertices[(n + 1) % nsize];
    Scal base[3];
    Scal edge[3];
    edge[0] = v1[0] - v0[0];
    edge[1] = v1[1] - v0[1];
    edge[2] = v1[2] - v0[2];

    base[0] = (v0[0] + v1[0]) / Scal(2);
    base[1] = (v0[1] + v1[1]) / Scal(2);
    base[2] = (v0[2] + v1[2]) / Scal(2);

    // Note: Only the sign of the projection is of interest, so this vector
    // does not have to be normalized.
    Scal s;
    s = (edge[1] * poly->normal[2] - poly->normal[1] * edge[2]) *
            (proj[0] - base[0]) +
        (poly->normal[0] * edge[2] - edge[0] * poly->normal[2]) *
            (proj[1] - base[1]) +
        (edge[0] * poly->normal[1] - poly->normal[0] * edge[1]) *
            (proj[2] - base[2]);

    // Check whether the projection of the point lies inside of the
    // polygon.
    if (s > 0) return false;
  }

  return true;
}

static bool contains(const Hexahedron* hex, const Scal* p) {
  size_t i;
  Scal s;
  for (i = 0; i < sizeof hex->faces / sizeof *hex->faces; i++) {
    s = hex->faces[i].normal[0] * (p[0] - hex->faces[i].center[0]) +
        hex->faces[i].normal[1] * (p[1] - hex->faces[i].center[1]) +
        hex->faces[i].normal[2] * (p[2] - hex->faces[i].center[2]);
    if (s > 0) return false;
  }
  return true;
}

static bool intersect0(
    const Sphere* s, const Scal* center, const Scal* normal) {
  Scal proj = normal[0] * (s->center[0] - center[0]) +
              normal[1] * (s->center[1] - center[1]) +
              normal[2] * (s->center[2] - center[2]);
  return proj * proj - s->radius * s->radius < 0;
}

static bool intersect(const Sphere* s, const Polygon* poly) {
  return intersect0(s, poly->center, poly->normal) && contains(poly, s->center);
}

void lineSphereIntersection(
    const Scal* origin, const Scal* direction, const Sphere* s, int* nsol,
    Scal* sol) {
  sol[0] = sol[1] = std::numeric_limits<Scal>::infinity();
  Scal originRel[3];
  originRel[0] = origin[0] - s->center[0];
  originRel[1] = origin[1] - s->center[1];
  originRel[2] = origin[2] - s->center[2];
  Scal a = util_squared_norm(direction[0], direction[1], direction[2]);
  if (a == 0.0) {
    *nsol = 0;
    return;
  }
  Scal b = 2 * (direction[0] * originRel[0] + direction[1] * originRel[1] +
                direction[2] * originRel[2]);
  Scal c = util_squared_norm(originRel[0], originRel[1], originRel[2]) -
           s->radius * s->radius;

  Scal discriminant = b * b - Scal(4) * a * c;
  if (discriminant > Scal(0)) {
    // Two real roots.
    Scal q = Scal(-0.5) * (b + std::copysign(std::sqrt(discriminant), b));

    sol[0] = q / a;
    sol[1] = c / q;

    if (sol[0] > sol[1]) std::swap(sol[0], sol[1]);
    *nsol = 2;
  } else if (std::abs(discriminant) == Scal(0)) {
    // Scal real root.
    sol[0] = (Scal(-0.5) * b) / a;
    sol[1] = sol[0];
    *nsol = 1;
  } else {
    // No real roots.
    *nsol = 0;
  }
}

// Calculate the volume of a regularized spherical wedge defined by the radius,
// the distance of the intersection point from the center of the sphere and the
// angle.
static Scal regularizedWedge(Scal r, Scal d, Scal alpha) {
#ifndef NDEBUG
  // Clamp slight deviations of the angle to valid range.
  if (alpha < Scal(0) && alpha > -tinyEpsilon) alpha = Scal(0);

  if (alpha > Scal(0.5 * pi) && alpha < Scal(0.5 * pi) + tinyEpsilon)
    alpha = Scal(0.5 * pi);
#endif

  // Check the parameters for validity (debug version only).
  assert(r > Scal(0));
  assert(d >= Scal(0) && d <= r);
  assert(alpha >= Scal(0) && alpha <= Scal(0.5 * pi));

  const Scal sinAlpha = std::sin(alpha);
  const Scal cosAlpha = std::cos(alpha);

  const Scal a = d * sinAlpha;
  const Scal b = std::sqrt(std::abs(r * r - d * d));
  const Scal c = d * cosAlpha;

  return Scal(1.0 / 3.0) * a * b * c +
         a * (Scal(1.0 / 3.0) * a * a - r * r) * std::atan2(b, c) +
         Scal(2.0 / 3.0) * r * r * r * std::atan2(sinAlpha * b, cosAlpha * r);
}

// Wrapper around the above function handling correctly handling the case of
// alpha > pi/2 and negative z.
static Scal regularizedWedge(Scal r, Scal d, Scal alpha, Scal z) {
  if (z >= Scal(0)) {
    if (alpha > Scal(0.5 * pi)) {
      Scal h = r - z;

      return Scal(pi / 3.0) * h * h * (Scal(3) * r - h) -
             regularizedWedge(r, d, pi - alpha);
    } else {
      return regularizedWedge(r, d, alpha);
    }
  } else {
    Scal vHem = Scal(2.0 / 3.0 * pi) * r * r * r;

    if (alpha > Scal(0.5 * pi)) {
      return vHem - regularizedWedge(r, d, pi - alpha);
    } else {
      Scal h = r + z;
      Scal vCap = Scal(pi / 3.0) * h * h * (Scal(3) * r - h);

      return vHem - (vCap - regularizedWedge(r, d, alpha));
    }
  }
}

static Scal generalWedge(
    const Sphere& s, const Plane& p0, const Plane& p1, const Scal* d) {
  Scal dist = util_stable_norm(d[0], d[1], d[2]);
  if (dist < tinyEpsilon) {
    // The wedge (almost) touches the center, the volume depends only on
    // the angle.
    Scal angle = pi - vangle(p0.normal, p1.normal);
    return Scal(2.0 / 3.0) * s.radius * s.radius * s.radius * angle;
  }

  Scal s0 = d[0] * p0.normal[0] + d[1] * p0.normal[1] + d[2] * p0.normal[2];
  Scal s1 = d[0] * p1.normal[0] + d[1] * p1.normal[1] + d[2] * p1.normal[2];

  // Detect degenerated general spherical wedge that can be treated as
  // a regularized spherical wedge.
  if (std::abs(s0) < tinyEpsilon || std::abs(s1) < tinyEpsilon) {
    Scal angle = pi - vangle(p0.normal, p1.normal);

    return regularizedWedge(
        s.radius, dist, angle, std::abs(s0) > std::abs(s1) ? s0 : s1);
  }

  Scal dUnit[3];
  dUnit[0] = d[0] * (Scal(1) / dist);
  dUnit[1] = d[1] * (Scal(1) / dist);
  dUnit[2] = d[2] * (Scal(1) / dist);

  if (dist < largeEpsilon) {
    Scal cross[3];
    Scal dUnit0[3];
    cross[0] = p0.normal[1] * p1.normal[2] - p1.normal[1] * p0.normal[2];
    cross[1] = p1.normal[0] * p0.normal[2] - p0.normal[0] * p1.normal[2];
    cross[2] = p0.normal[0] * p1.normal[1] - p1.normal[0] * p0.normal[1];
    gramSchmidt(cross, dUnit, dUnit0);
    dUnit[0] = dUnit0[0];
    dUnit[1] = dUnit0[1];
    dUnit[2] = dUnit0[2];
  }

  // Check the planes specify a valid setup (debug version only).
  assert(
      p0.normal[0] * (p1.center[0] - p0.center[0]) +
          p0.normal[1] * (p1.center[1] - p0.center[1]) +
          p0.normal[2] * (p1.center[2] - p0.center[2]) <=
      Scal(0));
  assert(
      p1.normal[0] * (p0.center[0] - p1.center[0]) +
          p1.normal[1] * (p0.center[1] - p1.center[1]) +
          p1.normal[2] * (p0.center[2] - p1.center[2]) <=
      Scal(0));

  // Calculate the angles between the vector from the sphere center
  // to the intersection line and the normal vectors of the two planes.
  Scal alpha0 = vangle(p0.normal, dUnit);
  Scal alpha1 = vangle(p1.normal, dUnit);
  Scal dir0 = dUnit[0] * (s.center[0] + d[0] - p0.center[0]) +
              dUnit[1] * (s.center[1] + d[1] - p0.center[1]) +
              dUnit[2] * (s.center[2] + d[2] - p0.center[2]);
  Scal dir1 = dUnit[0] * (s.center[0] + d[0] - p1.center[0]) +
              dUnit[1] * (s.center[1] + d[1] - p1.center[1]) +
              dUnit[2] * (s.center[2] + d[2] - p1.center[2]);
  if (s0 >= Scal(0) && s1 >= Scal(0)) {
    alpha0 = Scal(0.5 * pi) - std::copysign(alpha0, dir0);
    alpha1 = Scal(0.5 * pi) - std::copysign(alpha1, dir1);

    return regularizedWedge(s.radius, dist, alpha0, s0) +
           regularizedWedge(s.radius, dist, alpha1, s1);
  } else if (s0 < Scal(0) && s1 < Scal(0)) {
    alpha0 = Scal(0.5 * pi) + std::copysign(Scal(1), dir0) * (alpha0 - pi);

    alpha1 = Scal(0.5 * pi) + std::copysign(Scal(1), dir1) * (alpha1 - pi);

    return s.volume - (regularizedWedge(s.radius, dist, alpha0, -s0) +
                       regularizedWedge(s.radius, dist, alpha1, -s1));
  } else {
    alpha0 = Scal(0.5 * pi) - std::copysign(Scal(1), dir0 * s0) *
                                  (alpha0 - (s0 < Scal(0) ? pi : Scal(0)));

    alpha1 = Scal(0.5 * pi) - std::copysign(Scal(1), dir1 * s1) *
                                  (alpha1 - (s1 < Scal(0) ? pi : Scal(0)));

    Scal volume0 = regularizedWedge(s.radius, dist, alpha0, std::abs(s0));

    Scal volume1 = regularizedWedge(s.radius, dist, alpha1, std::abs(s1));

    return std::max(volume0, volume1) - std::min(volume0, volume1);
  }
}

// Depending on the dimensionality, either the volume or external surface area
// of the general wedge is computed.
static Scal generalWedge(
    const Sphere& sphere, const Hexahedron& element, size_t edge,
    const Scal intersections[][2][3]) {
  const Polygon* f0 = &element.faces[edge_mapping[edge][1][0]];
  const Polygon* f1 = &element.faces[edge_mapping[edge][1][1]];

  Scal edgeCenter[3];
  for (size_t i = 0; i < 3; i++)
    edgeCenter[i] = 0.5 * (intersections[edge][0][i] +
                           element.vertices[edge_mapping[edge][0][0]][i] +
                           intersections[edge][1][i] +
                           element.vertices[edge_mapping[edge][0][1]][i]);

  Plane p0(f0->center, f0->normal);
  Plane p1(f1->center, f1->normal);

  Scal diff[3];
  diff[0] = edgeCenter[0] - sphere.center[0];
  diff[1] = edgeCenter[1] - sphere.center[1];
  diff[2] = edgeCenter[2] - sphere.center[2];
  return generalWedge(sphere, p0, p1, diff);
}

static Scal overlap(const Sphere& sOrig, const Hexahedron& elementOrig) {
  // Construct AABBs and perform a coarse overlap detection.
  Scal a[3];
  Scal b[3];
  a[0] = sOrig.center[0] - sOrig.radius;
  a[1] = sOrig.center[1] - sOrig.radius;
  a[2] = sOrig.center[2] - sOrig.radius;
  b[0] = sOrig.center[0] + sOrig.radius;
  b[1] = sOrig.center[1] + sOrig.radius;
  b[2] = sOrig.center[2] + sOrig.radius;
  AABB sAABB(a, b);

  AABB eAABB;
  eAABB.include(
      sizeof elementOrig.vertices / sizeof *elementOrig.vertices,
      &elementOrig.vertices[0]);

  if (!sAABB.intersects(&eAABB)) return Scal(0);

  // Use scaled and shifted versions of the sphere and the element.
  Scal translation[3];
  translation[0] = -sOrig.center[0];
  translation[1] = -sOrig.center[1];
  translation[2] = -sOrig.center[2];
  Scal scaling = Scal(1) / sOrig.radius;

  Sphere s(Scal(0), Scal(0), Scal(0), Scal(1));

  Hexahedron element(elementOrig);
  element.apply(translation, scaling);

  size_t vOverlap = 0;
  // Check whether the vertices lie on or outside of the sphere.
  for (size_t i = 0; i < sizeof element.vertices / sizeof *element.vertices;
       i++) {
    const Scal* vertex = element.vertices[i];
    if (util_squared_norm(
            s.center[0] - vertex[0], s.center[1] - vertex[1],
            s.center[2] - vertex[2]) <= s.radius * s.radius)
      ++vOverlap;
  }

  // Check for trivial case: All vertices inside of the sphere, resulting in
  // a full overlap.
  if (vOverlap == nrVertices) return elementOrig.volume;

  // Sets of overlapping primitives.
  std::bitset<nrVertices> vMarked;
  std::bitset<nrEdges> eMarked;
  std::bitset<nrFaces> fMarked;

  // Initial value: Volume of the full sphere.
  Scal result = s.volume;

  // The intersection points between the single edges and the sphere, this
  // is needed later on.
  Scal eIntersections[nrEdges][2][3];

  // Process all edges of the element.
  for (size_t n = 0; n < nrEdges; ++n) {
    Scal start[3];
    Scal direction[3];
    start[0] = element.vertices[edge_mapping[n][0][0]][0];
    start[1] = element.vertices[edge_mapping[n][0][0]][1];
    start[2] = element.vertices[edge_mapping[n][0][0]][2];

    direction[0] = element.vertices[edge_mapping[n][0][1]][0] - start[0];
    direction[1] = element.vertices[edge_mapping[n][0][1]][1] - start[1];
    direction[2] = element.vertices[edge_mapping[n][0][1]][2] - start[2];

    int nsol;
    Scal sol[2];
    lineSphereIntersection(start, direction, &s, &nsol, sol);

    // No intersection between the edge and the sphere, where intersection
    // points close to the surface of the sphere are ignored.
    // Or:
    // The sphere cuts the edge twice, no vertex is inside of the
    // sphere, but the case of the edge only touching the sphere has to
    // be avoided.
    if (!nsol || (sol[0] >= Scal(1) - mediumEpsilon) ||
        sol[1] <= mediumEpsilon ||
        (sol[0] > Scal(0) && sol[1] < Scal(1) &&
         (sol[1] - sol[0] < largeEpsilon))) {
      continue;
    } else {
      vMarked[edge_mapping[n][0][0]] = sol[0] < Scal(0);
      vMarked[edge_mapping[n][0][1]] = sol[1] > Scal(1);
    }

    // Store the two intersection points of the edge with the sphere for
    // later usage.
    for (int i = 0; i < 3; i++) {
      eIntersections[n][0][i] =
          sol[0] * direction[i] +
          (start[i] - element.vertices[edge_mapping[n][0][0]][i]);
      eIntersections[n][1][i] =
          sol[1] * direction[i] +
          (start[i] - element.vertices[edge_mapping[n][0][1]][i]);
    }

    eMarked[n] = true;

    // If the edge is marked as having an overlap, the two faces forming it
    // have to be marked as well.
    fMarked[edge_mapping[n][1][0]] = true;
    fMarked[edge_mapping[n][1][1]] = true;
  }

  // Check whether the dependencies for a vertex intersection are fulfilled.
  for (size_t n = 0; n < nrVertices; ++n) {
    if (!vMarked[n]) continue;

    bool edgesValid = true;
    for (size_t eN = 0; eN < 3; ++eN) {
      size_t edgeId = vertex_mapping[n][0][eN];
      edgesValid &= eMarked[edgeId];
    }

    // If not all three edges intersecting at this vertex where marked, the
    // sphere is only touching.
    if (!edgesValid) vMarked[n] = false;
  }

  // Process all faces of the element, ignoring the edges as those where
  // already checked above.
  for (size_t n = 0; n < nrFaces; ++n)
    if (intersect(&s, &element.faces[n])) fMarked[n] = true;

  // Trivial case: The center of the sphere overlaps the element, but the
  // sphere does not intersect any of the faces of the element, meaning the
  // sphere is completely contained within the element.
  if (!fMarked.count() && contains(&element, s.center)) return sOrig.volume;

  // Spurious intersection: The initial intersection test was positive, but
  // the detailed checks revealed no overlap.
  if (!vMarked.count() && !eMarked.count() && !fMarked.count()) return Scal(0);

  // Iterate over all the marked faces and subtract the volume of the cap cut
  // off by the plane.
  for (size_t n = 0; n < nrFaces; ++n) {
    if (!fMarked[n]) continue;

    const Polygon* f = &element.faces[n];
    Scal dist = f->normal[0] * (s.center[0] - f->center[0]) +
                f->normal[1] * (s.center[1] - f->center[1]) +
                f->normal[2] * (s.center[2] - f->center[2]);

    Scal vCap = s.capVolume(s.radius + dist);

    result -= vCap;
  }

  // Handle the edges and add back the volume subtracted twice above in the
  // processing of the faces.
  for (size_t n = 0; n < nrEdges; ++n) {
    if (!eMarked[n]) continue;

    Scal edgeCorrection = generalWedge(s, element, n, eIntersections);

    result += edgeCorrection;
  }

  // Handle the vertices and subtract the volume added twice above in the
  // processing of the edges.
  for (size_t n = 0; n < nrVertices; ++n) {
    if (!vMarked[n]) continue;

    // Collect the points where the three edges intersecting at this
    // vertex intersect the sphere.
    // Both the relative and the absolute positions are required.
    Scal ip[3][3];
    Scal intersectionPoints[3][3];
    for (size_t e = 0; e < 3; ++e) {
      auto edgeIdx = vertex_mapping[n][0][e];
      ip[e][0] = eIntersections[edgeIdx][vertex_mapping[n][1][e]][0];
      ip[e][1] = eIntersections[edgeIdx][vertex_mapping[n][1][e]][1];
      ip[e][2] = eIntersections[edgeIdx][vertex_mapping[n][1][e]][2];
      intersectionPoints[e][0] = ip[e][0] + element.vertices[n][0];
      intersectionPoints[e][1] = ip[e][1] + element.vertices[n][1];
      intersectionPoints[e][2] = ip[e][2] + element.vertices[n][2];
    }

    // This triangle is constructed by hand to have more freedom of how
    // the normal vector is calculated.

    Scal center[3];
    center[0] =
        Scal(1.0 / 3.0) * (intersectionPoints[0][0] + intersectionPoints[1][0] +
                           intersectionPoints[2][0]);
    center[1] =
        Scal(1.0 / 3.0) * (intersectionPoints[0][1] + intersectionPoints[1][1] +
                           intersectionPoints[2][1]);
    center[2] =
        Scal(1.0 / 3.0) * (intersectionPoints[0][2] + intersectionPoints[1][2] +
                           intersectionPoints[2][2]);

    // Calculate the normal of the triangle defined by the intersection
    // points in relative coordinates to improve accuracy.
    // Also use Scal the normal precision to calculate this normal.
    Scal normal[3];
    triangleNormal(ip[0], ip[1], ip[2], normal);

    struct {
      Scal d;
      int i;
    } d0, d1, d2, dt;
    d0.d = util_squared_norm(ip[0][0], ip[0][1], ip[0][2]);
    d1.d = util_squared_norm(ip[1][0], ip[1][1], ip[1][2]);
    d2.d = util_squared_norm(ip[2][0], ip[2][1], ip[2][2]);
    d0.i = 0;
    d1.i = 1;
    d2.i = 2;
    if (d0.d > d2.d) {
      dt = d0; d0 = d2; d2 = dt;
    }
    if (d0.d > d1.d) {
      dt = d0; d0 = d1; d1 = dt;
    }
    if (d1.d > d2.d) {
      dt = d1; d1 = d2; d2 = dt;
    }
    assert(d0.d <= d1.d && d1.d <= d2.d);
    if (d1.d < d2.d * largeEpsilon) {
      // Use the general spherical wedge defined by the edge with the
      // non-degenerated intersection point and the normals of the
      // two faces forming it.
      Scal correction = generalWedge(
          s, element, vertex_mapping[n][0][d2.i],
          eIntersections);

      result -= correction;

      continue;
    }

    Scal tipTetVolume =
        Scal(1.0 / 6.0) *
        std::abs(util_diff_dot_cross(ip[2], ip[0], ip[2], ip[1], ip[2]));

    // Make sure the normal points in the right direction i.e. away from
    // the center of the element.
    if (normal[0] * (element.center[0] - center[0]) +
            normal[1] * (element.center[1] - center[1]) +
            normal[2] * (element.center[2] - center[2]) >
        0) {
      normal[0] = -normal[0];
      normal[1] = -normal[1];
      normal[2] = -normal[2];
    }
    Plane plane(center, normal);

    Scal dist = normal[0] * (s.center[0] - center[0]) +
                normal[1] * (s.center[1] - center[1]) +
                normal[2] * (s.center[2] - center[2]);
    Scal capVolume = s.capVolume(s.radius + dist);

    // The cap volume is tiny, so the corrections will be even smaller.
    // There is no way to actually calculate them with reasonable
    // precision, so just the volume of the tetrahedron at the tip is
    // used.
    if (capVolume < tinyEpsilon) {
      result -= tipTetVolume;
      continue;
    }

    // Calculate the volume of the three spherical segments between
    // the faces joining at the vertex and the plane through the
    // intersection points.
    Scal segmentVolume = 0;

    for (size_t e = 0; e < 3; ++e) {
      const auto& f = element.faces[vertex_mapping[n][2][e]];
      uint32_t e0 = face_mapping[e][0];
      uint32_t e1 = face_mapping[e][1];

      Scal c[3];
      c[0] =
          (intersectionPoints[e0][0] + intersectionPoints[e1][0]) / Scal(2);
      c[1] =
          (intersectionPoints[e0][1] + intersectionPoints[e1][1]) / Scal(2);
      c[2] =
          (intersectionPoints[e0][2] + intersectionPoints[e1][2]) / Scal(2);

      Scal neg_normal[3];
      Scal diff[3];
      neg_normal[0] = -f.normal[0];
      neg_normal[1] = -f.normal[1];
      neg_normal[2] = -f.normal[2];
      diff[0] = c[0] - s.center[0];
      diff[1] = c[1] - s.center[1];
      diff[2] = c[2] - s.center[2];
      Scal wedgeVolume =
          generalWedge(s, plane, Plane(f.center, neg_normal), diff);

      segmentVolume += wedgeVolume;
    }

    // Calculate the volume of the cone and clamp it to zero.
    Scal coneVolume =
        std::max(tipTetVolume + capVolume - segmentVolume, Scal(0));

    // Sanity check: detect negative cone volume.
    assert(coneVolume > -std::sqrt(tinyEpsilon));

    result -= coneVolume;

    // Sanity check: detect negative intermediate result.
    assert(result > -std::sqrt(tinyEpsilon));
  }

  // In case of different sized objects the error can become quite large,
  // so a relative limit is used.
  Scal maxOverlap = std::min(s.volume, element.volume);
  const Scal limit(
      std::sqrt(std::numeric_limits<Scal>::epsilon()) * maxOverlap);

  // Clamp tiny negative volumes to zero.
  if (result < Scal(0) && result > -limit) return Scal(0);

  // Clamp results slightly too large.
  if (result > maxOverlap && result - maxOverlap < limit)
    return std::min(sOrig.volume, elementOrig.volume);

  // Perform a sanity check on the final result (debug version only).
  assert(result >= Scal(0) && result <= maxOverlap);

  // Scale the overlap volume back for the original objects.
  result = (result / s.volume) * sOrig.volume;

  return result;
}
/*!
 * Exact calculation of the overlap volume of spheres and mesh elements.
 * http://dx.doi.org/10.1016/j.jcp.2016.02.003
 *
 * Copyright (C) 2015-2017 Severin Strobl <severin.strobl@fau.de>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */
