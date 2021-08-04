static const double pi = 3.141592653589793;
static const double tinyEpsilon = 2 * DBL_EPSILON;
static const double mediumEpsilon = 1e2 * tinyEpsilon;
static const double largeEpsilon = 1e-10;

// Ref: J.R. Shewchuk - Lecture Notes on Geometric Robustness
//      http://www.cs.berkeley.edu/~jrs/meshpapers/robnotes.pdf
typedef double_prec<double> real_t;
inline double orient2D(
    real_t a0, real_t a1, real_t b0, real_t b1, real_t c0, real_t c1) {
  real_t result = (a0 - c0) * (b1 - c1) - (a1 - c1) * (b0 - c0);
  return result.convert<double>();
}

// Numerically robust calculation of the normal of the triangle defined by
// the points a, b, and c.
// Ref: J.R. Shewchuk - Lecture Notes on Geometric Robustness
//      http://www.cs.berkeley.edu/~jrs/meshpapers/robnotes.pdf
inline vector_t triangleNormal(
    const vector_t& a, const vector_t& b, const vector_t& c) {
  double xy = orient2D(a[0], a[1], b[0], b[1], c[0], c[1]);
  double yz = orient2D(a[1], a[2], b[1], b[2], c[1], c[2]);
  double zx = orient2D(a[2], a[0], b[2], b[0], c[2], c[0]);
  double n = util_norm(yz, zx, xy);
  return n == 0 ? vector_t(yz, zx, xy) : vector_t(yz/n, zx/n, xy/n);
}

// Numerically robust routine to calculate the angle between normalized
// vectors.
// Ref: http://www.plunk.org/~hatch/rightway.php
inline double vangle(const vector_t& a, const vector_t& b) {
  if (a[0]*b[0] + a[1]*b[1] + a[2]*b[2] < 0)
    return pi - 2 * asin(util_stable_norm(a[0] + b[0],
					  a[1] + b[1],
					  a[2] + b[2])/2);
  else
    return 2 * asin(util_stable_norm(a[0] - b[0],
				     a[1] - b[1],
				     a[2] - b[2])/2);
}

inline vector_t gramSchmidt(const vector_t& x, const vector_t& y) {
  double s;
  vector_t c;
  s = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
  s = (s == 0.0) ? 0.0 : (x[0]*y[0] + x[1]*y[1] + x[2]*y[2])/s;
  c[0] = y[0] - s*x[0];
  c[1] = y[1] - s*x[1];
  c[2] = y[2] - s*x[2];
  s = util_norm(c[0], c[1], c[2]);
  if (s != 0.0) {
    c[0] /= s;
    c[1] /= s;
    c[2] /= s;
  }
  return c;
}

class Polygon {
 public:
  Polygon() { };
  Polygon(
      const vector_t& v0, const vector_t& v1, const vector_t& v2,
      const vector_t& v3)
      : vertices{v0, v1, v2, v3}, center(), normal() {
    double n, a, b, c, x, y, z;
    center[0] = (v0[0] + v1[0] + v2[0] + v3[0])/4;
    center[1] = (v0[1] + v1[1] + v2[1] + v3[1])/4;
    center[2] = (v0[2] + v1[2] + v2[2] + v3[2])/4;
    a = v2[0] - v0[0];
    b = v2[1] - v0[1];
    c = v2[2] - v0[2];
    x = v3[0] - v1[0];
    y = v3[1] - v1[1];
    z = v3[2] - v1[2];
    normal[0] = b*z-c*y;
    normal[1] = c*x-a*z;
    normal[2] = a*y-b*x;
    n = util_norm(normal[0], normal[1], normal[2]);
    if (n != 0.0) {
      normal[0] /= n;
      normal[1] /= n;
      normal[2] /= n;
    }
  }

  void apply(const vector_t& translation, double scaling) {
    size_t i;
    for (i = 0; i < sizeof vertices/sizeof *vertices; i++) {
      vertices[i][0] = scaling * (vertices[i][0] + translation[0]);
      vertices[i][1] = scaling * (vertices[i][1] + translation[1]);
      vertices[i][2] = scaling * (vertices[i][2] + translation[2]);
    }
    center[0] = scaling * (center[0] + translation[0]);
    center[1] = scaling * (center[1] + translation[1]);
    center[2] = scaling * (center[2] + translation[2]);
  }

 public:
  vector_t vertices[4];
  vector_t center;
  vector_t normal;
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
      const vector_t& v0, const vector_t& v1, const vector_t& v2,
      const vector_t& v3, const vector_t& v4, const vector_t& v5,
      const vector_t& v6, const vector_t& v7)
    : vertices{v0, v1, v2, v3, v4, v5, v6, v7},
      faces{Polygon(v3, v2, v1, v0),
	    Polygon(v0, v1, v5, v4),
	    Polygon(v1, v2, v6, v5),
	    Polygon(v2, v3, v7, v6),
	    Polygon(v3, v0, v4, v7),
	    Polygon(v4, v5, v6, v7)}
  {
    center[0] = (v0[0] + v1[0] + v2[0] + v3[0] + v4[0] + v5[0] + v6[0] + v7[0])/8;
    center[1] = (v0[1] + v1[1] + v2[1] + v3[1] + v4[1] + v5[1] + v6[1] + v7[1])/8;
    center[2] = (v0[2] + v1[2] + v2[2] + v3[2] + v4[2] + v5[2] + v6[2] + v7[2])/8;
    volume = calcVolume();
  }

  void apply(const vector_t& translation, double scaling) {
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
	 vertices[4][i] + vertices[5][i] + vertices[6][i] + vertices[7][i])/8;
    volume = calcVolume();
  }

  double calcVolume() const {
    vector_t d(vertices[6] - vertices[0]);
    return
      (util_diff_dot_cross(d, vertices[1], vertices[0], vertices[2], vertices[5]) +
       util_diff_dot_cross(d, vertices[4], vertices[0], vertices[5], vertices[7]) +
       util_diff_dot_cross(d, vertices[3], vertices[0], vertices[7], vertices[2]))/6;
  }

 public:
  vector_t vertices[nrVertices];
  Polygon faces[nrFaces];
  vector_t center;
  double volume;
};

class Sphere {
 public:
  Sphere(double x, double y, double z, double r)
      : center(vector_t(x, y, z))
      , radius(r)
      , volume(double(4.0 / 3.0 * pi) * r * r * r) {}

  double capVolume(double h) const {
    if (h <= double(0))
      return double(0);
    else if (h >= double(2) * radius)
      return volume;
    else
      return double(pi / 3.0) * h * h * (double(3) * radius - h);
  }

 public:
  vector_t center;
  double radius;
  double volume;
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
       const vector_t& min = vector_t(DBL_MAX, DBL_MAX, DBL_MAX),
       const vector_t& max = vector_t(-DBL_MAX, -DBL_MAX, -DBL_MAX))
      : min(min), max(max) {}

  bool intersects(const AABB& aabb) const {
    return
      !(min[0] > aabb.max[0] ||
	min[1] > aabb.max[1] ||
	min[2] > aabb.max[2] ||
	max[0] < aabb.min[0] ||
	max[1] < aabb.min[1] ||
	max[2] < aabb.min[2]);
  }

  void include(int n, const vector_t *p) {
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
  vector_t min, max;
};

// The (convex!) polygon is assumed to be planar, making this a 2D problem.
// Check the projection of the point onto the plane of the polygon for
// containment within the polygon.
bool contains(const Polygon& poly, const vector_t& point) {
  const vector_t proj(
      point - poly.normal.dot(point - poly.center) * poly.normal);
  const int nsize = sizeof poly.vertices / sizeof poly.vertices[0];
  for (size_t n = 0; n < nsize; ++n) {
    const auto& v0 = poly.vertices[n];
    const auto& v1 = poly.vertices[(n + 1) % nsize];
    vector_t base(double(0.5) * (v0 + v1));
    vector_t edge(v1 - v0);

    // Note: Only the sign of the projection is of interest, so this vector
    // does not have to be normalized.
    vector_t dir(edge.cross(poly.normal));

    // Check whether the projection of the point lies inside of the
    // polygon.
    if (dir.dot(proj - base) > double(0)) return false;
  }

  return true;
}

inline bool contains(const Hexahedron& hex, const vector_t& p) {
  for (const auto& f : hex.faces)
    if (f.normal.dot(p - f.center) > double(0)) return false;

  return true;
}

inline bool intersect0(const Sphere& s, const vector_t& center, const vector_t &normal) {
  double proj = normal.dot(s.center - center);
  return proj * proj - s.radius * s.radius < 0;
}

inline bool intersect(const Sphere& s, const Polygon& poly) {
  return intersect0(s, poly.center, poly.normal) &&
	 contains(poly, s.center);
}

inline std::pair<std::array<double, 2>, size_t> lineSphereIntersection(
    const vector_t& origin, const vector_t& direction, const Sphere& s) {
  std::array<double, 2> solutions = {{DBL_MAX, DBL_MAX}};

  vector_t originRel(origin - s.center);
  double a = util_squared_norm(direction[0], direction[1], direction[2]);

  if (a == 0.0) return std::make_pair(solutions, 0);

  double b = double(2) * direction.dot(originRel);
  double c = originRel.squaredNorm() - s.radius * s.radius;

  double discriminant = b * b - double(4) * a * c;
  if (discriminant > double(0)) {
    // Two real roots.
    double q = double(-0.5) * (b + std::copysign(std::sqrt(discriminant), b));

    solutions[0] = q / a;
    solutions[1] = c / q;

    if (solutions[0] > solutions[1]) std::swap(solutions[0], solutions[1]);

    return std::make_pair(solutions, 2);
  } else if (std::abs(discriminant) == double(0)) {
    // Double real root.
    solutions[0] = (double(-0.5) * b) / a;
    solutions[1] = solutions[0];

    return std::make_pair(solutions, 1);
  } else {
    // No real roots.
    return std::make_pair(solutions, 0);
  }
}

// Calculate the volume of a regularized spherical wedge defined by the radius,
// the distance of the intersection point from the center of the sphere and the
// angle.
inline double regularizedWedge(double r, double d, double alpha) {
#ifndef NDEBUG
  // Clamp slight deviations of the angle to valid range.
  if (alpha < double(0) && alpha > -tinyEpsilon) alpha = double(0);

  if (alpha > double(0.5 * pi) && alpha < double(0.5 * pi) + tinyEpsilon)
    alpha = double(0.5 * pi);
#endif

  // Check the parameters for validity (debug version only).
  assert(r > double(0));
  assert(d >= double(0) && d <= r);
  assert(alpha >= double(0) && alpha <= double(0.5 * pi));

  const double sinAlpha = std::sin(alpha);
  const double cosAlpha = std::cos(alpha);

  const double a = d * sinAlpha;
  const double b = std::sqrt(std::abs(r * r - d * d));
  const double c = d * cosAlpha;

  return double(1.0 / 3.0) * a * b * c +
	 a * (double(1.0 / 3.0) * a * a - r * r) * std::atan2(b, c) +
	 double(2.0 / 3.0) * r * r * r * std::atan2(sinAlpha * b, cosAlpha * r);
}

// Wrapper around the above function handling correctly handling the case of
// alpha > pi/2 and negative z.
inline double regularizedWedge(double r, double d, double alpha, double z) {
  if (z >= double(0)) {
    if (alpha > double(0.5 * pi)) {
      double h = r - z;

      return double(pi / 3.0) * h * h * (double(3) * r - h) -
	     regularizedWedge(r, d, pi - alpha);
    } else {
      return regularizedWedge(r, d, alpha);
    }
  } else {
    double vHem = double(2.0 / 3.0 * pi) * r * r * r;

    if (alpha > double(0.5 * pi)) {
      return vHem - regularizedWedge(r, d, pi - alpha);
    } else {
      double h = r + z;
      double vCap = double(pi / 3.0) * h * h * (double(3) * r - h);

      return vHem - (vCap - regularizedWedge(r, d, alpha));
    }
  }
}

inline double generalWedge(
    const Sphere& s, const Plane& p0, const Plane& p1, const vector_t& d) {
  double dist = util_stable_norm(d[0], d[1], d[2]);
  if (dist < tinyEpsilon) {
    // The wedge (almost) touches the center, the volume depends only on
    // the angle.
    double angle = pi - vangle(p0.normal, p1.normal);
    return double(2.0 / 3.0) * s.radius * s.radius * s.radius * angle;
  }

  double s0 = d.dot(p0.normal);
  double s1 = d.dot(p1.normal);

  // Detect degenerated general spherical wedge that can be treated as
  // a regularized spherical wedge.
  if (std::abs(s0) < tinyEpsilon || std::abs(s1) < tinyEpsilon) {
    double angle = pi - vangle(p0.normal, p1.normal);

    return regularizedWedge(
	s.radius, dist, angle, std::abs(s0) > std::abs(s1) ? s0 : s1);
  }

  vector_t dUnit(d * (double(1) / dist));
  if (dist < largeEpsilon)
    dUnit = gramSchmidt(p0.normal.cross(p1.normal), dUnit);

  // Check the planes specify a valid setup (debug version only).
  assert(p0.normal.dot(p1.center - p0.center) <= double(0));
  assert(p1.normal.dot(p0.center - p1.center) <= double(0));

  // Calculate the angles between the vector from the sphere center
  // to the intersection line and the normal vectors of the two planes.
  double alpha0 = vangle(p0.normal, dUnit);
  double alpha1 = vangle(p1.normal, dUnit);

  double dir0 = dUnit.dot((s.center + d) - p0.center);
  double dir1 = dUnit.dot((s.center + d) - p1.center);

  if (s0 >= double(0) && s1 >= double(0)) {
    alpha0 = double(0.5 * pi) - std::copysign(alpha0, dir0);
    alpha1 = double(0.5 * pi) - std::copysign(alpha1, dir1);

    return regularizedWedge(s.radius, dist, alpha0, s0) +
	   regularizedWedge(s.radius, dist, alpha1, s1);
  } else if (s0 < double(0) && s1 < double(0)) {
    alpha0 = double(0.5 * pi) + std::copysign(double(1), dir0) * (alpha0 - pi);

    alpha1 = double(0.5 * pi) + std::copysign(double(1), dir1) * (alpha1 - pi);

    return s.volume - (regularizedWedge(s.radius, dist, alpha0, -s0) +
		       regularizedWedge(s.radius, dist, alpha1, -s1));
  } else {
    alpha0 =
	double(0.5 * pi) - std::copysign(double(1), dir0* s0) *
			       (alpha0 - (s0 < double(0) ? pi : double(0)));

    alpha1 =
	double(0.5 * pi) - std::copysign(double(1), dir1* s1) *
			       (alpha1 - (s1 < double(0) ? pi : double(0)));

    double volume0 = regularizedWedge(s.radius, dist, alpha0, std::abs(s0));

    double volume1 = regularizedWedge(s.radius, dist, alpha1, std::abs(s1));

    return std::max(volume0, volume1) - std::min(volume0, volume1);
  }
}

// Depending on the dimensionality, either the volume or external surface area
// of the general wedge is computed.
double generalWedge(
    const Sphere& sphere, const Hexahedron& element, size_t edge,
    const std::array<std::array<vector_t, 2>, nrEdges>& intersections) {
  const auto& f0 = element.faces[edge_mapping[edge][1][0]];
  const auto& f1 = element.faces[edge_mapping[edge][1][1]];

  vector_t edgeCenter(
      double(0.5) *
      ((intersections[edge][0] + element.vertices[edge_mapping[edge][0][0]]) +
       (intersections[edge][1] + element.vertices[edge_mapping[edge][0][1]])));

  Plane p0(f0.center, f0.normal);
  Plane p1(f1.center, f1.normal);

  return generalWedge(sphere, p0, p1, edgeCenter - sphere.center);
}

double overlap(const Sphere& sOrig, const Hexahedron& elementOrig) {
  // Construct AABBs and perform a coarse overlap detection.
  AABB sAABB(
      sOrig.center - vector_t::Constant(sOrig.radius),
      sOrig.center + vector_t::Constant(sOrig.radius));

  AABB eAABB;
  eAABB.include(sizeof elementOrig.vertices / sizeof *elementOrig.vertices,
		elementOrig.vertices);

  if (!sAABB.intersects(eAABB)) return double(0);

  // Use scaled and shifted versions of the sphere and the element.
  vector_t translation = -sOrig.center;
  double scaling = 1/sOrig.radius;

  Sphere s(double(0), double(0), double(0), double(1));

  Hexahedron element(elementOrig);
  element.apply(translation, scaling);

  size_t vOverlap = 0;
  // Check whether the vertices lie on or outside of the sphere.
  for (const auto& vertex : element.vertices)
    if ((s.center - vertex).squaredNorm() <= s.radius * s.radius) ++vOverlap;

  // Check for trivial case: All vertices inside of the sphere, resulting in
  // a full overlap.
  if (vOverlap == nrVertices) return elementOrig.volume;

  // Sets of overlapping primitives.
  std::bitset<nrVertices> vMarked;
  std::bitset<nrEdges> eMarked;
  std::bitset<nrFaces> fMarked;

  // Initial value: Volume of the full sphere.
  double result = s.volume;

  // The intersection points between the single edges and the sphere, this
  // is needed later on.
  std::array<std::array<vector_t, 2>, nrEdges> eIntersections;

  // Process all edges of the element.
  for (size_t n = 0; n < nrEdges; ++n) {
    vector_t start(element.vertices[edge_mapping[n][0][0]]);
    vector_t direction(element.vertices[edge_mapping[n][0][1]] - start);

    auto solutions = lineSphereIntersection(start, direction, s);

    // No intersection between the edge and the sphere, where intersection
    // points close to the surface of the sphere are ignored.
    // Or:
    // The sphere cuts the edge twice, no vertex is inside of the
    // sphere, but the case of the edge only touching the sphere has to
    // be avoided.
    if (!solutions.second ||
	(solutions.first[0] >= double(1) - mediumEpsilon) ||
	solutions.first[1] <= mediumEpsilon ||
	(solutions.first[0] > double(0) && solutions.first[1] < double(1) &&
	 (solutions.first[1] - solutions.first[0] < largeEpsilon))) {
      continue;
    } else {
      vMarked[edge_mapping[n][0][0]] = solutions.first[0] < double(0);

      vMarked[edge_mapping[n][0][1]] = solutions.first[1] > double(1);
    }

    // Store the two intersection points of the edge with the sphere for
    // later usage.
    eIntersections[n][0] = solutions.first[0] * direction +
			   (start - element.vertices[edge_mapping[n][0][0]]);
    eIntersections[n][1] = solutions.first[1] * direction +
			   (start - element.vertices[edge_mapping[n][0][1]]);

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
    if (intersect(s, element.faces[n])) fMarked[n] = true;

  // Trivial case: The center of the sphere overlaps the element, but the
  // sphere does not intersect any of the faces of the element, meaning the
  // sphere is completely contained within the element.
  if (!fMarked.count() && contains(element, s.center)) return sOrig.volume;

  // Spurious intersection: The initial intersection test was positive, but
  // the detailed checks revealed no overlap.
  if (!vMarked.count() && !eMarked.count() && !fMarked.count())
    return double(0);

  // Iterate over all the marked faces and subtract the volume of the cap cut
  // off by the plane.
  for (size_t n = 0; n < nrFaces; ++n) {
    if (!fMarked[n]) continue;

    const auto& f = element.faces[n];
    double dist = f.normal.dot(s.center - f.center);
    double vCap = s.capVolume(s.radius + dist);

    result -= vCap;
  }

  // Handle the edges and add back the volume subtracted twice above in the
  // processing of the faces.
  for (size_t n = 0; n < nrEdges; ++n) {
    if (!eMarked[n]) continue;

    double edgeCorrection = generalWedge(s, element, n, eIntersections);

    result += edgeCorrection;
  }

  // Handle the vertices and subtract the volume added twice above in the
  // processing of the edges.
  for (size_t n = 0; n < nrVertices; ++n) {
    if (!vMarked[n]) continue;

    // Collect the points where the three edges intersecting at this
    // vertex intersect the sphere.
    // Both the relative and the absolute positions are required.
    vector_t intersectionPointsRelative[3];
    vector_t intersectionPoints[3];
    for (size_t e = 0; e < 3; ++e) {
      auto edgeIdx = vertex_mapping[n][0][e];
      intersectionPointsRelative[e] =
	  eIntersections[edgeIdx][vertex_mapping[n][1][e]];

      intersectionPoints[e] =
	  intersectionPointsRelative[e] + element.vertices[n];
    }

    // This triangle is constructed by hand to have more freedom of how
    // the normal vector is calculated.

    vector_t center = double(1.0 / 3.0) * (intersectionPoints[0] +
					   intersectionPoints[1] +
					   intersectionPoints[2]);

    // Calculate the normal of the triangle defined by the intersection
    // points in relative coordinates to improve accuracy.
    // Also use double the normal precision to calculate this normal.
    vector_t normal = triangleNormal(
	intersectionPointsRelative[0], intersectionPointsRelative[1],
	intersectionPointsRelative[2]);

    std::array<std::pair<size_t, double>, 3> distances;
    for (size_t i = 0; i < 3; ++i)
      distances[i] =
	std::make_pair(i, util_squared_norm(intersectionPointsRelative[i][0],
					    intersectionPointsRelative[i][1],
					    intersectionPointsRelative[i][2]));
    std::sort(
	distances.begin(), distances.end(),
	[](const std::pair<size_t, double>& a,
	   const std::pair<size_t, double>& b) -> bool {
	  return a.second < b.second;
	});

    if (distances[1].second < distances[2].second * largeEpsilon) {
      // Use the general spherical wedge defined by the edge with the
      // non-degenerated intersection point and the normals of the
      // two faces forming it.
      double correction = generalWedge(
	  s, element, vertex_mapping[n][0][distances[2].first], eIntersections);

      result -= correction;

      continue;
    }

    double tipTetVolume =
	double(1.0 / 6.0) *
	std::abs(-intersectionPointsRelative[2].dot(
	    (intersectionPointsRelative[0] - intersectionPointsRelative[2])
		.cross(
		    intersectionPointsRelative[1] -
		    intersectionPointsRelative[2])));

    // Make sure the normal points in the right direction i.e. away from
    // the center of the element.
    if (normal.dot(element.center - center) > double(0)) {
      normal = -normal;
    }

    Plane plane(center, normal);

    double dist = normal.dot(s.center - center);
    double capVolume = s.capVolume(s.radius + dist);

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
    double segmentVolume = 0;

    for (size_t e = 0; e < 3; ++e) {
      const auto& f = element.faces[vertex_mapping[n][2][e]];
      uint32_t e0 = face_mapping[e][0];
      uint32_t e1 = face_mapping[e][1];

      vector_t center(
	  double(0.5) * (intersectionPoints[e0] + intersectionPoints[e1]));

      double wedgeVolume =
	  generalWedge(s, plane, Plane(f.center, -f.normal), center - s.center);

      segmentVolume += wedgeVolume;
    }

    // Calculate the volume of the cone and clamp it to zero.
    double coneVolume =
	std::max(tipTetVolume + capVolume - segmentVolume, double(0));

    // Sanity check: detect negative cone volume.
    assert(coneVolume > -std::sqrt(tinyEpsilon));

    result -= coneVolume;

    // Sanity check: detect negative intermediate result.
    assert(result > -std::sqrt(tinyEpsilon));
  }

  // In case of different sized objects the error can become quite large,
  // so a relative limit is used.
  double maxOverlap = std::min(s.volume, element.volume);
  const double limit(std::sqrt(DBL_EPSILON) * maxOverlap);

  // Clamp tiny negative volumes to zero.
  if (result < double(0) && result > -limit) return double(0);

  // Clamp results slightly too large.
  if (result > maxOverlap && result - maxOverlap < limit)
    return std::min(sOrig.volume, elementOrig.volume);

  // Perform a sanity check on the final result (debug version only).
  assert(result >= double(0) && result <= maxOverlap);

  // Scale the overlap volume back for the original objects.
  result = (result / s.volume) * sOrig.volume;

  return result;
}
