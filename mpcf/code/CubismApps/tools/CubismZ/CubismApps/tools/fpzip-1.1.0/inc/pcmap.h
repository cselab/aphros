#ifndef PC_MAP_H
#define PC_MAP_H

#include <climits>
#if !defined WITH_REINTERPRET_CAST && !defined WITH_UNION
#include <cstring>
#endif

#define bitsizeof(t) ((unsigned)(CHAR_BIT * sizeof(t)))

template <typename T, unsigned width = bitsizeof(T), typename U = void>
struct PCmap;

// specialized for integer-to-integer map
template <typename T, unsigned width>
struct PCmap<T, width, void> {
  typedef T DOMAIN;
  typedef T RANGE;
  static const unsigned bits = width;                    // RANGE bits
  static const unsigned shift = bitsizeof(RANGE) - bits; // DOMAIN\RANGE bits
  RANGE forward(DOMAIN d) const { return d >> shift; }
  DOMAIN inverse(RANGE r) const { return r << shift; }
  DOMAIN identity(DOMAIN d) const { return inverse(forward(d)); }
};

// specialized for float type
template <unsigned width>
struct PCmap<float, width, void> {
  typedef float    DOMAIN;
  typedef unsigned RANGE;
  union UNION {
    UNION(DOMAIN d) : d(d) {}
    UNION(RANGE r) : r(r) {}
    DOMAIN d;
    RANGE r;
  };
  static const unsigned bits = width;                    // RANGE bits
  static const unsigned shift = bitsizeof(RANGE) - bits; // DOMAIN\RANGE bits
  RANGE fcast(DOMAIN d) const;
  DOMAIN icast(RANGE r) const;
  RANGE forward(DOMAIN d) const;
  DOMAIN inverse(RANGE r) const;
  DOMAIN identity(DOMAIN d) const;
};

// specialized for double type
template <unsigned width>
struct PCmap<double, width, void> {
  typedef double             DOMAIN;
  typedef unsigned long long RANGE;
  union UNION {
    UNION(DOMAIN d) : d(d) {}
    UNION(RANGE r) : r(r) {}
    DOMAIN d;
    RANGE r;
  };
  static const unsigned bits = width;                    // RANGE bits
  static const unsigned shift = bitsizeof(RANGE) - bits; // DOMAIN\RANGE bits
  RANGE fcast(DOMAIN d) const;
  DOMAIN icast(RANGE r) const;
  RANGE forward(DOMAIN d) const;
  DOMAIN inverse(RANGE r) const;
  DOMAIN identity(DOMAIN d) const;
};

#include "pcmap.inl"

#endif
