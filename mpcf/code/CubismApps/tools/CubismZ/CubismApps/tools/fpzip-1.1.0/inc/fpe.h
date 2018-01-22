#ifndef FPE_H
#define FPE_H

// Approximate floating-point arithmetic for normalized inputs and outputs.
// NOTE: This code has been designed for speed and simplicity, and does not
// fully implement IEEE floating-point arithmetic.  In particular, it does
// not properly handle:
//
//   (1) Denormalized numbers, infinities, or NaNs.
//   (2) Overflow or underflow outside the range of normalized numbers.
//   (3) Proper IEEE rounding.
//
// For normalized input and output, the error of a single operation should
// be no more than one ulp.

template <typename T>
struct FPEtraits;

template <>
struct FPEtraits<float> {
  typedef unsigned U;
  static const unsigned mbits = 23;
  static const unsigned ebits = 8;
};

template <>
struct FPEtraits<double> {
  typedef unsigned long long U;
  static const unsigned mbits = 52;
  static const unsigned ebits = 11;
};

template <
  typename T,                            // floating type to implement
  typename U = typename FPEtraits<T>::U, // corresponding integer type
  unsigned mbits = FPEtraits<T>::mbits,  // number of bits in significand
  unsigned ebits = FPEtraits<T>::ebits   // number of bits in exponent
>
class FPE {
public:
  FPE(U u = 0) : u(u) {}
  FPE(T t) : u((U&)t) {}
  FPE(int i) : u(i) {}
  operator T() const { return (T&)u; }
  bool equalsign(FPE x) const { return !((u ^ x.u) & s); }
  FPE operator -() const { return FPE(u ^ s); }
  FPE operator +(FPE x) const;
  FPE operator -(FPE x) const;
private:
  FPE add(U g) const;                 // effective addition
  FPE sub(U g) const;                 // effective subtraction
  static const unsigned m = mbits;    // number of bits in significand
  static const U e = U(1) << mbits;   // exponent least significant bit mask
  static const U s = e << ebits;      // sign bit mask
  U u;                                // binary representation of float
};

#include "fpe.inl"

#endif
