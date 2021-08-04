// This implementation of double_prec is based on:
// T.J. Dekker, A floating-point technique for extending the available
// precision, http://dx.doi.org/10.1007/BF01397083

template <typename T>
struct double_prec_constant;

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
    const double_prec<T>& rhs) ENFORCE_EXACT_FPMATH_ATTR;

template <typename T>
class double_prec {
 private:
  static const uint32_t c = double_prec_constant<T>::value;

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
