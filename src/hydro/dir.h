#pragma once

namespace geom {

template <size_t dim_>
class GDir {
 public:
  using MIdx = GVect<IntIdx, dim_>;
  static constexpr size_t dim = dim_;

  GDir() {}
  explicit GDir(size_t d)
    : d_(d) {}
  char GetLetter() const {
    return "xyz"[d_];
  }
  explicit operator size_t() const {
    return d_;
  }
  explicit operator MIdx() const {
    MIdx r(0);
    ++r[d_];
    return r;
  }
  bool operator==(const GDir& o) const {
    return d_ == o.d_;
  }
  bool operator!=(const GDir& o) const {
    return !((*this) == o);
  }
  bool operator<(const GDir& o) const {
    return d_ < o.d_;
  }
  static const GDir i;
  static const GDir j;
  static const GDir k;

 private:
  size_t d_;
};

template <size_t dim>
const GDir<dim> GDir<dim>::i(0);
template <size_t dim>
const GDir<dim> GDir<dim>::j(1);
template <size_t dim>
const GDir<dim> GDir<dim>::k(2);

} // namespace geom
