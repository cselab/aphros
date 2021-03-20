// Created by Petr Karnakov on 07.03.2021
// Copyright 2021 ETH Zurich

#pragma once

#include <iosfwd>
#include <string>
#include <vector>

#include "parse/vars.h"

namespace util {

using Float3 = generic::Vect<float, 3>;
using Byte3 = generic::Vect<unsigned char, 3>;
using MIdx2 = generic::MIdx<2>;
using Pixel = uint32_t;

struct Canvas {
  Canvas(MIdx2 size_) : size(size_), buf(size.prod()) {}
  const MIdx2 size;
  std::vector<Pixel> buf;
};

struct CanvasView {
  CanvasView(MIdx2 size_, MIdx2 start_, MIdx2 end_, Pixel* buf_)
      : size(size_), start(start_), end(end_), buf(buf_) {}
  CanvasView(Canvas& canvas, MIdx2 start_, MIdx2 end_)
      : size(canvas.size), start(start_), end(end_), buf(canvas.buf.data()) {}
  explicit CanvasView(Canvas& canvas)
      : CanvasView(canvas, MIdx2(0), canvas.size) {}
  MIdx2 size; // full canvas size
  MIdx2 start; // start of selection
  MIdx2 end; // end of selection
  Pixel* buf;
  auto& operator()(int x, int y) {
    return buf[size[0] * (size[1] - y - 1) + x];
  }
  const auto& operator()(int x, int y) const {
    return buf[size[0] * (size[1] - y - 1) + x];
  }
};

struct Colormap {
  using Scal = double;
  template <class T>
  static T Piecewise(
      Scal x, const std::vector<Scal>& xx, const std::vector<T>& yy) {
    if (x <= xx.front()) {
      return yy.front();
    }
    if (x >= xx.back()) {
      return yy.back();
    }
    // invariant: xx[left] <= x < xx[right]
    size_t left = 0;
    size_t right = xx.size() - 1;
    while (left + 1 < right) {
      size_t mid = (left + right) / 2;
      if (xx[mid] <= x) {
        left = mid;
      } else { // x < xx[mid]
        right = mid;
      }
    }
    return xx[left] < xx[right]
               ? yy[left] + (yy[right] - yy[left]) *
                                ((x - xx[left]) / (xx[right] - xx[left]))
               : yy[left];
  }
  Float3 operator()(Float3 src, Scal u) const {
    const auto color = Piecewise(u, values, colors);
    const auto opacity = Piecewise(u, values, opacities);
    return src * (1 - opacity) + color * opacity;
  }
  void Check() const {
    fassert(values.size() >= 1);
    fassert(values.size() <= colors.size());
    fassert(values.size() <= opacities.size());
  }

  std::vector<Scal> values; // values of field in ascending order
  std::vector<Float3> colors; // colors for each field value
  std::vector<Scal> opacities; // opacities for each field value
};

struct Entry {
  std::string field;
  Vars var;
};

std::vector<Entry> ParseEntries(std::istream& in);
std::vector<Entry> ParseEntries(std::string path);

Colormap GetColormap(const Vars& var);

void WritePpm(std::ostream& out, const CanvasView& view, bool binary = true);
void WritePpm(std::string path, const CanvasView& view, bool binary = true);
Canvas ReadPpm(std::istream& in);
Canvas ReadPpm(std::string path);

template <class M_>
struct Visual {
  using M = M_;
  using Scal = typename M::Scal;
  using MIdx = typename M::MIdx;
  using Vect = typename M::Vect;
  using Vect3 = generic::Vect<Scal, 3>;

  static void RenderToField(
      FieldCell<Float3>& fc_color, const FieldCell<Scal>& fcu,
      const Colormap& cmap, const M& m);

  static void RenderEntriesToField(
      FieldCell<Float3>& fc_color, const std::vector<Entry>& entries,
      const std::function<FieldCell<Scal>(std::string)>& get_field, const M& m);

  static void RenderToCanvasNearest(
      CanvasView& view, const FieldCell<Float3>& fc_color, const M& m);

  static void RenderToCanvasBilinear(
      CanvasView& view, const FieldCell<Float3>& fc_color, const M& m);

  struct Imp;
};

} // namespace util
