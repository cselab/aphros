// Created by Petr Karnakov on 07.03.2021
// Copyright 2021 ETH Zurich

#include <fstream>

#include "geom/mesh.h"

namespace util {

template <class M_>
struct Visual {
  using M = M_;
  using Scal = typename M::Scal;
  using MIdx = typename M::MIdx;
  using MIdx2 = generic::MIdx<2>;
  using Vect = typename M::Vect;
  using Vect3 = generic::Vect<Scal, 3>;
  using Float3 = generic::Vect<float, 3>;
  using Byte3 = generic::Vect<unsigned char, 3>;

  using Pixel = uint32_t;

  static Scal Clamp(Scal f, Scal min, Scal max) {
    return f < min ? min : f > max ? max : f;
  }

  template <class T>
  static T Clamp(T v) {
    return v.max(T(0)).min(T(1));
  }

  static Scal Clamp(Scal f) {
    return Clamp(f, 0, 1);
  }

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
        : CanvasView(canvas, MIdx(0), canvas.size) {}
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
    Scal min = 0;
    Scal max = 1;
    Float3 color_min{1};
    Float3 color_max{0};
    Scal opacity_min = 1;
    Scal opacity_max = 1;
    Float3 operator()(Float3 src, Scal u) const {
      u = Clamp((u - min) / (max - min));
      const Scal a = opacity_min * (1 - u) + opacity_max * u;
      return src * (1 - a) + (color_min * (1 - u) + color_max * u) * a;
    }
  };

  static void AppendField(
      FieldCell<Float3>& fc_color, const FieldCell<Scal>& fcu,
      const Colormap& cmap, const M& m) {
    for (auto c : m.Cells()) {
      fc_color[c] = cmap(fc_color[c], fcu[c]);
    }
  }

  static void RenderColorField(
      const FieldCell<Float3>& fc_color, CanvasView& view, const M& m) {
    const auto msize = m.GetGlobalSize();
    for (auto c : m.CellsM()) {
      const MIdx w(c);
      const MIdx start = w * view.size / msize;
      const MIdx end = (w + MIdx(1)) * view.size / msize;
      const Byte3 q(Clamp(fc_color[c]) * 255);
      const Pixel v = 0xff000000 | (q[0] << 0) | (q[1] << 8) | (q[2] << 16);
      for (int y = start[1]; y < end[1]; y++) {
        for (int x = start[0]; x < end[0]; x++) {
          view(x, y) = v;
        }
      }
    }
  }

  static void WritePpm(std::ostream& out, const CanvasView& view) {
    out << "P3\n";
    out << view.size[0] << ' ' << view.size[1] << '\n';
    out << 255 << '\n';
    for (int y = view.size[1]; y > 0;) {
      --y;
      for (int x = 0; x < view.size[0]; ++x) {
        out << (view(x, y) & 0xFF) << ' ';
        out << ((view(x, y) >> 8) & 0xFF) << ' ';
        out << ((view(x, y) >> 16) & 0xFF) << ' ';
      }
      out << '\n';
    }
  }

  static void WritePpm(std::string path, const CanvasView& view) {
    std::ofstream f(path);
    fassert(f.good(), "Can't open file '" + path + "' for writing");
    WritePpm(f, view);
  }

  static Canvas ReadPpm(std::istream& in) {
    std::string line;

    std::getline(in, line);
    fassert_equal(line, "P3");

    MIdx2 size;
    in >> size[0] >> size[1];

    int max_color;
    in >> max_color;
    fassert(max_color == 255);

    Canvas canvas(size);
    CanvasView view(canvas);

    for (int y = canvas.size[1]; y > 0;) {
      --y;
      for (int x = 0; x < canvas.size[0]; ++x) {
        generic::MIdx<3> t;
        in >> t;
        Byte3 q(t);
        view(x, y) = 0xFF000000 | (q[0] << 0) | (q[1] << 8) | (q[2] << 16);
      }
    }
    return canvas;
  }

  static Canvas ReadPpm(std::string path) {
    std::ifstream f(path);
    fassert(f.good(), "Can't open file '" + path + "' for reading");
    return ReadPpm(f);
  }
};

} // namespace util
