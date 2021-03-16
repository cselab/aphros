// Created by Petr Karnakov on 07.03.2021
// Copyright 2021 ETH Zurich

#pragma once

#include <fstream>
#include <functional>
#include <string>

#include "geom/mesh.h"
#include "parse/codeblocks.h"
#include "parse/parser.h"
#include "parse/vars.h"
#include "util/logger.h"

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
    std::vector<Scal> values; // values of field in ascending order
    std::vector<Float3> colors; // colors for each field value
    std::vector<Scal> opacities; // opacities for each field value
    template <class T>
    static auto Piecewise(
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
  };

  static void RenderToField(
      FieldCell<Float3>& fc_color, const FieldCell<Scal>& fcu,
      const Colormap& cmap, const M& m) {
    cmap.Check();
    for (auto c : m.SuCells()) {
      fc_color[c] = cmap(fc_color[c], fcu[c]);
    }
  }

  struct Entry {
    std::string field;
    Vars var;
  };
  static std::vector<Entry> ParseEntries(std::istream& in) {
    std::vector<Entry> entries;
    auto blocks = ParseCodeBlocks(in);
    for (auto& b : blocks) {
      Entry entry;
      entry.field = b.name;
      std::stringstream ss(b.content);
      Parser(entry.var).ParseStream(ss);
      entries.push_back(entry);
    }
    return entries;
  }
  static auto ParseEntries(std::string path) {
    std::ifstream f(path);
    fassert(f.good(), "Can't open file '" + path + "' for reading");
    return ParseEntries(f);
  }
  static Colormap GetColormap(const Vars& var) {
    Colormap cmap;
    cmap.values = var.Vect["values"];
    cmap.opacities = var.Vect("opacities", {});
    const auto& colors = var.Vect("colors", {});
    cmap.colors.resize(cmap.values.size());
    for (size_t i = 0; i < std::min(cmap.values.size(), colors.size() / 3);
         ++i) {
      cmap.colors[i][0] = colors[3 * i + 0];
      cmap.colors[i][1] = colors[3 * i + 1];
      cmap.colors[i][2] = colors[3 * i + 2];
    }
    // fill the tail with default values
    for (size_t i = cmap.opacities.size(); i < cmap.values.size(); ++i) {
      cmap.opacities.push_back(1);
    }
    for (size_t i = cmap.colors.size(); i < cmap.values.size(); ++i) {
      cmap.colors.push_back(Float3(0, 0, 0));
    }
    return cmap;
  }
  static void RenderEntriesToField(
      FieldCell<Float3>& fc_color, const std::vector<Entry>& entries,
      const std::function<FieldCell<Scal>(std::string)>& get_field,
      const M& m) {
    for (const auto& entry : entries) {
      auto cmap = GetColormap(entry.var);
      RenderToField(fc_color, get_field(entry.field), cmap, m);
    }
  }

  static void RenderToCanvasNearest(
      CanvasView& view, const FieldCell<Float3>& fc_color, const M& m) {
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

  // Evaluates bilinear interpolant on points (0,0), (1,0), (0,1) and (1,1).
  // x,y: target point
  // u,ux,uy,uyx:  values of function for (x,y) = (0,0), (1,0), (0,1), (1,1)
  // FIXME: this is a copy from approx_eb.ipp
  template <class T, class Scal>
  static T Bilinear(Scal x, Scal y, T u, T ux, T uy, T uyx) {
    //                      //
    //   y                  //
    //   |                  //
    //   |*uy    *uyx       //
    //   |                  //
    //   |                  //
    //   |*u     *ux        //
    //   |-------------x    //
    //                      //
    const auto v = u * (1 - x) + ux * x;
    const auto vy = uy * (1 - x) + uyx * x;
    return v * (1 - y) + vy * y;
  }

  static void RenderToCanvasBilinear(
      CanvasView& view, const FieldCell<Float3>& fc_color, const M& m) {
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

    for (auto c : m.SuCellsM()) {
      const MIdx w(c);
      const MIdx bs = view.size / msize;
      const MIdx start = (w * view.size / msize + bs / 2).max(MIdx(0));
      const MIdx end =
          ((w + MIdx(1)) * view.size / msize + bs / 2).min(view.size);
      const auto dx = m.direction(0);
      const auto dy = m.direction(1);
      const auto q = fc_color[c];
      const auto qx = fc_color[c + dx];
      const auto qy = fc_color[c + dy];
      const auto qyx = fc_color[c + dx + dy];
      for (int y = start[1]; y < end[1]; y++) {
        const Scal fy = Scal(y - start[1]) / (end[1] - start[1]);
        for (int x = start[0]; x < end[0]; x++) {
          const Scal fx = Scal(x - start[0]) / (end[0] - start[0]);
          auto qb = Bilinear(std::abs(fx), std::abs(fy), q, qx, qy, qyx);
          const Byte3 mq(qb * 255);
          Pixel v = 0xff000000 | (mq[0] << 0) | (mq[1] << 8) | (mq[2] << 16);
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
