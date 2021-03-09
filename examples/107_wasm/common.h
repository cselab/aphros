void ErrorHandler(int code, const char* str) {
  fputs(str, stderr);
  fputs("\n", stderr);
}

auto ToMulti(const std::vector<Scal>& v) {
  Multi<Scal> w(v.size());
  w.data() = v;
  return w;
};

void CopyToCanvas(uint32_t* buf, int w, int h) {
  EM_ASM_(
      {
        let data = Module.HEAPU8.slice($0, $0 + $1 * $2 * 4);
        let tctx = g_tmp_canvas.getContext("2d");
        let image = tctx.getImageData(0, 0, $1, $2);
        image.data.set(data);
        tctx.putImageData(image, 0, 0);
      },
      buf, w, h);
}

struct Canvas {
  Canvas(MIdx size_) : size(size_), buf(size.prod()) {}
  MIdx size;
  Scal scale = 1;
  int xmargin = 0;
  std::vector<uint32_t> buf;
};

inline Scal Clamp(Scal f, Scal min, Scal max) {
  return f < min ? min : f > max ? max : f;
}

template <class T>
inline T Clamp(T v) {
  return v.max(T(0)).min(T(1));
}

inline Scal Clamp(Scal f) {
  return Clamp(f, 0, 1);
}

void GetCircle(FieldCell<Scal>& fcu, Vect c, Scal r, const M& m) {
  Vars par;
  par.String.Set("init_vf", "circlels");
  par.Vect.Set("circle_c", c);
  par.Double.Set("circle_r", r);
  par.Int.Set("dim", 2);
  auto func = CreateInitU<M>(par, false);
  func(fcu, m);
}

template <class M>
static MIdx GetCanvasCoords(Vect x, const Canvas& canvas, const M& m) {
  auto scaledsize = canvas.size * canvas.scale;
  MIdx w = MIdx(x * Vect(scaledsize) / m.GetGlobalLength());
  w = w.max(MIdx(0)).min(scaledsize - MIdx(1));
  w[1] = scaledsize[1] - w[1] - 1;
  w[0] += canvas.xmargin;
  return w;
}

template <class MEB>
FieldCell<typename MEB::Scal> GetVortScal(
    const FieldCell<typename MEB::Vect>& fcvel,
    const MapEmbed<BCond<typename MEB::Vect>>& me_vel, MEB& eb) {
  using M = typename MEB::M;
  constexpr auto dim = M::dim;
  auto& m = eb.GetMesh();
  using Scal = typename MEB::Scal;
  using Vect = typename MEB::Vect;
  using UEB = UEmbed<M>;

  std::array<FieldCell<Vect>, dim> grad;
  for (size_t d = 0; d < dim; ++d) {
    grad[d].Reinit(m, Vect(0));
    const auto mebc = GetScalarCond(me_vel, d, m);
    const FieldCell<Scal> fcu = GetComponent(fcvel, d);
    const FieldFace<Scal> ffg = UEB::Gradient(fcu, mebc, m);
    grad[d] = UEB::AverageGradient(ffg, m);
  }

  FieldCell<Scal> res(m, 0);
  for (auto c : m.Cells()) {
    res[c] = grad[1][c][0] - grad[0][c][1];
  }
  return res;
}
