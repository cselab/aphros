// specialization for small alphabets -----------------------------------------

template <typename T, class M>
class PCdecoder<T, M, false> {
public:
  PCdecoder(RCdecoder* rd, RCmodel*const* rm) : rd(rd), rm(rm) {}
  ~PCdecoder() {}
  T decode(T pred, unsigned context = 0);
  static const unsigned symbols = 2 * (1 << M::bits) - 1;
private:
  static const unsigned bias = (1 << M::bits) - 1;
  M                     map;            // maps T to some unsigned int type
  RCdecoder*const       rd;             // entropy decoder
  RCmodel*const*        rm;             // probability modeler(s)
};

// decode narrow range type
template <typename T, class M>
T PCdecoder<T, M, false>::decode(T pred, unsigned context)
{
  // map type T to unsigned integer type
  typedef typename M::RANGE U;
  U p = map.forward(pred);
  // entropy decode d = r - p
  U r = p + rd->decode(rm[context]) - bias;
  return map.inverse(r);
}

// specialization for large alphabets -----------------------------------------

template <typename T, class M>
class PCdecoder<T, M, true> {
public:
  PCdecoder(RCdecoder* rd, RCmodel*const* rm) : rd(rd), rm(rm) {}
  ~PCdecoder() {}
  T decode(T pred, unsigned context = 0);
  static const unsigned symbols = 2 * M::bits + 1;
private:
  static const unsigned bias = M::bits;
  M                     map;            // maps T to some unsigned int type
  RCdecoder*const       rd;             // entropy decoder
  RCmodel*const*        rm;             // probability modeler(s)
};

// decode wide range type
template <typename T, class M>
T PCdecoder<T, M, true>::decode(T pred, unsigned context)
{
  typedef typename M::RANGE U;
  unsigned s = rd->decode(rm[context]);
  if (s > bias) {      // underprediction
    unsigned k = s - bias - 1;
    U d = (U(1) << k) + rd->template decode<U>(k);
    U p = map.forward(pred);
    U r = p + d;
    return map.inverse(r);
  }
  else if (s < bias) { // overprediction
    unsigned k = bias - 1 - s;
    U d = (U(1) << k) + rd->template decode<U>(k);
    U p = map.forward(pred);
    U r = p - d;
    return map.inverse(r);
  }
  else                 // perfect prediction
    return map.identity(pred);
}
