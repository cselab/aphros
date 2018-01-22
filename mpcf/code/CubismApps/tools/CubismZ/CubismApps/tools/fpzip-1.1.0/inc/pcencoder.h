#ifndef PC_ENCODER_H
#define PC_ENCODER_H

#include "pccodec.h"
#include "pcmap.h"
#include "rcencoder.h"
#include "rcmodel.h"

template <typename T, class M = PCmap<T>, bool wide = (M::bits > PC_BIT_MAX)>
class PCencoder {
public:
  PCencoder(RCencoder* re, RCmodel*const* rm);

  // encode a value with prediction and optional context
  T encode(T real, T pred, unsigned context = 0);

  // number of symbols (needed by probability modeler)
  static const unsigned symbols;
};

#include "pcencoder.inl"

#endif
