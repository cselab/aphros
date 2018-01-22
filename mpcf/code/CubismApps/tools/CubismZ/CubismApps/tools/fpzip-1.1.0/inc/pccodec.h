#ifndef PC_CODEC_H
#define PC_CODEC_H

#ifndef PC_BIT_MAX
  #define PC_BIT_MAX 8 // maximum bit width of directly encodable integers
#endif

namespace PC {
  template <typename U>
  unsigned bsr(U x);
  #include "pccodec.inl"
}

#endif
