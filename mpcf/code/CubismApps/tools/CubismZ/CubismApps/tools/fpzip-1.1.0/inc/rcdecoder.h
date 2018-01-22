#ifndef RC_DECODER_H
#define RC_DECODER_H

#include "rcmodel.h"

class RCdecoder {
public:
  RCdecoder() : error(false), low(0), range(-1u), code(0) {}
  virtual ~RCdecoder() {}

  // initialize decoding
  void init();

  // decode a bit
  bool decode();

  // decode a number s : 0 <= s < 2^n
  template <typename UINT>
  UINT decode(unsigned n);

  // decode a number s : l <= s < h
  template <typename UINT>
  UINT decode(UINT l, UINT h);

  // decode a symbol using probability modeling
  unsigned decode(RCmodel* rm);

  // virtual function for reading byte stream
  virtual unsigned getbyte() = 0;

  // number of bytes read
  virtual size_t bytes() const = 0;

  bool error;

private:
  unsigned decode_shift(unsigned n);
  unsigned decode_ratio(unsigned n);
  void get(unsigned n);
  void normalize();

  unsigned low;   // low end of interval
  unsigned range; // length of interval
  unsigned code;  // incoming data
};

#include "rcdecoder.inl"

#endif
