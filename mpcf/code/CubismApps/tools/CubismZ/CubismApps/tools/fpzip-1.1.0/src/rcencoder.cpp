#include <cstddef>
#include "rcencoder.h"

// finalize encoder
void RCencoder::finish()
{
  put(4);
  flush();
}

// encode a bit s
void RCencoder::encode(bool s)
{
  range >>= 1;
  if (s)
    low += range;
  normalize();
}

// encode a symbol s using probability modeling
void RCencoder::encode(unsigned s, RCmodel* rm)
{
  unsigned l, r;
  rm->encode(s, l, r);
  rm->normalize(range);
  low += range * l;
  range *= r;
  normalize();
}

// encode a number s : 0 <= s < 2^n <= 2^16
void RCencoder::encode_shift(unsigned s, unsigned n)
{
  range >>= n;
  low += range * s;
  normalize();
}

// encode a number s : 0 <= s < n <= 2^16
void RCencoder::encode_ratio(unsigned s, unsigned n)
{
  range /= n;
  low += range * s;
  normalize();
}

// normalize the range and output data
void RCencoder::normalize()
{
  while (!((low ^ (low + range)) >> 24)) {
    // top 8 bits are fixed; output them
    put(1);
    range <<= 8;
  }
  if (!(range >> 16)) {
    // top 8 bits are not fixed but range is small;
    // fudge range to avoid carry and output 16 bits
    put(2);
    range = -low;
  }
}
