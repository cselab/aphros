// get frequencies for a symbol s
inline void RCqsmodel::encode(unsigned s, unsigned& l, unsigned& r)
{
  l = cumf[s];
  r = cumf[s + 1] - l;
  update(s);
}

inline void RCqsmodel::normalize(unsigned& r)
{
  r >>= bits;
}
