template <unsigned width>
unsigned
PCmap<float, width, void>::fcast(float d) const
{
#ifdef WITH_REINTERPRET_CAST
  return reinterpret_cast<const RANGE&>(d);
#elif defined WITH_UNION
  UNION shared(d);
  return shared.r;
#else
  RANGE r;
  memcpy(&r, &d, sizeof(r));
  return r;
#endif
}

template <unsigned width>
float
PCmap<float, width, void>::icast(unsigned r) const
{
#ifdef WITH_REINTERPRET_CAST
  return reinterpret_cast<const DOMAIN&>(r);
#elif defined WITH_UNION
  UNION shared(r);
  return shared.d;
#else
  DOMAIN d;
  memcpy(&d, &r, sizeof(d));
  return d;
#endif
}

template <unsigned width>
unsigned
PCmap<float, width, void>::forward(float d) const
{
  RANGE r = fcast(d);
  r = ~r;
  r >>= shift;
  r ^= -(r >> (bits - 1)) >> (shift + 1);
  return r;
}

template <unsigned width>
float
PCmap<float, width, void>::inverse(unsigned r) const
{
  r ^= -(r >> (bits - 1)) >> (shift + 1);
  r = ~r;
  r <<= shift;
  return icast(r);
}

template <unsigned width>
float
PCmap<float, width, void>::identity(float d) const
{
  RANGE r = fcast(d);
  r >>= shift;
  r <<= shift;
  return icast(r);
}

template <unsigned width>
unsigned long long
PCmap<double, width, void>::fcast(double d) const
{
#ifdef WITH_REINTERPRET_CAST
  return reinterpret_cast<const RANGE&>(d);
#elif defined WITH_UNION
  UNION shared(d);
  return shared.r;
#else
  RANGE r;
  memcpy(&r, &d, sizeof(r));
  return r;
#endif
}

template <unsigned width>
double
PCmap<double, width, void>::icast(unsigned long long r) const
{
#ifdef WITH_REINTERPRET_CAST
  return reinterpret_cast<const DOMAIN&>(r);
#elif defined WITH_UNION
  UNION shared(r);
  return shared.d;
#else
  DOMAIN d;
  memcpy(&d, &r, sizeof(d));
  return d;
#endif
}

template <unsigned width>
unsigned long long
PCmap<double, width, void>::forward(double d) const
{
  RANGE r = fcast(d);
  r = ~r;
  r >>= shift;
  r ^= -(r >> (bits - 1)) >> (shift + 1);
  return r;
}

template <unsigned width>
double
PCmap<double, width, void>::inverse(unsigned long long r) const
{
  r ^= -(r >> (bits - 1)) >> (shift + 1);
  r = ~r;
  r <<= shift;
  return icast(r);
}

template <unsigned width>
double
PCmap<double, width, void>::identity(double d) const
{
  RANGE r = fcast(d);
  r >>= shift;
  r <<= shift;
  return icast(r);
}
