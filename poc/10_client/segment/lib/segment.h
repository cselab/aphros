// Created by Sergey Litvinov on 11.03.2019
// Copyright 2019 ETH Zurich

#ifdef __cplusplus
#define EXTERNC extern "C"
#else
#define EXTERNC
#endif

typedef double Scal;
EXTERNC int segment_get(const Scal*, /**/ Scal** normal, Scal** a, Scal** seg);
EXTERNC int segment_norm(int, int, const Scal*, /**/ Scal*, Scal*);
EXTERNC int segment_ends(Scal nx, Scal ny, Scal a, /**/ Scal s[4]);
EXTERNC Scal segment_line(Scal nx, Scal ny, Scal u);
