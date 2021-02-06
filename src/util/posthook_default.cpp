// Created by Petr Karnakov on 25.02.2020
// Copyright 2020 ETH Zurich

#include "posthook.h"

template <class M>
void PostHook(const Vars&, const FieldCell<typename M::Vect>&, M&) {}

template <class M>
void PostHook(
    const Vars&, const FieldCell<typename M::Vect>&, M&, const Embed<M>&) {}

template <class M>
void InitVelHook(FieldCell<typename M::Vect>&, const Vars&, const M&) {}

template <class M>
void InitVelHook(
    FieldCell<typename M::Vect>&, const Vars&, const M&, const Embed<M>&) {}

template <class M>
void InitEmbedHook(FieldNode<typename M::Scal>&, const Vars&, const M&) {}

template <class M>
void FluidDummyHook(
    FieldCell<typename M::Vect>&, FieldFace<typename M::Scal>&,
    typename M::Scal, typename M::Scal, const Vars&, const M&) {}

template <class M>
void StepHook(Hydro<M>*) {}

#define XX(M)                                                                  \
  template void PostHook(const Vars&, const FieldCell<typename M::Vect>&, M&); \
  template void PostHook(                                                      \
      const Vars&, const FieldCell<typename M::Vect>&, M&, const Embed<M>&);   \
  template void InitVelHook(                                                   \
      FieldCell<typename M::Vect>&, const Vars&, const M&);                    \
  template void InitVelHook(                                                   \
      FieldCell<typename M::Vect>&, const Vars&, const M&, const Embed<M>&);   \
  template void InitEmbedHook(                                                 \
      FieldNode<typename M::Scal>&, const Vars&, const M&);                    \
  template void FluidDummyHook(                                                \
      FieldCell<typename M::Vect>&, FieldFace<typename M::Scal>&,              \
      typename M::Scal t, typename M::Scal dt, const Vars&, const M&);         \
  template void StepHook(Hydro<M>*);

#define COMMA ,
#define X(dim) XX(MeshCartesian<double COMMA dim>)
MULTIDIMX
