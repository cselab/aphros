// Created by Petr Karnakov on 11.02.2020
// Copyright 2020 ETH Zurich

#include "approx_eb.ipp"
#include "linear/linear.h"

#define XX(M)                                                               \
  template struct UEmbed<M>;                                                \
  template struct ULinearFit<typename M::Vect>;                             \
                                                                            \
  template FieldCell<typename M::Scal> UEmbed<M>::Interpolate(              \
      const FieldEmbed<typename M::Scal>& feu, const Embed<M>& eb);         \
  template FieldCell<typename M::Scal> UEmbed<M>::Interpolate(              \
      const FieldFace<typename M::Scal>& feu, const M&);                    \
                                                                            \
  template FieldCell<typename M::Scal> UEmbed<M>::AverageCutCells(          \
      const FieldCell<typename M::Scal>& fcu, const Embed<M>& eb);          \
  template FieldCell<typename M::Vect> UEmbed<M>::AverageCutCells(          \
      const FieldCell<typename M::Vect>& fcu, const Embed<M>& eb);          \
                                                                            \
  template FieldCell<typename M::Scal> UEmbed<M>::RedistributeCutCells(     \
      const FieldCell<typename M::Scal>& fcu, const Embed<M>& eb);          \
  template FieldCell<typename M::Scal> UEmbed<M>::RedistributeCutCells(     \
      const FieldCell<typename M::Scal>& fcu, const M& m);                  \
                                                                            \
  template FieldFace<typename M::Scal> UEmbed<M>::InterpolateBilinearFaces( \
      const FieldFace<typename M::Scal>& ffu, const Embed<M>& eb);          \
  template FieldFace<typename M::Scal> UEmbed<M>::InterpolateBilinearFaces( \
      const FieldFace<typename M::Scal>& ffu, const M& m);                  \
                                                                            \
  template FieldEmbed<typename M::Scal> UEmbed<M>::Interpolate(             \
      const FieldCell<typename M::Scal>& fcu,                               \
      const MapEmbed<BCond<typename M::Scal>>& mebc, const Embed<M>& eb);   \
  template FieldEmbed<typename M::Vect> UEmbed<M>::Interpolate(             \
      const FieldCell<typename M::Vect>& fcu,                               \
      const MapEmbed<BCond<typename M::Vect>>& mebc, const Embed<M>& eb);   \
  template FieldEmbed<typename M::Scal> UEmbed<M>::Gradient(                \
      const FieldCell<typename M::Scal>& fcu,                               \
      const MapEmbed<BCond<typename M::Scal>>& mebc, const Embed<M>& eb);   \
  template FieldEmbed<typename M::Vect> UEmbed<M>::Gradient(                \
      const FieldCell<typename M::Vect>& fcu,                               \
      const MapEmbed<BCond<typename M::Vect>>& mebc, const Embed<M>& eb);   \
                                                                            \
  template FieldFace<typename M::Scal> UEmbed<M>::Interpolate(              \
      const FieldCell<typename M::Scal>& fcu,                               \
      const MapEmbed<BCond<typename M::Scal>>& mebc, const M& m);           \
  template FieldFace<typename M::Vect> UEmbed<M>::Interpolate(              \
      const FieldCell<typename M::Vect>& fcu,                               \
      const MapEmbed<BCond<typename M::Vect>>& mebc, const M& m);           \
  template FieldFace<typename M::Scal> UEmbed<M>::Gradient(                 \
      const FieldCell<typename M::Scal>& fcu,                               \
      const MapEmbed<BCond<typename M::Scal>>& mebc, const M& m);           \
  template FieldFace<typename M::Vect> UEmbed<M>::Gradient(                 \
      const FieldCell<typename M::Vect>& fcu,                               \
      const MapEmbed<BCond<typename M::Vect>>& mebc, const M& m);           \
                                                                            \
  template FieldCell<typename M::Vect> UEmbed<M>::GetVort(                  \
      const FieldCell<typename M::Vect>& fcv,                               \
      const MapEmbed<BCond<typename M::Vect>>& mebc, const M& m);           \
  template FieldCell<typename M::Vect> UEmbed<M>::GetVort(                  \
      const FieldCell<typename M::Vect>& fcv,                               \
      const MapEmbed<BCond<typename M::Vect>>& mebc, const Embed<M>& eb);   \
  template FieldCell<typename M::Scal> UEmbed<M>::GetVortScal(              \
      const FieldCell<typename M::Vect>& fcv,                               \
      const MapEmbed<BCond<typename M::Vect>>& mebc, const M& m);           \
  template FieldCell<typename M::Scal> UEmbed<M>::GetVortScal(              \
      const FieldCell<typename M::Vect>& fcv,                               \
      const MapEmbed<BCond<typename M::Vect>>& mebc, const Embed<M>& eb);   \
                                                                            \
  template void Smoothen(                                                   \
      FieldCell<typename M::Scal>& fc,                                      \
      const MapEmbed<BCond<typename M::Scal>>& mfc, Embed<M>& eb,           \
      size_t iters);                                                        \
  template void Smoothen(                                                   \
      FieldCell<typename M::Scal>& fc,                                      \
      const MapEmbed<BCond<typename M::Scal>>& mfc, M& eb, size_t iters);   \
                                                                            \
  template MapEmbed<BCond<typename M::Scal>> GetScalarCond(                 \
      const MapEmbed<BCond<typename M::Vect>>& mev, size_t d, const M& m);

#define COMMA ,
#define X(dim) XX(MeshCartesian<double COMMA dim>)
MULTIDIMX
