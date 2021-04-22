// Created by Petr Karnakov on 20.11.2019
// Copyright 2019 ETH Zurich

#include "hydro.ipp"

#define XX(M)                                                                  \
  template void InitVel(                                                       \
      FieldCell<typename M::Vect>& fcv, const Vars& var, const M& m);          \
                                                                               \
  template std::tuple<                                                         \
      MapEmbed<BCondFluid<typename M::Vect>>,                                  \
      MapEmbed<BCondAdvection<typename M::Scal>>, MapEmbed<size_t>,            \
      std::vector<std::string>,                                                \
      std::vector<std::map<std::string, typename M::Scal>>>                    \
  InitBc(                                                                      \
      const Vars& var, const M& eb, std::set<std::string> known_keys,          \
      const FieldCell<bool>& fc_innermask);                                    \
                                                                               \
  template std::tuple<                                                         \
      MapEmbed<BCondFluid<typename M::Vect>>,                                  \
      MapEmbed<BCondAdvection<typename M::Scal>>, MapEmbed<size_t>,            \
      std::vector<std::string>,                                                \
      std::vector<std::map<std::string, typename M::Scal>>>                    \
  InitBc(                                                                      \
      const Vars& var, const Embed<M>& eb, std::set<std::string> known_keys,   \
      const FieldCell<bool>& fc_innermask);                                    \
                                                                               \
  template void DumpBcPoly(                                                    \
      const std::string filename, const MapEmbed<size_t>& me_group,            \
      const MapEmbed<typename M::Scal>& me_contang, const M&, M&);             \
                                                                               \
  template void DumpBcPoly(                                                    \
      const std::string filename, const MapEmbed<size_t>& me_group,            \
      const MapEmbed<typename M::Scal>& me_contang, const Embed<M>&, M&);      \
                                                                               \
  template void GetFluidCellCond(                                              \
      const Vars& var, M& m, MapCell<std::shared_ptr<CondCellFluid>>& mcvel);  \
                                                                               \
  template void InitVort(                                                      \
      const FieldCell<typename M::Vect>& fcvort,                               \
      FieldCell<typename M::Vect>& fcvel,                                      \
      const MapEmbed<BCondFluid<typename M::Vect>>& mebc_fluid,                \
      std::shared_ptr<linear::Solver<M>> linsolver, M& m);                     \
                                                                               \
  template void DumpTraj(                                                      \
      M& m, bool dm, const Vars& var, size_t frame, typename M::Scal t,        \
      const GRange<size_t>& layers,                                            \
      const Multi<const FieldCell<typename M::Scal>*>& fcvf,                   \
      const Multi<const FieldCell<typename M::Scal>*>& fccl,                   \
      const Multi<const FieldCell<typename M::MIdx>*>& fcim,                   \
      const FieldCell<typename M::Scal>& fcp,                                  \
      const FieldCell<typename M::Vect>& fcvel,                                \
      const FieldCell<typename M::Vect>& fcvelm, typename M::Scal dt);         \
  template void DumpTraj(                                                      \
      Embed<M>& m, bool dm, const Vars& var, size_t frame, typename M::Scal t, \
      const GRange<size_t>& layers,                                            \
      const Multi<const FieldCell<typename M::Scal>*>& fcvf,                   \
      const Multi<const FieldCell<typename M::Scal>*>& fccl,                   \
      const Multi<const FieldCell<typename M::MIdx>*>& fcim,                   \
      const FieldCell<typename M::Scal>& fcp,                                  \
      const FieldCell<typename M::Vect>& fcvel,                                \
      const FieldCell<typename M::Vect>& fcvelm, typename M::Scal dt);         \
                                                                               \
  template void CalcTraj(                                                      \
      M& m, const GRange<size_t>& layers,                                      \
      const Multi<const FieldCell<typename M::Scal>*>& fcvf,                   \
      const Multi<const FieldCell<typename M::Scal>*>& fccl,                   \
      const Multi<const FieldCell<typename M::MIdx>*>& fcim,                   \
      const FieldCell<typename M::Scal>& fcp,                                  \
      const FieldCell<typename M::Vect>& fcvel,                                \
      std::vector<std::string>& column_names,                                  \
      std::vector<typename M::Scal>& row_colors,                               \
      std::vector<std::vector<typename M::Scal>>& table);                      \
  template void CalcTraj(                                                      \
      Embed<M>& eb, const GRange<size_t>& layers,                              \
      const Multi<const FieldCell<typename M::Scal>*>& fcvf,                   \
      const Multi<const FieldCell<typename M::Scal>*>& fccl,                   \
      const Multi<const FieldCell<typename M::MIdx>*>& fcim,                   \
      const FieldCell<typename M::Scal>& fcp,                                  \
      const FieldCell<typename M::Vect>& fcvel,                                \
      std::vector<std::string>& column_names,                                  \
      std::vector<typename M::Scal>& row_colors,                               \
      std::vector<std::vector<typename M::Scal>>& table);                      \
                                                                               \
  template void CalcSurfaceTension(                                            \
      const M& m, const GRange<size_t>& layers, const Vars& var,               \
      FieldCell<typename M::Vect>& fc_force,                                   \
      FieldFace<typename M::Scal>& ff_force,                                   \
      const FieldCell<typename M::Scal>& fc_sig,                               \
      const MapEmbed<BCond<typename M::Scal>>& mf_sig,                         \
      const Multi<const FieldCell<typename M::Scal>*> fck,                     \
      const FieldCell<typename M::Scal>& fcvf,                                 \
      const FieldFace<typename M::Scal>& ffvfsm,                               \
      const AdvectionSolver<M>* asb);                                          \
                                                                               \
  template void ProjectVolumeFlux(                                             \
      FieldFace<typename M::Scal>& ffv,                                        \
      const MapEmbed<BCondFluid<typename M::Vect>>& mfc,                       \
      std::shared_ptr<linear::Solver<M>> linsolver, M& m);                     \
                                                                               \
  template std::map<typename M::Scal, typename M::Scal> CalcArea(              \
      const GRange<size_t>& layers,                                            \
      const Multi<const FieldCell<typename M::Vect>*> fcn,                     \
      const Multi<const FieldCell<typename M::Scal>*> fca,                     \
      const Multi<const FieldCell<typename M::Scal>*> fccl, M& m);             \
                                                                               \
  template std::map<typename M::Scal, typename M::Scal> CalcVolume(            \
      const GRange<size_t>& layers,                                            \
      const Multi<const FieldCell<typename M::Scal>*> fcu,                     \
      const Multi<const FieldCell<typename M::Scal>*> fccl, M& m);

#define COMMA ,
#define X(dim) XX(MeshCartesian<double COMMA dim>)
MULTIDIMX
