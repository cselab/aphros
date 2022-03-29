// Created by Petr Karnakov on 20.03.2021
// Copyright 2021 ETH Zurich

#include <array>
#include <iomanip>
#include <iostream>
#include <stack>
#include <string>

#include "dump/dumper.h"
#include "parse/util.h"

#include "hydro_post.h"

template <class M>
struct HydroPost<M>::Imp {
  static bool GetFieldByName(
      const Hydro<M>* hydro, std::string name, FieldCell<Scal>& fc_out,
      const M& m) {
    if (name == "p" || name == "pressure") {
      fc_out = hydro->fs_->GetPressure();
      return true;
    }
    if (name == "omz" || name == "vorticity") {
      if (hydro->eb_) {
        fc_out = UEmbed<M>::GetVortScal(
            hydro->fs_->GetVelocity(), hydro->fs_->GetVelocityCond(),
            *hydro->eb_);
      }
      fc_out = UEmbed<M>::GetVortScal(
          hydro->fs_->GetVelocity(), hydro->fs_->GetVelocityCond(), m);
      return true;
    }
    if (name == "ebvf" || name == "embed fraction") {
      FieldCell<Scal> fc(m, 1);
      if (hydro->eb_) {
        auto& eb = *hydro->eb_;
        for (auto c : m.SuCells()) {
          fc[c] = eb.GetVolumeFraction(c);
        }
      }
      fc_out = fc;
      return true;
    }
    if (name == "vf" || name == "volume fraction") {
      fc_out = hydro->as_->GetField();
      return true;
    }
    if (name == "vmagn" || name == "velocity magnitude") {
      const auto& fcvel = hydro->fs_->GetVelocity();
      FieldCell<Scal> fc(m, 0);
      for (auto c : m.SuCells()) {
        fc[c] = fcvel[c].norm();
      }
      fc_out = fc;
      return true;
    }
    if (name == "vx" || name == "velocity x") {
      fc_out = GetComponent(hydro->fs_->GetVelocity(), 0);
      return true;
    }
    if (name == "vy" || name == "velocity y") {
      fc_out = GetComponent(hydro->fs_->GetVelocity(), 1);
      return true;
    }
    if (auto& electro = hydro->electro_) {
      if (name == "elpot" || name == "electric potential") {
        fc_out = electro->GetPotential();
        return true;
      }
      if (name == "elcurx" || name == "electric current x") {
        fc_out = GetComponent(electro->GetCurrent(), 0);
        return true;
      }
      if (name == "elcury" || name == "electric current y") {
        fc_out = GetComponent(electro->GetCurrent(), 1);
        return true;
      }
      if (name == "elcurmagn" || name == "electric current magnitude") {
        const auto& fccur = electro->GetCurrent();
        FieldCell<Scal> fc(m, 0);
        for (auto c : m.SuCells()) {
          fc[c] = fccur[c].norm();
        }
        fc_out = fc;
        return true;
      }
    }
    if (auto& tracer = hydro->tracer_) {
      for (auto l : tracer->GetView().layers) {
        const auto sl = std::to_string(l);
        if (name == "tu" + sl || name == "tracer" + sl) {
          fc_out = tracer->GetVolumeFraction()[l];
          return true;
        }
      }
    }
    return false;
  }
  static FieldCell<Scal> GetField(
      const Hydro<M>* hydro, std::string name, const M& m) {
    FieldCell<Scal> fc_res(m, 0);
    const bool status = GetFieldByName(hydro, name, fc_res, m);
    if (!status && m.IsRoot()) {
      std::cerr << "Unknown field '" + name + "'\n";
    }
    return fc_res;
  };
  static void DumpFields(Hydro<M>* hydro, M& m) {
    auto sem = m.GetSem("dumpfields");
    struct {
      // containers that do not invalidate pointers on insertion
      std::stack<FieldCell<Scal>> stack_fc;
      std::stack<FieldCell<Vect>> stack_fcv;
    } * ctx(sem);
    auto& t = *ctx;
    const auto& var = hydro->var;
    if (sem("dump")) {
      if (m.IsRoot()) {
        hydro->dumper_.Report(std::cerr);
      }

      auto dl = GetWords(var.String["dumplist"]);

      auto dump = [&dl, &m](const FieldCell<Scal>& fc, std::string name) {
        if (dl.count(name)) {
          m.Dump(&fc, name);
        }
      };
      auto dumpv = [&dl, &m](
                       const FieldCell<Vect>& fc, size_t i, std::string name) {
        if (dl.count(name) && i < M::dim) {
          m.Dump(&fc, i, name);
        }
      };

      auto& fc_vel = hydro->fs_->GetVelocity();
      dumpv(fc_vel, 0, "vx");
      dumpv(fc_vel, 1, "vy");
      dumpv(fc_vel, 2, "vz");
      if (dl.count("vm")) {
        t.stack_fc.emplace(m, 0);
        auto& fc_vel_magn = t.stack_fc.top();
        for (auto c : m.Cells()) {
          fc_vel_magn[c] = fc_vel[c].norm();
        }
        dump(fc_vel_magn, "vm");
      }
      dump(hydro->fs_->GetPressure(), "p");
      dump(hydro->as_->GetField(), "vf");
      dump(hydro->fc_rho_, "rho");
      dump(hydro->fc_mu_, "mu");
      dump(hydro->fc_sig_, "sig");
      dump(hydro->fc_contang_, "contang");
      dumpv(hydro->fc_wall_dist_, 0, "walldistx");
      dumpv(hydro->fc_wall_dist_, 1, "walldisty");
      dumpv(hydro->fc_wall_dist_, 2, "walldistz");
      if (dl.count("cellcond")) {
        t.stack_fc.emplace(m, 0);
        auto& fc_cellcond = t.stack_fc.top();
        for (auto& it : hydro->mc_velcond_) {
          fc_cellcond[it.first] = 1;
        }
        m.Dump(&fc_cellcond, "cellcond");
      }
      if (dl.count("blockid")) {
        t.stack_fc.emplace(m, m.GetId());
        auto& fc_blockid = t.stack_fc.top();
        m.Dump(&fc_blockid, "blockid");
      }
      if (dl.count("omx") || dl.count("omy") || dl.count("omz") ||
          dl.count("omm") || dl.count("omcalc")) {
        auto& fcom = hydro->fcom_;
        if (M::dim == 3) {
          hydro->CalcVort();
          dumpv(fcom, 0, "omx");
          dumpv(fcom, 1, "omy");
          dumpv(fcom, 2, "omz");
        } else {
          fcom.Reinit(m, Vect(0));
          if (hydro->eb_) {
            SetComponent(
                fcom, 0,
                UEmbed<M>::GetVortScal(
                    hydro->fs_->GetVelocity(), hydro->fs_->GetVelocityCond(),
                    *hydro->eb_));
          } else {
            SetComponent(
                fcom, 0,
                UEmbed<M>::GetVortScal(
                    hydro->fs_->GetVelocity(), hydro->fs_->GetVelocityCond(),
                    m));
          }
          dumpv(fcom, 0, "omz");
        }
        dump(hydro->fcomm_, "omm");
      }
      if (dl.count("fluxx") || dl.count("fluxy") || dl.count("fluxz")) {
        t.stack_fcv.emplace(m, Vect(0));
        auto& fc_flux = t.stack_fcv.top();
        auto& ffv = hydro->fs_->GetVolumeFlux();
        for (auto c : m.Cells()) {
          for (auto d : m.dirs) {
            fc_flux[c][d] = ffv[m.GetFace(c, IdxNci(d * 2))];
          }
        }
        dumpv(fc_flux, 0, "fluxx");
        dumpv(fc_flux, 1, "fluxy");
        dumpv(fc_flux, 2, "fluxz");
      }
      if (dl.count("dis") || dl.count("strain")) {
        hydro->fc_strain_ = hydro->CalcStrain(hydro->fs_->GetVelocity());
        if (dl.count("strain")) {
          m.Dump(&hydro->fc_strain_, "strain");
        }
        if (dl.count("dis")) {
          t.stack_fc.push(hydro->fc_strain_);
          auto& fc_dis = t.stack_fc.top();
          for (auto c : m.Cells()) {
            fc_dis[c] *= 2. * hydro->fc_mu_[c];
          }
          m.Dump(&fc_dis, "dis");
        }
      }
      if (dl.count("div")) {
        t.stack_fc.push(hydro->GetDiv());
        m.Dump(&t.stack_fc.top(), "div");
      }
      using ASV = typename Hydro<M>::ASV;
      using ASVEB = typename Hydro<M>::ASVEB;
      using ASVM = typename Hydro<M>::ASVM;
      if (auto as = dynamic_cast<ASV*>(hydro->as_.get())) {
        dumpv(as->GetNormal(), 0, "nx");
        dumpv(as->GetNormal(), 1, "ny");
        dumpv(as->GetNormal(), 2, "nz");
        dump(as->GetColor(), "cls");
        dump(hydro->fck_[0], "k");
      }
      // TODO reuse ASV code
      if (auto as = dynamic_cast<ASVEB*>(hydro->as_.get())) {
        dumpv(as->GetNormal(), 0, "nx");
        dumpv(as->GetNormal(), 1, "ny");
        dumpv(as->GetNormal(), 2, "nz");
        dump(as->GetColor(), "cls");
        dump(hydro->fck_[0], "k");
      }
      if (auto as = dynamic_cast<ASVM*>(hydro->as_.get())) {
        for (auto l : hydro->layers) {
          auto sl = std::to_string(l);
          dump(*as->GetFieldM()[l], "vf" + sl);
          dump(*as->GetColor()[l], "cl" + sl);
          dumpv(*as->GetNormal()[l], 0, "nx" + sl);
          dumpv(*as->GetNormal()[l], 1, "ny" + sl);
          dumpv(*as->GetNormal()[l], 2, "nz" + sl);
          dump(hydro->fck_[l], "k" + sl);
        }

        // combined colors
        dump(as->GetColorSum(), "cls");

        // image vector
        for (auto l : hydro->layers) {
          t.stack_fcv.emplace(m, Vect(0));
          auto& fcim = t.stack_fcv.top();
          auto& fcim_midx = *as->GetImage()[l];
          for (auto c : m.Cells()) {
            fcim[c] = Vect(fcim_midx[c]);
          }
          auto sl = std::to_string(l);
          dumpv(fcim, 0, "imx" + sl);
          dumpv(fcim, 1, "imy" + sl);
          dumpv(fcim, 2, "imz" + sl);
        }
      }
      // TODO add ASVMEB

      if (hydro->eb_) {
        auto& eb = *hydro->eb_;
        if (dl.count("ebvf")) {
          t.stack_fc.emplace(m, 0);
          auto& fc_ebvf = t.stack_fc.top();
          for (auto c : eb.Cells()) {
            fc_ebvf[c] = eb.GetVolumeFraction(c);
          }
          m.Dump(&fc_ebvf, "ebvf");
        }
        // Face area.
        if (dl.count("ebsx") || dl.count("ebsy") || dl.count("ebsz")) {
          t.stack_fcv.emplace(m, Vect(0));
          auto& fcs = t.stack_fcv.top();
          for (auto c : m.Cells()) {
            for (auto d : m.dirs) {
              const IdxFace f = m.GetFace(c, IdxNci(d * 2));
              fcs[c][d] = eb.GetArea(f);
            }
          }
          dumpv(fcs, 0, "ebsx");
          dumpv(fcs, 1, "ebsy");
          dumpv(fcs, 2, "ebsz");
        }
      }

      if (auto& tracer = hydro->tracer_) {
        for (auto l : tracer->GetView().layers) {
          dump(tracer->GetVolumeFraction()[l], "tu" + std::to_string(l));
        }
        if (dl.count("tusum")) {
          t.stack_fc.emplace(m, 0);
          auto& fc_tracer_sum = t.stack_fc.top();
          for (auto l : tracer->GetView().layers) {
            if (l > 0) {
              const auto& fc = tracer->GetVolumeFraction()[l];
              for (auto c : m.Cells()) {
                fc_tracer_sum[c] += fc[c];
              }
            }
          }
          dump(fc_tracer_sum, "tusum");
        }
      }
      if (auto& electro = hydro->electro_) {
        dump(electro->GetPotential(), "elpot");
        dumpv(electro->GetCurrent(), 0, "elcurx");
        dumpv(electro->GetCurrent(), 1, "elcury");
        dumpv(electro->GetCurrent(), 2, "elcurz");
        if (dl.count("elcurfx") || dl.count("elcurfy") || dl.count("elcurfz")) {
          t.stack_fcv.emplace(m, Vect(0));
          auto& fccur = t.stack_fcv.top();
          auto& ffcur = electro->GetFaceCurrent();
          for (auto c : m.Cells()) {
            for (auto d : m.dirs) {
              fccur[c][d] = ffcur[m.GetFace(c, IdxNci(d * 2))];
            }
          }
          dumpv(fccur, 0, "elcurfx");
          dumpv(fccur, 1, "elcurfy");
          dumpv(fccur, 2, "elcurfz");
        }
      }
      if (hydro->particles_) {
        dumpv(hydro->fc_momentum_part_, 0, "mompartx");
        dumpv(hydro->fc_momentum_part_, 1, "momparty");
        dumpv(hydro->fc_momentum_part_, 2, "mompartz");
      }
    }
    if (var.Int["enable_advection"]) {
      if (var.Int["dumppoly"] && sem.Nested()) {
        std::vector<Multi<const FieldCell<Scal>*>> extra_fields;
        std::vector<std::string> extra_names;
        if (var.Int("dumppoly_curv", false)) {
          extra_fields.push_back(hydro->fck_);
          extra_names.push_back("k");
        }
        hydro->as_->DumpInterface(
            GetDumpName("s", ".vtk", hydro->dumper_.GetN()), extra_fields,
            extra_names);
      }
      if (var.Int["dumppolymarch"] && sem.Nested()) {
        hydro->as_->DumpInterfaceMarch(
            GetDumpName("sm", ".vtk", hydro->dumper_.GetN()));
      }
    }
    if (hydro->particles_ && var.Int["dump_particles"]) {
      const bool dumpvtk = var.Int["particles_dumpvtk"];
      const std::string path =
          GetDumpName("part", dumpvtk ? ".vtk" : ".csv", hydro->dumper_.GetN());
      if (sem()) {
        if (m.IsRoot()) {
          std::cerr << std::fixed << std::setprecision(8) << "dump"
                    << " t=" << hydro->particles_->GetTime() << " to " << path
                    << std::endl;
        }
      }
      if (sem.Nested()) {
        hydro->particles_->DumpParticles(
            path, GetWords(var.String["particles_dumplist"]), dumpvtk);
      }
    }
    if (sem()) {
      if (m.IsRoot() && var.Int("create_dumpdone", false)) {
        std::ofstream(GetDumpName(".dumpdone", "", hydro->dumper_.GetN()));
      }
      // XXX: empty stage, otherwise ctx is destroyed before dump
    }
  }
};

template <class M>
auto HydroPost<M>::GetField(const Hydro<M>* hydro, std::string name, const M& m)
    -> FieldCell<Scal> {
  return Imp::GetField(hydro, name, m);
}

template <class M>
void HydroPost<M>::DumpFields(Hydro<M>* hydro, M& m) {
  Imp::DumpFields(hydro, m);
}
