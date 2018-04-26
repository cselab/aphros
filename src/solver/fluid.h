#pragma once

#include <cmath>
#include <sstream>
#include <stdexcept>
#include <memory>
#include <string>

#include "solver.h"

namespace solver {

template <class Mesh>
class FluidSolver : public UnsteadyIterativeSolver {
  static constexpr size_t dim = Mesh::dim;
  using Scal = typename Mesh::Scal;
  using Vect = typename Mesh::Vect;
  using Expr = Expression<Scal, IdxCell, 1 + dim * 2>;

 protected:
  FieldCell<Scal>* p_fc_density_;
  FieldCell<Scal>* p_fc_viscosity_;
  FieldFace<Scal>* p_ff_force_; 
  FieldCell<Scal>* p_fc_volume_source_;
  FieldCell<Scal>* p_fc_mass_source_;

 public:
  FluidSolver(double time, double time_step,
              FieldCell<Scal>* p_fc_density,
              FieldCell<Scal>* p_fc_viscosity,
              FieldFace<Scal>* p_ff_force,
              FieldCell<Scal>* p_fc_volume_source,
              FieldCell<Scal>* p_fc_mass_source)
      : UnsteadyIterativeSolver(time, time_step)
      , p_fc_density_(p_fc_density)
      , p_fc_viscosity_(p_fc_viscosity)
      , p_ff_force_(p_ff_force)
      , p_fc_volume_source_(p_fc_volume_source)
      , p_fc_mass_source_(p_fc_mass_source)
  {

  }
  virtual const FieldCell<Vect>& GetVelocity() = 0;
  virtual const FieldCell<Vect>& GetVelocity(Layers layer) = 0;
  virtual const FieldCell<Scal>& GetPressure() = 0;
  virtual const FieldCell<Scal>& GetPressure(Layers layer) = 0;
  virtual const FieldFace<Scal>& GetVolumeFlux() = 0;
  virtual const FieldFace<Scal>& GetVolumeFlux(Layers layer) = 0;
  virtual double GetAutoTimeStep() { return GetTimeStep(); }
};

class ConditionFaceFluid : public ConditionFace {
 public:
  ConditionFaceFluid(size_t nci) : ConditionFace(nci) {}
};

class ConditionCellFluid : public ConditionCell {};

namespace fluid_condition {

template <class Mesh>
class NoSlipWall : public ConditionFaceFluid {
  using Vect = typename Mesh::Vect;
 public:
  NoSlipWall(size_t nci) : ConditionFaceFluid(nci) {}
  virtual Vect GetVelocity() const = 0;
};

template <class Mesh>
class NoSlipWallFixed : public NoSlipWall<Mesh> {
  using Vect = typename Mesh::Vect;
  Vect velocity_;
 public:
  NoSlipWallFixed(Vect velocity, size_t nci)
      : NoSlipWall<Mesh>(nci), velocity_(velocity)
  {}
  Vect GetVelocity() const override {
    return velocity_;
  }
};

template <class Mesh>
class Inlet : public ConditionFaceFluid {
  using Vect = typename Mesh::Vect;
 public:
  Inlet(size_t nci) : ConditionFaceFluid(nci) {}
  virtual Vect GetVelocity() const = 0;
  virtual void SetVelocity(Vect velocity) = 0;
};

template <class Mesh>
class InletFixed : public Inlet<Mesh> {
  using Vect = typename Mesh::Vect;
  Vect velocity_;

 public:
  InletFixed(Vect velocity, size_t nci)
      : Inlet<Mesh>(nci)
      , velocity_(velocity)
  {}
  Vect GetVelocity() const override {
    return velocity_;
  }
  void SetVelocity(Vect velocity) override {
    velocity_ = velocity;
  }
};

template <class Mesh>
class InletFlux : public Inlet<Mesh> {
  using Vect = typename Mesh::Vect;
  Vect vel_;
  size_t id_;

 public:
  InletFlux(Vect vel, size_t id, size_t nci)
      : Inlet<Mesh>(nci)
      , vel_(vel)
      , id_(id)
  {}
  Vect GetVelocity() const override {
    return vel_;
  }
  void SetVelocity(Vect vel) override {
    vel_ = vel;
  }
  size_t GetId() {
    return id_;
  }
};

template <class Mesh>
class Outlet : public ConditionFaceFluid {
  using Vect = typename Mesh::Vect;
 public:
  Outlet(size_t nci) : ConditionFaceFluid(nci) {}
  virtual Vect GetVelocity() const = 0;
  virtual void SetVelocity(Vect velocity) = 0;
};

template <class Mesh>
class OutletAuto : public Outlet<Mesh> {
  using Vect = typename Mesh::Vect;
  Vect velocity_;

 public:
  OutletAuto(size_t nci)
      : Outlet<Mesh>(nci)
      , velocity_(Vect::kZero)
  {}
  Vect GetVelocity() const override {
    return velocity_;
  }
  void SetVelocity(Vect velocity) override {
    velocity_ = velocity;
  }
};

template <class Mesh>
class GivenVelocityAndPressure : public ConditionCellFluid {
  using Vect = typename Mesh::Vect;
  using Scal = typename Mesh::Scal;
 public:
  virtual Vect GetVelocity() const = 0;
  virtual Scal GetPressure() const = 0;
};

template <class Mesh>
class GivenVelocityAndPressureFixed : public GivenVelocityAndPressure<Mesh> {
  using Vect = typename Mesh::Vect;
  using Scal = typename Mesh::Scal;
  Vect velocity_;
  Scal pressure_;

 public:
  GivenVelocityAndPressureFixed(Vect velocity, Scal pressure)
      : velocity_(velocity)
      , pressure_(pressure)
  {}
  Vect GetVelocity() const override {
    return velocity_;
  }
  Scal GetPressure() const override {
    return pressure_;
  }
};

template <class Mesh>
class GivenPressure : public ConditionCellFluid {
  using Scal = typename Mesh::Scal;
 public:
  virtual Scal GetPressure() const = 0;
};

template <class Mesh>
class GivenPressureFixed : public GivenPressure<Mesh> {
  using Scal = typename Mesh::Scal;
  Scal pressure_;

 public:
  GivenPressureFixed(Scal pressure)
      : pressure_(pressure)
  {}
  Scal GetPressure() const override {
    return pressure_;
  }
};

} // namespace fluid_condition


template <class Mesh>
std::shared_ptr<ConditionFaceFluid> Parse(std::string argstr,
                                          IdxFace idxface,
                                          size_t nc, // neighbour cell id
                                          const Mesh& mesh) {
  using namespace fluid_condition;
  using Vect=  typename Mesh::Vect;
  std::stringstream arg(argstr);

  std::string name;
  arg >> name;

  if (name == "wall") {
    // No-slip wall.
    // wall <velocity>
    Vect vel;
    arg >> vel;
    return std::make_shared<NoSlipWallFixed<Mesh>>(vel, nc);
  } else if (name == "inlet") {
    // Fixed velocity inlet.
    // inlet <velocity>
    Vect vel;
    arg >> vel;
    return std::make_shared<InletFixed<Mesh>>(vel, nc);
  } else if (name == "inletflux") {
    // Fixed flux inlet. Flux defined by given velocity is redistributed
    // over all faces with same id.
    // inletflux <velocity> <id>
    Vect vel;
    int id;
    arg >> vel >> id;
    return std::make_shared<InletFlux<Mesh>>(vel, id, nc);
  } else if (name == "outlet") {
    // Outlet. Velocity is extrapolated from neighbour cells and corrected
    // to yield zero total flux over outlet and inlet faces.
    return std::make_shared<OutletAuto<Mesh>>(nc);
  } else {
    throw std::runtime_error("Parse: Unknown boundary condition type");
  }
}

} // namespace solver

