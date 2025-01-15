#pragma once

#include "BFInputStencils.hpp"
#include "DistanceStencil.hpp"
#include "FGHStencil.hpp"
#include "MaxUStencil.hpp"
#include "MaxViscStencil.hpp"
#include "MovingWallStencils.hpp"
#include "NeumannBoundaryStencils.hpp"
#include "ObstacleCoordinatesStencil.hpp"
#include "ObstacleStencil.hpp"
#include "PeriodicBoundaryStencils.hpp"
#include "RHSStencil.hpp"
#include "TurbulentFGHStencil.hpp"
#include "VelocityStencil.hpp"
#include "ViscosityStencil.hpp"

enum StencilType { FGH, RHS, TURBFGH, VELOCITY, VISCOSITY, MAXU, OBSTACLE, WALLFGH, WALLVELOCITY };

class StencilDelegate {
public:
  StencilType stencilType;

  Stencils::FGHStencil          fghStencil;
  Stencils::RHSStencil          rhsStencil;
  Stencils::TurbulentFGHStencil turbulentFGHStencil;
  Stencils::VelocityStencil     velocityStencil;
  Stencils::ViscosityStencil    viscosityStencil;
  Stencils::MaxUStencil         maxUStencil;
  Stencils::ObstacleStencil     obstacleStencil;

  Stencils::MovingWallVelocityStencil movingWallVelocityStencil;
  Stencils::MovingWallFGHStencil      movingWallFGHStencil;

  Stencils::NeumannVelocityBoundaryStencil neumannVelocityBoundaryStencil;
  Stencils::NeumannFGHBoundaryStencil      neumannFGHBoundaryStencil;

  Stencils::BFInputVelocityStencil bfInputVelocityStencil;
  Stencils::BFInputFGHStencil      bfInputFGHStencil;


  StencilDelegate(const StencilType type);
  ~StencilDelegate() = default;

  inline void apply(const Parameters& parameters, FlowField& flowField, int i, int j) {
    switch (stencilType) {
    case FGH:
      fghStencil.apply(parameters, flowField, i, j);
      break;

    case RHS:
      rhsStencil.apply(parameters, flowField, i, j);
      break;

    case VELOCITY:
      velocityStencil.apply(parameters, flowField, i, j);
      break;

    case OBSTACLE:
      obstacleStencil.apply(parameters, flowField, i, j);
      break;

    default:
      assert(false);
    }
  }

  inline void apply(const Parameters& parameters, FlowField& flowField, int i, int j, int k) {
    switch (stencilType) {
    case FGH:
      fghStencil.apply(parameters, flowField, i, j, k);
      break;

    case RHS:
      rhsStencil.apply(parameters, flowField, i, j, k);
      break;

    case VELOCITY:
      velocityStencil.apply(parameters, flowField, i, j, k);
      break;

    case OBSTACLE:
      obstacleStencil.apply(parameters, flowField, i, j, k);
      break;

    default:
      assert(false);
    }
  }

  inline void apply(const Parameters& parameters, TurbulentFlowField& flowField, int i, int j) {
    switch (stencilType) {
    case VISCOSITY:
      viscosityStencil.apply(parameters, flowField, i, j);
      break;

    case TURBFGH:
      turbulentFGHStencil.apply(parameters, flowField, i, j);
      break;

    default:
      assert(false);
    }
  }

  inline void apply(const Parameters& parameters, TurbulentFlowField& flowField, int i, int j, int k) {
    switch (stencilType) {
    case VISCOSITY:
      viscosityStencil.apply(parameters, flowField, i, j, k);
      break;

    case TURBFGH:
      turbulentFGHStencil.apply(parameters, flowField, i, j, k);
      break;

    default:
      assert(false);
    }
  }

  inline void applyLeftWall(const Parameters& parameters, FlowField& flowField, int i, int j, int k) {
    ScenarioType scenario = parameters.simulation.scenario;
    switch (stencilType) {
    case WALLFGH:
      if (scenario == CAVITY) {
        movingWallFGHStencil.applyLeftWall(parameters, flowField, i, j, k);
      } else if (scenario == CHANNEL) {
        bfInputFGHStencil.applyLeftWall(parameters, flowField, i, j, k);
      }
      break;
    case WALLVELOCITY:
      if (scenario == CAVITY) {
        movingWallVelocityStencil.applyLeftWall(parameters, flowField, i, j, k);
      } else if (scenario == CHANNEL) {
        bfInputVelocityStencil.applyLeftWall(parameters, flowField, i, j, k);
      }
      break;
    default:
      assert(false);
    }
  }

  inline void applyRightWall(const Parameters& parameters, FlowField& flowField, int i, int j, int k) {
    ScenarioType scenario = parameters.simulation.scenario;
    switch (stencilType) {
    case WALLFGH:
      if (scenario == CAVITY) {
        movingWallFGHStencil.applyRightWall(parameters, flowField, i, j, k);
      } else if (scenario == CHANNEL) {
        neumannFGHBoundaryStencil.applyRightWall(parameters, flowField, i, j, k);
      }
      break;
    case WALLVELOCITY:
      if (scenario == CAVITY) {
        movingWallVelocityStencil.applyRightWall(parameters, flowField, i, j, k);
      } else if (scenario == CHANNEL) {
        neumannVelocityBoundaryStencil.applyRightWall(parameters, flowField, i, j, k);
      }
      break;
    default:
      assert(false);
    }
  }

  inline void applyTopWall(const Parameters& parameters, FlowField& flowField, int i, int j, int k) {
    ScenarioType scenario = parameters.simulation.scenario;
    switch (stencilType) {
    case WALLFGH:
      movingWallFGHStencil.applyTopWall(parameters, flowField, i, j, k);
      break;
    case WALLVELOCITY:
      movingWallVelocityStencil.applyTopWall(parameters, flowField, i, j, k);
      break;
    default:
      assert(false);
    }
  }

  inline void applyBottomWall(const Parameters& parameters, FlowField& flowField, int i, int j, int k) {
    ScenarioType scenario = parameters.simulation.scenario;
    switch (stencilType) {
    case WALLFGH:
      movingWallFGHStencil.applyBottomWall(parameters, flowField, i, j, k);
      break;
    case WALLVELOCITY:
      movingWallVelocityStencil.applyBottomWall(parameters, flowField, i, j, k);
      break;
    default:
      assert(false);
    }
  }

  inline void applyFrontWall(const Parameters& parameters, FlowField& flowField, int i, int j, int k) {
    ScenarioType scenario = parameters.simulation.scenario;
    switch (stencilType) {
    case WALLFGH:
      movingWallFGHStencil.applyFrontWall(parameters, flowField, i, j, k);
      break;
    case WALLVELOCITY:
      movingWallVelocityStencil.applyFrontWall(parameters, flowField, i, j, k);
      break;
    default:
      assert(false);
    }
  }

  inline void applyBackWall(const Parameters& parameters, FlowField& flowField, int i, int j, int k) {
    ScenarioType scenario = parameters.simulation.scenario;
    switch (stencilType) {
    case WALLFGH:
      movingWallFGHStencil.applyBackWall(parameters, flowField, i, j, k);
      break;
    case WALLVELOCITY:
      movingWallVelocityStencil.applyBackWall(parameters, flowField, i, j, k);
      break;
    default:
      assert(false);
    }
  }

  inline void applyLeftWall(const Parameters& parameters, FlowField& flowField, int i, int j) {
    ScenarioType scenario = parameters.simulation.scenario;
    switch (stencilType) {
    case WALLFGH:
      if (scenario == CAVITY) {
        movingWallFGHStencil.applyLeftWall(parameters, flowField, i, j);
      } else if (scenario == CHANNEL) {
        bfInputFGHStencil.applyLeftWall(parameters, flowField, i, j);
      }
      break;
    case WALLVELOCITY:
      if (scenario == CAVITY) {
        movingWallVelocityStencil.applyLeftWall(parameters, flowField, i, j);
      } else if (scenario == CHANNEL) {
        bfInputVelocityStencil.applyLeftWall(parameters, flowField, i, j);
      }
      break;
    default:
      assert(false);
    }
  }

  inline void applyRightWall(const Parameters& parameters, FlowField& flowField, int i, int j) {
    ScenarioType scenario = parameters.simulation.scenario;
    switch (stencilType) {
    case WALLFGH:
      if (scenario == CAVITY) {
        movingWallFGHStencil.applyRightWall(parameters, flowField, i, j);
      } else if (scenario == CHANNEL) {
        neumannFGHBoundaryStencil.applyRightWall(parameters, flowField, i, j);
      }
      break;
    case WALLVELOCITY:
      if (scenario == CAVITY) {
        movingWallVelocityStencil.applyRightWall(parameters, flowField, i, j);
      } else if (scenario == CHANNEL) {
        neumannVelocityBoundaryStencil.applyRightWall(parameters, flowField, i, j);
      }
      break;
    default:
      assert(false);
    }
  }

  inline void applyTopWall(const Parameters& parameters, FlowField& flowField, int i, int j) {
    ScenarioType scenario = parameters.simulation.scenario;
    switch (stencilType) {
    case WALLFGH:
      movingWallFGHStencil.applyTopWall(parameters, flowField, i, j);
      break;
    case WALLVELOCITY:
      movingWallVelocityStencil.applyTopWall(parameters, flowField, i, j);
      break;
    default:
      assert(false);
    }
  }

  inline void applyBottomWall(const Parameters& parameters, FlowField& flowField, int i, int j) {
    ScenarioType scenario = parameters.simulation.scenario;
    switch (stencilType) {
    case WALLFGH:
      movingWallFGHStencil.applyBottomWall(parameters, flowField, i, j);
      break;
    case WALLVELOCITY:
      movingWallVelocityStencil.applyBottomWall(parameters, flowField, i, j);
      break;
    default:
      assert(false);
    }
  }

  inline void applyFrontWall(const Parameters& parameters, FlowField& flowField, int i, int j) {
    ScenarioType scenario = parameters.simulation.scenario;
    switch (stencilType) {
    case WALLFGH:
      movingWallFGHStencil.applyFrontWall(parameters, flowField, i, j);
      break;
    case WALLVELOCITY:
      movingWallVelocityStencil.applyFrontWall(parameters, flowField, i, j);
      break;
    default:
      assert(false);
    }
  }

  inline void applyBackWall(const Parameters& parameters, FlowField& flowField, int i, int j) {
    ScenarioType scenario = parameters.simulation.scenario;
    switch (stencilType) {
    case WALLFGH:
      movingWallFGHStencil.applyBackWall(parameters, flowField, i, j);
      break;
    case WALLVELOCITY:
      movingWallVelocityStencil.applyBackWall(parameters, flowField, i, j);
      break;
    default:
      assert(false);
    }
  }

  inline StencilType getType() const { return stencilType; }

  inline void setType(StencilType newType) { stencilType = newType; }
};
#pragma omp declare mapper(StencilDelegate fsd) map(fsd) map(fsd.stencilType \
) map(fsd.fghStencil, fsd.fghStencil.localVelocity_[0 : 81], fsd.fghStencil.localMeshsize_[0 : 81]) map(fsd.rhsStencil) map(fsd.velocityStencil \
) map(fsd.viscosityStencil, fsd.viscosityStencil.localVelocity_[0 : 81], fsd.viscosityStencil.localMeshsize_[0 : 81], fsd.viscosityStencil.coords[0 : 3]) \
  map(fsd.turbulentFGHStencil, fsd.turbulentFGHStencil.localVelocity_[0 : 81], fsd.turbulentFGHStencil.localMeshsize_[0 : 81], fsd.turbulentFGHStencil.localViscosity_[0 : 81])
