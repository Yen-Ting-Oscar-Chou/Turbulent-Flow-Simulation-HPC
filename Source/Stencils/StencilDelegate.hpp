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

enum StencilType { FGH, RHS, TURBFGH, VELOCITY, VISCOSITY, MAXU, MAXV, OBSTACLE, WALLFGH, WALLVELOCITY };

class StencilDelegate {
public:
  StencilType stencilType = FGH;

  Stencils::FGHStencil          fghStencil_;
  Stencils::RHSStencil          rhsStencil_;
  Stencils::TurbulentFGHStencil turbulentFGHStencil_;
  Stencils::VelocityStencil     velocityStencil_;
  Stencils::ViscosityStencil    viscosityStencil_;
  Stencils::MaxUStencil         maxUStencil_;
  Stencils::ObstacleStencil     obstacleStencil_;
  Stencils::MaxViscStencil      maxViscStencil_;

  Stencils::MovingWallVelocityStencil movingWallVelocityStencil_;
  Stencils::MovingWallFGHStencil      movingWallFGHStencil_;

  Stencils::NeumannVelocityBoundaryStencil neumannVelocityBoundaryStencil_;
  Stencils::NeumannFGHBoundaryStencil      neumannFGHBoundaryStencil_;

  Stencils::BFInputVelocityStencil bfInputVelocityStencil_;
  Stencils::BFInputFGHStencil      bfInputFGHStencil_;

  StencilDelegate();
  ~StencilDelegate() = default;

  inline void apply(const Parameters& parameters, FlowField& flowField, int i, int j) {
    switch (stencilType) {
    case FGH:
      fghStencil_.apply(parameters, flowField, i, j);
      break;

    case RHS:
      rhsStencil_.apply(parameters, flowField, i, j);
      break;

    case VELOCITY:
      velocityStencil_.apply(parameters, flowField, i, j);
      break;

    case OBSTACLE:
      obstacleStencil_.apply(parameters, flowField, i, j);
      break;

    default:
      assert(false);
    }
  }

  inline void apply(const Parameters& parameters, FlowField& flowField, int i, int j, int k) {
    switch (stencilType) {
    case FGH:
      fghStencil_.apply(parameters, flowField, i, j, k);
      break;

    case RHS:
      rhsStencil_.apply(parameters, flowField, i, j, k);
      break;

    case VELOCITY:
      velocityStencil_.apply(parameters, flowField, i, j, k);
      break;

    case OBSTACLE:
      obstacleStencil_.apply(parameters, flowField, i, j, k);
      break;

    default:
      assert(false);
    }
  }

  inline void apply(const Parameters& parameters, TurbulentFlowField& flowField, int i, int j) {
    switch (stencilType) {
    case VISCOSITY:
      viscosityStencil_.apply(parameters, flowField, i, j);
      break;

    case TURBFGH:
      turbulentFGHStencil_.apply(parameters, flowField, i, j);
      break;

    case MAXV:
      maxViscStencil_.apply(parameters, flowField, i, j);
      break;

    default:
      assert(false);
    }
  }

  inline void apply(const Parameters& parameters, TurbulentFlowField& flowField, int i, int j, int k) {
    switch (stencilType) {
    case VISCOSITY:
      viscosityStencil_.apply(parameters, flowField, i, j, k);
      break;

    case TURBFGH:
      turbulentFGHStencil_.apply(parameters, flowField, i, j, k);
      break;

    case MAXV:
      maxViscStencil_.apply(parameters, flowField, i, j, k);
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
        movingWallFGHStencil_.applyLeftWall(parameters, flowField, i, j, k);
      } else if (scenario == CHANNEL) {
        bfInputFGHStencil_.applyLeftWall(parameters, flowField, i, j, k);
      }
      break;
    case WALLVELOCITY:
      if (scenario == CAVITY) {
        movingWallVelocityStencil_.applyLeftWall(parameters, flowField, i, j, k);
      } else if (scenario == CHANNEL) {
        bfInputVelocityStencil_.applyLeftWall(parameters, flowField, i, j, k);
      }
      break;
    case MAXU:
      maxUStencil_.applyLeftWall(parameters, flowField, i, j, k);
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
        movingWallFGHStencil_.applyRightWall(parameters, flowField, i, j, k);
      } else if (scenario == CHANNEL) {
        neumannFGHBoundaryStencil_.applyRightWall(parameters, flowField, i, j, k);
      }
      break;
    case WALLVELOCITY:
      if (scenario == CAVITY) {
        movingWallVelocityStencil_.applyRightWall(parameters, flowField, i, j, k);
      } else if (scenario == CHANNEL) {
        neumannVelocityBoundaryStencil_.applyRightWall(parameters, flowField, i, j, k);
      }
      break;
    case MAXU:
      maxUStencil_.applyRightWall(parameters, flowField, i, j, k);
      break;
    default:
      assert(false);
    }
  }

  inline void applyTopWall(const Parameters& parameters, FlowField& flowField, int i, int j, int k) {
    ScenarioType scenario = parameters.simulation.scenario;
    switch (stencilType) {
    case WALLFGH:
      movingWallFGHStencil_.applyTopWall(parameters, flowField, i, j, k);
      break;
    case WALLVELOCITY:
      movingWallVelocityStencil_.applyTopWall(parameters, flowField, i, j, k);
      break;
    case MAXU:
      maxUStencil_.applyTopWall(parameters, flowField, i, j, k);
      break;
    default:
      assert(false);
    }
  }

  inline void applyBottomWall(const Parameters& parameters, FlowField& flowField, int i, int j, int k) {
    ScenarioType scenario = parameters.simulation.scenario;
    switch (stencilType) {
    case WALLFGH:
      movingWallFGHStencil_.applyBottomWall(parameters, flowField, i, j, k);
      break;
    case WALLVELOCITY:
      movingWallVelocityStencil_.applyBottomWall(parameters, flowField, i, j, k);
      break;
    case MAXU:
      maxUStencil_.applyBottomWall(parameters, flowField, i, j, k);
      break;
    default:
      assert(false);
    }
  }

  inline void applyFrontWall(const Parameters& parameters, FlowField& flowField, int i, int j, int k) {
    ScenarioType scenario = parameters.simulation.scenario;
    switch (stencilType) {
    case WALLFGH:
      movingWallFGHStencil_.applyFrontWall(parameters, flowField, i, j, k);
      break;
    case WALLVELOCITY:
      movingWallVelocityStencil_.applyFrontWall(parameters, flowField, i, j, k);
      break;
    case MAXU:
      maxUStencil_.applyFrontWall(parameters, flowField, i, j, k);
      break;
    default:
      assert(false);
    }
  }

  inline void applyBackWall(const Parameters& parameters, FlowField& flowField, int i, int j, int k) {
    ScenarioType scenario = parameters.simulation.scenario;
    switch (stencilType) {
    case WALLFGH:
      movingWallFGHStencil_.applyBackWall(parameters, flowField, i, j, k);
      break;
    case WALLVELOCITY:
      movingWallVelocityStencil_.applyBackWall(parameters, flowField, i, j, k);
      break;
    case MAXU:
      maxUStencil_.applyBackWall(parameters, flowField, i, j, k);
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
        movingWallFGHStencil_.applyLeftWall(parameters, flowField, i, j);
      } else if (scenario == CHANNEL) {
        bfInputFGHStencil_.applyLeftWall(parameters, flowField, i, j);
      }
      break;
    case WALLVELOCITY:
      if (scenario == CAVITY) {
        movingWallVelocityStencil_.applyLeftWall(parameters, flowField, i, j);
      } else if (scenario == CHANNEL) {
        bfInputVelocityStencil_.applyLeftWall(parameters, flowField, i, j);
      }
      break;
    case MAXU:
      maxUStencil_.applyLeftWall(parameters, flowField, i, j);
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
        movingWallFGHStencil_.applyRightWall(parameters, flowField, i, j);
      } else if (scenario == CHANNEL) {
        neumannFGHBoundaryStencil_.applyRightWall(parameters, flowField, i, j);
      }
      break;
    case WALLVELOCITY:
      if (scenario == CAVITY) {
        movingWallVelocityStencil_.applyRightWall(parameters, flowField, i, j);
      } else if (scenario == CHANNEL) {
        neumannVelocityBoundaryStencil_.applyRightWall(parameters, flowField, i, j);
      }
      break;
    case MAXU:
      maxUStencil_.applyRightWall(parameters, flowField, i, j);
      break;
    default:
      assert(false);
    }
  }

  inline void applyTopWall(const Parameters& parameters, FlowField& flowField, int i, int j) {
    ScenarioType scenario = parameters.simulation.scenario;
    switch (stencilType) {
    case WALLFGH:
      movingWallFGHStencil_.applyTopWall(parameters, flowField, i, j);
      break;
    case WALLVELOCITY:
      movingWallVelocityStencil_.applyTopWall(parameters, flowField, i, j);
      break;
    case MAXU:
      maxUStencil_.applyTopWall(parameters, flowField, i, j);
      break;
    default:
      assert(false);
    }
  }

  inline void applyBottomWall(const Parameters& parameters, FlowField& flowField, int i, int j) {
    ScenarioType scenario = parameters.simulation.scenario;
    switch (stencilType) {
    case WALLFGH:
      movingWallFGHStencil_.applyBottomWall(parameters, flowField, i, j);
      break;
    case WALLVELOCITY:
      movingWallVelocityStencil_.applyBottomWall(parameters, flowField, i, j);
      break;
    case MAXU:
      maxUStencil_.applyBottomWall(parameters, flowField, i, j);
      break;
    default:
      assert(false);
    }
  }

  inline void applyLeftWall(const Parameters& parameters, TurbulentFlowField& flowField, int i, int j) {
    switch (stencilType) {
    case MAXV:
      maxViscStencil_.applyLeftWall(parameters, flowField, i, j);
      break;
    default:
      assert(false);
    }
  }

  inline void applyRightWall(const Parameters& parameters, TurbulentFlowField& flowField, int i, int j) {
    switch (stencilType) {
    case MAXV:
      maxViscStencil_.applyRightWall(parameters, flowField, i, j);
      break;
    default:
      assert(false);
    }
  }

  inline void applyTopWall(const Parameters& parameters, TurbulentFlowField& flowField, int i, int j) {
    switch (stencilType) {
    case MAXV:
      maxViscStencil_.applyTopWall(parameters, flowField, i, j);
      break;
    default:
      assert(false);
    }
  }

  inline void applyBottomWall(const Parameters& parameters, TurbulentFlowField& flowField, int i, int j) {
    switch (stencilType) {
    case MAXV:
      maxViscStencil_.applyBottomWall(parameters, flowField, i, j);
      break;
    default:
      assert(false);
    }
  }

  inline void applyLeftWall(const Parameters& parameters, TurbulentFlowField& flowField, int i, int j, int k) {
    switch (stencilType) {
    case MAXV:
      maxViscStencil_.applyLeftWall(parameters, flowField, i, j, k);
      break;
    default:
      assert(false);
    }
  }

  inline void applyRightWall(const Parameters& parameters, TurbulentFlowField& flowField, int i, int j, int k) {
    switch (stencilType) {
    case MAXV:
      maxViscStencil_.applyRightWall(parameters, flowField, i, j, k);
      break;
    default:
      assert(false);
    }
  }

  inline void applyFrontWall(const Parameters& parameters, TurbulentFlowField& flowField, int i, int j, int k) {
    switch (stencilType) {
    case MAXV:
      maxViscStencil_.applyFrontWall(parameters, flowField, i, j, k);
      break;
    default:
      assert(false);
    }
  }

  inline void applyBackWall(const Parameters& parameters, TurbulentFlowField& flowField, int i, int j, int k) {
    switch (stencilType) {
    case MAXV:
      maxViscStencil_.applyBackWall(parameters, flowField, i, j, k);
      break;
    default:
      assert(false);
    }
  }

  inline void applyBottomWall(const Parameters& parameters, TurbulentFlowField& flowField, int i, int j, int k) {
    switch (stencilType) {
    case MAXV:
      maxViscStencil_.applyBottomWall(parameters, flowField, i, j, k);
      break;
    default:
      assert(false);
    }
  }

  inline void applyTopWall(const Parameters& parameters, TurbulentFlowField& flowField, int i, int j, int k) {
    switch (stencilType) {
    case MAXV:
      maxViscStencil_.applyTopWall(parameters, flowField, i, j, k);
      break;
    default:
      assert(false);
    }
  }

  inline StencilType getType() const { return stencilType; }

  inline void setType(StencilType newType) { stencilType = newType; }
};
#pragma omp declare mapper(StencilDelegate fsd) map(fsd) map(fsd.stencilType \
) map(fsd.fghStencil_, fsd.fghStencil_.localVelocity_[0 : 81], fsd.fghStencil_.localMeshsize_[0 : 81]) map(fsd.rhsStencil_) map(fsd.velocityStencil_ \
) map(fsd.viscosityStencil_, fsd.viscosityStencil_.localVelocity_[0 : 81], fsd.viscosityStencil_.localMeshsize_[0 : 81], fsd.viscosityStencil_.coords[0 : 3]) \
  map(fsd.turbulentFGHStencil_, fsd.turbulentFGHStencil_.localVelocity_[0 : 81], fsd.turbulentFGHStencil_.localMeshsize_[0 : 81], fsd.turbulentFGHStencil_.localViscosity_[0 : 81] \
  )
