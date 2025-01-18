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


  static StencilDelegate* mapToGPU(int hostDevice, int targetDevice, StencilDelegate& stencilDelegate) {
    size_t    stencilDelegateSize = sizeof(stencilDelegate);

    StencilDelegate* stencilDelegateGPU = static_cast<StencilDelegate*>(omp_target_alloc(stencilDelegateSize, targetDevice));
    if (!stencilDelegateGPU) {
      std::cout << "Error: Allocation failed for stencilDelegateGPU pointer." << std::endl;
    }

    bool associatedStencilDelegate = omp_target_associate_ptr(&stencilDelegate, stencilDelegateGPU, stencilDelegateSize, 0, targetDevice) == 0;
    if (!associatedStencilDelegate) {
      std::cout << "Error: StencilDelegate could not be associated to GPU pointer." << std::endl;
    }

    omp_target_memcpy(&(stencilDelegateGPU->stencilType), &(stencilDelegate.stencilType), sizeof(StencilType), 0, 0, targetDevice, hostDevice);

    // Max Velocity
    { omp_target_memcpy(&(stencilDelegateGPU->maxUStencil_.maxValue_), &(stencilDelegate.maxUStencil_.maxValue_), sizeof(RealType), 0, 0, targetDevice, hostDevice); }

    // Max Visc
    { omp_target_memcpy(&(stencilDelegateGPU->maxViscStencil_.maxValue_), &(stencilDelegate.maxViscStencil_.maxValue_), sizeof(RealType), 0, 0, targetDevice, hostDevice); }

    // BFInputStencils
    {
      omp_target_memcpy(&(stencilDelegateGPU->bfInputFGHStencil_.stepSize_), &(stencilDelegate.bfInputFGHStencil_.stepSize_), sizeof(RealType), 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(
        &(stencilDelegateGPU->bfInputVelocityStencil_.stepSize_), &(stencilDelegate.bfInputVelocityStencil_.stepSize_), sizeof(RealType), 0, 0, targetDevice, hostDevice
      );
    }

    return stencilDelegateGPU;
  }

  static void freeGPU(int hostDevice, int targetDevice, StencilDelegate& stencilDelegate, StencilDelegate* stencilDelegateGPU) {

    bool disassociatedStencilDelegate = omp_target_disassociate_ptr(&stencilDelegate, targetDevice) == 0;
    if (!disassociatedStencilDelegate) {
      std::cout << "Error: Could not disassociate StencilDelegate pointer." << std::endl;
    }

    omp_target_free(stencilDelegateGPU, targetDevice);
  }
};
