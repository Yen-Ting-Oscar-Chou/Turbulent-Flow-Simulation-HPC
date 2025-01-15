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

class StencilDelegate;
struct StencilDelegatePrts {
  StencilDelegate* stencilDelegateGPU_;

  RealType* fghLocalVelocityGPU_;
  RealType* fghLocalMeshsizeGPU_;

  RealType* turbFGHLocalVelocityGPU_;
  RealType* turbFGHLocalMeshsizeGPU_;
  RealType* turbFGHLocalViscosityGPU_;

  RealType* viscLocalVelocityGPU_;
  RealType* viscLocalMeshsizeGPU_;
  RealType* viscCoordsGPU_;

  RealType* maxUMaxValuesGPU_;

  StencilDelegatePrts(
    StencilDelegate* stencilDelegateGPU,
    RealType*        fghLocalVelocityGPU,
    RealType*        fghLocalMeshsizeGPU,
    RealType*        turbFGHLocalVelocityGPU,
    RealType*        turbFGHLocalMeshsizeGPU,
    RealType*        turbFGHLocalViscosityGPU,
    RealType*        viscLocalVelocityGPU,
    RealType*        viscLocalMeshsizeGPU,
    RealType*        viscCoordsGPU,
    RealType*        maxUMaxValuesGPU
  ):
    stencilDelegateGPU_(stencilDelegateGPU),
    fghLocalVelocityGPU_(fghLocalVelocityGPU),
    fghLocalMeshsizeGPU_(fghLocalMeshsizeGPU),
    turbFGHLocalVelocityGPU_(turbFGHLocalVelocityGPU),
    turbFGHLocalMeshsizeGPU_(turbFGHLocalMeshsizeGPU),
    turbFGHLocalViscosityGPU_(turbFGHLocalViscosityGPU),
    viscLocalVelocityGPU_(viscLocalVelocityGPU),
    viscLocalMeshsizeGPU_(viscLocalMeshsizeGPU),
    viscCoordsGPU_(viscCoordsGPU),
    maxUMaxValuesGPU_(maxUMaxValuesGPU) {}
};

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


  static StencilDelegatePrts mapToGPU(int hostDevice, int targetDevice, StencilDelegate& stencilDelegate) {
    size_t    stencilDelegateSize = sizeof(stencilDelegate);
    RealType* fghLocalVelocityGPU;
    RealType* fghLocalMeshsizeGPU;
    RealType* turbFGHLocalVelocityGPU;
    RealType* turbFGHLocalMeshsizeGPU;
    RealType* turbFGHLocalViscosityGPU;
    RealType* viscLocalVelocityGPU;
    RealType* viscLocalMeshsizeGPU;
    RealType* viscCoordsGPU;
    RealType* maxUMaxValuesGPU;

    StencilDelegate* stencilDelegateGPU = static_cast<StencilDelegate*>(omp_target_alloc(stencilDelegateSize, targetDevice));
    if (!stencilDelegateGPU) {
      std::cout << "Error: Allocation failed for stencilDelegateGPU pointer." << std::endl;
    }

    bool associatedStencilDelegate = omp_target_associate_ptr(&stencilDelegate, stencilDelegateGPU, stencilDelegateSize, 0, targetDevice) == 0;
    if (!associatedStencilDelegate) {
      std::cout << "Error: StencilDelegate could not be associated to GPU pointer." << std::endl;
    }

    omp_target_memcpy(&(stencilDelegateGPU->stencilType), &(stencilDelegate.stencilType), sizeof(StencilType), 0, 0, targetDevice, hostDevice);

    // FGH
    {
      size_t arraySize = 27 * 3;

      fghLocalVelocityGPU = static_cast<RealType*>(omp_target_alloc(arraySize, targetDevice));
      if (!fghLocalVelocityGPU) {
        std::cout << "Error: Allocation failed for fghLocalVelocityGPU pointer." << std::endl;
      }
      fghLocalMeshsizeGPU = static_cast<RealType*>(omp_target_alloc(arraySize, targetDevice));
      if (!fghLocalMeshsizeGPU) {
        std::cout << "Error: Allocation failed for fghLocalMeshsizeGPU pointer." << std::endl;
      }

      bool copiedLocaVelocity = omp_target_memcpy(&stencilDelegateGPU->fghStencil_.localVelocity_, &fghLocalVelocityGPU, sizeof(RealType*), 0, 0, targetDevice, hostDevice) == 0;
      if (!copiedLocaVelocity) {
        std::cout << "Error: Could not copy localVelocity pointer to GPU (FGH)." << std::endl;
      }
      bool copiedLocalMeshsize = omp_target_memcpy(&stencilDelegateGPU->fghStencil_.localMeshsize_, &fghLocalMeshsizeGPU, sizeof(RealType*), 0, 0, targetDevice, hostDevice) == 0;
      if (!copiedLocalMeshsize) {
        std::cout << "Error: Could not copy localMeshsize pointer to GPU (FGH)." << std::endl;
      }
    }

    // TurbFGH
    {
      size_t arraySize = 27 * 3;

      turbFGHLocalVelocityGPU = static_cast<RealType*>(omp_target_alloc(arraySize, targetDevice));
      if (!turbFGHLocalVelocityGPU) {
        std::cout << "Error: Allocation failed for turbFGHLocalVelocityGPU pointer." << std::endl;
      }
      turbFGHLocalMeshsizeGPU = static_cast<RealType*>(omp_target_alloc(arraySize, targetDevice));
      if (!turbFGHLocalMeshsizeGPU) {
        std::cout << "Error: Allocation failed for turbFGHLocalMeshsizeGPU pointer." << std::endl;
      }

      turbFGHLocalViscosityGPU = static_cast<RealType*>(omp_target_alloc(arraySize, targetDevice));
      if (!turbFGHLocalViscosityGPU) {
        std::cout << "Error: Allocation failed for turbFGHLocalViscosityGPU pointer." << std::endl;
      }

      bool copiedLocaVelocity = omp_target_memcpy(
                                  &stencilDelegateGPU->turbulentFGHStencil_.localVelocity_, &turbFGHLocalVelocityGPU, sizeof(RealType*), 0, 0, targetDevice, hostDevice
                                )
                                == 0;
      if (!copiedLocaVelocity) {
        std::cout << "Error: Could not copy localVelocity pointer to GPU (turbFGH)." << std::endl;
      }
      bool copiedLocalMeshsize = omp_target_memcpy(
                                   &stencilDelegateGPU->turbulentFGHStencil_.localMeshsize_, &turbFGHLocalMeshsizeGPU, sizeof(RealType*), 0, 0, targetDevice, hostDevice
                                 )
                                 == 0;
      if (!copiedLocalMeshsize) {
        std::cout << "Error: Could not copy localMeshsize pointer to GPU (turbFGH)." << std::endl;
      }
      bool copiedLocalViscosity = omp_target_memcpy(
                                    &stencilDelegateGPU->turbulentFGHStencil_.localViscosity_, &turbFGHLocalViscosityGPU, sizeof(RealType*), 0, 0, targetDevice, hostDevice
                                  )
                                  == 0;
      if (!copiedLocalViscosity) {
        std::cout << "Error: Could not copy localViscosity pointer to GPU (turbFGH)." << std::endl;
      }
    }

    // Viscosity
    {
      size_t localArraySize = 27 * 3;
      size_t coordsSize     = 3;
      viscLocalVelocityGPU  = static_cast<RealType*>(omp_target_alloc(localArraySize, targetDevice));
      if (!fghLocalVelocityGPU) {
        std::cout << "Error: Allocation failed for viscLocalVelocityGPU pointer." << std::endl;
      }
      viscLocalMeshsizeGPU = static_cast<RealType*>(omp_target_alloc(localArraySize, targetDevice));
      if (!viscLocalMeshsizeGPU) {
        std::cout << "Error: Allocation failed for viscLocalMeshsizeGPU pointer." << std::endl;
      }
      viscCoordsGPU = static_cast<RealType*>(omp_target_alloc(coordsSize, targetDevice));
      if (!viscLocalMeshsizeGPU) {
        std::cout << "Error: Allocation failed for viscCoordsGPU pointer." << std::endl;
      }

      bool copiedLocaVelocity = omp_target_memcpy(&stencilDelegateGPU->viscosityStencil_.localVelocity_, &viscLocalVelocityGPU, sizeof(RealType*), 0, 0, targetDevice, hostDevice)
                                == 0;
      if (!copiedLocaVelocity) {
        std::cout << "Error: Could not copy localVelocity pointer to GPU (Viscosity)." << std::endl;
      }
      bool copiedLocalMeshsize = omp_target_memcpy(&stencilDelegateGPU->viscosityStencil_.localMeshsize_, &viscLocalMeshsizeGPU, sizeof(RealType*), 0, 0, targetDevice, hostDevice)
                                 == 0;
      if (!copiedLocalMeshsize) {
        std::cout << "Error: Could not copy localMeshsize pointer to GPU (Viscosity)." << std::endl;
      }
      bool copiedCoords = omp_target_memcpy(&stencilDelegateGPU->viscosityStencil_.coords, &viscCoordsGPU, sizeof(RealType*), 0, 0, targetDevice, hostDevice) == 0;
      if (!copiedLocalMeshsize) {
        std::cout << "Error: Could not copy coords pointer to GPU (Viscosity)." << std::endl;
      }
    }

    // MaxU
    {
      size_t arraySize = 3;
      // RealType* maxUMaxValuesGPU;
      maxUMaxValuesGPU = static_cast<RealType*>(omp_target_alloc(arraySize, targetDevice));
      if (!maxUMaxValuesGPU) {
        std::cout << "Error: Allocation failed for maxUMaxValuesGPU pointer." << std::endl;
      }
      bool copiedCoords = omp_target_memcpy(&stencilDelegateGPU->maxUStencil_.maxValues_, &maxUMaxValuesGPU, sizeof(RealType*), 0, 0, targetDevice, hostDevice) == 0;
      if (!copiedCoords) {
        std::cout << "Error: Could not copy max values pointer to GPU (MaxU)." << std::endl;
      }
    }

    // Max Visc
    {
      omp_target_memcpy(&(stencilDelegateGPU->maxViscStencil_.maxValue_), &(stencilDelegate.maxViscStencil_.maxValue_), sizeof(RealType), 0, 0, targetDevice, hostDevice);
    }

    // BFInputStencils
    {
      omp_target_memcpy(&(stencilDelegateGPU->bfInputFGHStencil_.stepSize_), &(stencilDelegate.bfInputFGHStencil_.stepSize_), sizeof(RealType), 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&(stencilDelegateGPU->bfInputVelocityStencil_.stepSize_), &(stencilDelegate.bfInputVelocityStencil_.stepSize_), sizeof(RealType), 0, 0, targetDevice, hostDevice);
    }

    return StencilDelegatePrts(
      stencilDelegateGPU,
      fghLocalVelocityGPU,
      fghLocalMeshsizeGPU,
      turbFGHLocalVelocityGPU,
      turbFGHLocalMeshsizeGPU,
      turbFGHLocalViscosityGPU,
      viscLocalVelocityGPU,
      viscLocalMeshsizeGPU,
      viscCoordsGPU,
      maxUMaxValuesGPU
    );
  }

  static void mapToCPU(int hostDevice, int targetDevice, StencilDelegatePrts& simulationPtrs) {
    return;
  }

  static void mapToCPUAndFree(int hostDevice, int targetDevice, StencilDelegate& stencilDelegate, StencilDelegatePrts& simulationPtrs) {

    bool disassociatedStencilDelegate = omp_target_disassociate_ptr(&stencilDelegate, targetDevice) == 0;
    if (!disassociatedStencilDelegate) {
      std::cout << "Error: Could not disassociate StencilDelegate pointer." << std::endl;
    }

    omp_target_free(simulationPtrs.fghLocalVelocityGPU_, targetDevice);
    omp_target_free(simulationPtrs.fghLocalMeshsizeGPU_, targetDevice);
    omp_target_free(simulationPtrs.turbFGHLocalVelocityGPU_, targetDevice);
    omp_target_free(simulationPtrs.turbFGHLocalMeshsizeGPU_, targetDevice);
    omp_target_free(simulationPtrs.turbFGHLocalViscosityGPU_, targetDevice);
    omp_target_free(simulationPtrs.viscLocalVelocityGPU_, targetDevice);
    omp_target_free(simulationPtrs.viscLocalMeshsizeGPU_, targetDevice);
    omp_target_free(simulationPtrs.viscCoordsGPU_, targetDevice);
    omp_target_free(simulationPtrs.maxUMaxValuesGPU_, targetDevice);
    omp_target_free(simulationPtrs.stencilDelegateGPU_, targetDevice);
  }
};
