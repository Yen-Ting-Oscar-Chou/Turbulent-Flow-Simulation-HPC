#pragma once

#include "DistanceStencil.hpp"
#include "FGHStencil.hpp"
#include "MaxUStencil.hpp"
#include "MaxViscStencil.hpp"
#include "ObstacleCoordinatesStencil.hpp"
#include "ObstacleStencil.hpp"
#include "RHSStencil.hpp"
#include "TurbulentFGHStencil.hpp"
#include "VelocityStencil.hpp"
#include "ViscosityStencil.hpp"

enum StencilType { DISTANCE, FGH, MAXVELO, MAXVISC, OBSTACLE, OBSTCOORD, RHS, TURBFGH, VELOCITY, VISCOSITY };

class FieldStencilDelegate {
public:
  StencilType stencilType;
  // Stencils::DistanceStencil            distanceStencil;
  Stencils::FGHStencil fghStencil;
  // Stencils::MaxUStencil    maxUStencil;
  // Stencils::MaxViscStencil maxViscStencil;
  // Stencils::ObstacleCoordinatesStencil obstacleCoordinatesStencil;
  Stencils::ObstacleStencil     obstacleStencil;
  Stencils::RHSStencil          rhsStencil;
  Stencils::TurbulentFGHStencil turbulentFGHStencil;
  Stencils::VelocityStencil     velocityStencil;
  Stencils::ViscosityStencil    viscosityStencil;

  FieldStencilDelegate(const StencilType type);
  ~FieldStencilDelegate() = default;

  inline void apply(const Parameters& parameters, FlowField& flowField, int i, int j) {
    switch (stencilType) {
    case FGH:
      fghStencil.apply(parameters, flowField, i, j);
      break;

    case MAXVELO:
      // maxUStencil.apply(parameters, flowField, i, j);
      assert(false);
      break;

    case OBSTACLE:
      obstacleStencil.apply(parameters, flowField, i, j);
      break;

    case RHS:
      rhsStencil.apply(parameters, flowField, i, j);
      break;

    case VELOCITY:
      velocityStencil.apply(parameters, flowField, i, j);
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

    case MAXVELO:
      // maxUStencil.apply(parameters, flowField, i, j, k);
      assert(false);
      break;

    case OBSTACLE:
      obstacleStencil.apply(parameters, flowField, i, j, k);
      break;

    case RHS:
      rhsStencil.apply(parameters, flowField, i, j, k);
      break;

    case VELOCITY:
      velocityStencil.apply(parameters, flowField, i, j, k);
      break;

    default:
      assert(false);
    }
  }

  inline void apply(const Parameters& parameters, TurbulentFlowField& flowField, int i, int j) {
    switch (stencilType) {
    case DISTANCE:
      //   distanceStencil.apply(parameters, flowField, i, j);
      assert(false);
      break;

    case MAXVISC:
      // maxViscStencil.apply(parameters, flowField, i, j);
      assert(false);
      break;

    case VISCOSITY:
      viscosityStencil.apply(parameters, flowField, i, j);
      break;

    case TURBFGH:
      turbulentFGHStencil.apply(parameters, flowField, i, j);
      break;

    case OBSTCOORD:
      assert(false);
      break;

    default:
      assert(false);
    }
  }

  inline void apply(const Parameters& parameters, TurbulentFlowField& flowField, int i, int j, int k) {
    switch (stencilType) {
    case DISTANCE:
      // distanceStencil.apply(parameters, flowField, i, j, k);
      assert(false);
      break;

    case MAXVISC:
      // maxViscStencil.apply(parameters, flowField, i, j, k);
      assert(false);
      break;

    case VISCOSITY:
      viscosityStencil.apply(parameters, flowField, i, j, k);
      break;

    case OBSTCOORD:
      assert(false);
      break;

    case TURBFGH:
      turbulentFGHStencil.apply(parameters, flowField, i, j, k);
      break;

    default:
      assert(false);
    }
  }

  inline StencilType getType() const { return stencilType; }
};
#pragma omp declare mapper(FieldStencilDelegate fsd) \
  map(fsd) \
  map(fsd.stencilType) \
  map(fsd.fghStencil) \
  map(fsd.rhsStencil) \
  map(fsd.velocityStencil) \
  map(fsd.viscosityStencil) \
  map(fsd.obstacleStencil) 
  // map(fsd.maxUStencil) \
  // map(fsd.maxViscStencil) \
  // map(fsd.turbulentFGHStencil) \