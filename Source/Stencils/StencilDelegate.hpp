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

enum StencilType { FGH, RHS, TURBFGH, VELOCITY, VISCOSITY };

class StencilDelegate {
public:
  StencilType                   stencilType;
  Stencils::FGHStencil          fghStencil;
  Stencils::RHSStencil          rhsStencil;
  Stencils::TurbulentFGHStencil turbulentFGHStencil;
  Stencils::VelocityStencil     velocityStencil;
  Stencils::ViscosityStencil    viscosityStencil;

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

  inline StencilType getType() const { return stencilType; }

  inline void setType(StencilType newType) { stencilType = newType; }
};
#pragma omp declare mapper(StencilDelegate fsd) map(fsd) map(fsd.stencilType \
) map(fsd.fghStencil, fsd.fghStencil.localVelocity_[0 : 81], fsd.fghStencil.localMeshsize_[0 : 81]) map(fsd.rhsStencil) map(fsd.velocityStencil \
) map(fsd.viscosityStencil, fsd.viscosityStencil.localVelocity_[0 : 81], fsd.viscosityStencil.localMeshsize_[0 : 81], fsd.viscosityStencil.coords[0 : 3]) \
  map(fsd.turbulentFGHStencil, fsd.turbulentFGHStencil.localVelocity_[0 : 81], fsd.turbulentFGHStencil.localMeshsize_[0 : 81], fsd.turbulentFGHStencil.localViscosity_[0 : 81])
