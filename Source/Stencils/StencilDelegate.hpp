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

  inline void applyRHS(const Parameters& parameters, ScalarField& rhs, VectorField& fgh, int i, int j) {
    rhsStencil.applyGPU(parameters, rhs, fgh, i, j);
  }

  inline void applyRHS(const Parameters& parameters, ScalarField& rhs, VectorField& fgh, int i, int j, int k) {
    rhsStencil.applyGPU(parameters, rhs, fgh, i, j, k);
  }

  inline void applyFGH(const Parameters& parameters, VectorField& velocity, VectorField& fgh, IntScalarField& flags, int i, int j) {
    fghStencil.applyGPU(parameters, velocity, fgh, flags, i, j);
  }

  inline void applyFGH(const Parameters& parameters, VectorField& velocity, VectorField& fgh, IntScalarField& flags, int i, int j, int k) {
    fghStencil.applyGPU(parameters, velocity, fgh, flags, i, j, k);
  }

  inline void applyVelocity(const Parameters& parameters, VectorField& velocity, VectorField& fgh, ScalarField& pressure, IntScalarField& flags, int i, int j) {
    velocityStencil.applyGPU(parameters, velocity, fgh, pressure, flags, i, j);
  }

  inline void applyVelocity(const Parameters& parameters, VectorField& velocity, VectorField& fgh, ScalarField& pressure, IntScalarField& flags, int i, int j, int k) {
    velocityStencil.applyGPU(parameters, velocity, fgh, pressure, flags, i, j, k);
  }

  inline void applyViscosity(const Parameters& parameters, VectorField& velocity, ScalarField& viscosity, ScalarField& distance, IntScalarField& flags, int i, int j) {
    viscosityStencil.applyGPU(parameters, velocity, viscosity, distance, flags, i, j);
  }

  inline void applyViscosity(const Parameters& parameters, VectorField& velocity, ScalarField& viscosity, ScalarField& distance, IntScalarField& flags, int i, int j, int k) {
    viscosityStencil.applyGPU(parameters, velocity, viscosity, distance, flags, i, j, k);
  }

  inline void applyTurbulentFGH(const Parameters& parameters, VectorField& fgh, ScalarField& viscosity, VectorField& velocity, IntScalarField& flags, int i, int j) {
    turbulentFGHStencil.applyGPU(parameters, fgh, viscosity, velocity, flags, i, j);
  }

  inline void applyTurbulentFGH(const Parameters& parameters, VectorField& fgh, ScalarField& viscosity, VectorField& velocity, IntScalarField& flags, int i, int j, int k) {
    turbulentFGHStencil.applyGPU(parameters, fgh, viscosity, velocity, flags, i, j, k);
  }

  inline StencilType getType() const { return stencilType; }
};
#pragma omp declare mapper(StencilDelegate fsd) map(fsd) \
  map(fsd.stencilType) \
  map(fsd.fghStencil, fsd.fghStencil.localVelocity_[0 : 81], fsd.fghStencil.localMeshsize_[0 : 81]) \
  map(fsd.rhsStencil) map(fsd.velocityStencil) \
  map(fsd.viscosityStencil, fsd.viscosityStencil.localVelocity_[0 : 81], fsd.viscosityStencil.localMeshsize_[0 : 81], fsd.viscosityStencil.coords[0 : 3]) \
  map(fsd.turbulentFGHStencil, fsd.turbulentFGHStencil.localVelocity_[0 : 81], fsd.turbulentFGHStencil.localMeshsize_[0 : 81], fsd.turbulentFGHStencil.localViscosity_[0 : 81])
