#include "StdAfx.hpp"

#include "TurbulentFGHStencil.hpp"

#include "Definitions.hpp"
#include "StencilFunctions.hpp"

void Stencils::TurbulentFGHStencil::apply(const Parameters& parameters, TurbulentFlowField& flowField, int i, int j) {
  VectorField& velocity  = flowField.getVelocity();
  ScalarField& viscosity = flowField.getViscosity();
  // Now the localVelocity array should contain lexicographically ordered elements around the given index
  flowField.getFGH().getVectorElement(i, j, 0) = computeTurbulentF2D(velocity, viscosity, parameters, parameters.timestep.dt, i, j);
  flowField.getFGH().getVectorElement(i, j, 1) = computeTurbulentG2D(velocity, viscosity, parameters, parameters.timestep.dt, i, j);
}

void Stencils::TurbulentFGHStencil::apply(const Parameters& parameters, TurbulentFlowField& flowField, int i, int j, int k) {
  // The same as in 2D, with slight modifications.

  const int obstacle = flowField.getFlags().getValue(i, j, k);

  if ((obstacle & OBSTACLE_SELF) == 0) { // If the cell is fluid
    VectorField& velocity  = flowField.getVelocity();
    ScalarField& viscosity = flowField.getViscosity();

    if ((obstacle & OBSTACLE_RIGHT) == 0) { // If the right cell is fluid
      flowField.getFGH().getVectorElement(i, j, k, 0) = computeTurbulentF3D(velocity, viscosity, parameters, parameters.timestep.dt, i, j, k);
    }
    if ((obstacle & OBSTACLE_TOP) == 0) {
      flowField.getFGH().getVectorElement(i, j, k, 1) = computeTurbulentG3D(velocity, viscosity, parameters, parameters.timestep.dt, i, j, k);
    }
    if ((obstacle & OBSTACLE_BACK) == 0) {
      flowField.getFGH().getVectorElement(i, j, k, 2) = computeTurbulentH3D(velocity, viscosity, parameters, parameters.timestep.dt, i, j, k);
    }
  }
}
