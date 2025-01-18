#include "StdAfx.hpp"

#include "FGHStencil.hpp"

#include "Definitions.hpp"
#include "StencilFunctions.hpp"

void Stencils::FGHStencil::apply(const Parameters& parameters, FlowField& flowField, int i, int j) {
  // Now the localVelocity array should contain lexicographically ordered elements around the given index
  flowField.getFGH().getVectorElement(i, j, 0) = computeF2D(flowField.getVelocity(), parameters, parameters.timestep.dt, i, j);
  flowField.getFGH().getVectorElement(i, j, 1) = computeG2D(flowField.getVelocity(), parameters, parameters.timestep.dt, i, j);
}

void Stencils::FGHStencil::apply(const Parameters& parameters, FlowField& flowField, int i, int j, int k) {
  // Load local velocities into the center layer of the local array
  const int obstacle = flowField.getFlags().getValue(i, j, k);

  if ((obstacle & OBSTACLE_SELF) == 0) { // If the cell is fluid

    if ((obstacle & OBSTACLE_RIGHT) == 0) { // If the right cell is fluid
      flowField.getFGH().getVectorElement(i, j, 0) = computeF3D(flowField.getVelocity(), parameters, parameters.timestep.dt, i, j, k);
    }
    if ((obstacle & OBSTACLE_TOP) == 0) {
      flowField.getFGH().getVectorElement(i, j, 1) = computeG3D(flowField.getVelocity(), parameters, parameters.timestep.dt, i, j, k);
    }
    if ((obstacle & OBSTACLE_BACK) == 0) {
      flowField.getFGH().getVectorElement(i, j, 2) = computeH3D(flowField.getVelocity(), parameters, parameters.timestep.dt, i, j, k);
    }
  }
}
