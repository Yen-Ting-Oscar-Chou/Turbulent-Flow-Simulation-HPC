#include "StdAfx.hpp"

#include "TurbulentFGHStencil.hpp"

#include "Definitions.hpp"
#include "StencilFunctions.hpp"

void Stencils::TurbulentFGHStencil::apply(const Parameters& parameters, TurbulentFlowField& flowField, int i, int j) {
  // Load local velocities into the center layer of the local array
  loadLocalVelocity2D(flowField, localVelocity_, i, j);
  loadLocalMeshsize2D(parameters, localMeshsize_, i, j);
  loadLocalViscosity2D(parameters, flowField, localViscosity_, i, j);

  RealType* const values = flowField.getFGH().getVector(i, j);

  // Now the localVelocity array should contain lexicographically ordered elements around the given index
  values[0] = computeTurbulentF2D(localVelocity_, localViscosity_, localMeshsize_, parameters, parameters.timestep.dt);
  values[1] = computeTurbulentG2D(localVelocity_, localViscosity_, localMeshsize_, parameters, parameters.timestep.dt);
}

void Stencils::TurbulentFGHStencil::applyGPU(const Parameters& parameters, VectorField& fgh, ScalarField& viscosity, VectorField& velocity, IntScalarField& flags, int i, int j) {
  // Load local velocities into the center layer of the local array
  loadLocalVelocity2D(velocity, localVelocity_, i, j);
  loadLocalMeshsize2D(parameters, localMeshsize_, i, j);
  loadLocalViscosity2D(parameters, viscosity, localViscosity_, i, j);

  // Now the localVelocity array should contain lexicographically ordered elements around the given index
  fgh.getVectorElement(i, j, 0) = computeTurbulentF2D(localVelocity_, localViscosity_, localMeshsize_, parameters, parameters.timestep.dt);
  fgh.getVectorElement(i, j, 1) = computeTurbulentG2D(localVelocity_, localViscosity_, localMeshsize_, parameters, parameters.timestep.dt);
}

void Stencils::TurbulentFGHStencil::apply(const Parameters& parameters, TurbulentFlowField& flowField, int i, int j, int k) {
  // The same as in 2D, with slight modifications.

  const int       obstacle = flowField.getFlags().getValue(i, j, k);
  RealType* const values   = flowField.getFGH().getVector(i, j, k);

  if ((obstacle & OBSTACLE_SELF) == 0) { // If the cell is fluid
    loadLocalVelocity3D(flowField, localVelocity_, i, j, k);
    loadLocalMeshsize3D(parameters, localMeshsize_, i, j, k);
    loadLocalViscosity3D(parameters, flowField, localViscosity_, i, j, k);

    if ((obstacle & OBSTACLE_RIGHT) == 0) { // If the right cell is fluid
      values[0] = computeTurbulentF3D(localVelocity_, localViscosity_, localMeshsize_, parameters, parameters.timestep.dt);
    }
    if ((obstacle & OBSTACLE_TOP) == 0) {
      values[1] = computeTurbulentG3D(localVelocity_, localViscosity_, localMeshsize_, parameters, parameters.timestep.dt);
    }
    if ((obstacle & OBSTACLE_BACK) == 0) {
      values[2] = computeTurbulentH3D(localVelocity_, localViscosity_, localMeshsize_, parameters, parameters.timestep.dt);
    }
  }
}

void Stencils::TurbulentFGHStencil::applyGPU(const Parameters& parameters, VectorField& fgh, ScalarField& viscosity, VectorField& velocity, IntScalarField& flags, int i, int j, int k) {
  const int       obstacle = flags.getValue(i, j, k);

  if ((obstacle & OBSTACLE_SELF) == 0) { // If the cell is fluid
    loadLocalVelocity3D(velocity, localVelocity_, i, j, k);
    loadLocalMeshsize3D(parameters, localMeshsize_, i, j, k);
    loadLocalViscosity3D(parameters, viscosity, localViscosity_, i, j, k);

    if ((obstacle & OBSTACLE_RIGHT) == 0) { // If the right cell is fluid
      fgh.getVectorElement(i, j, k, 0) = computeTurbulentF3D(localVelocity_, localViscosity_, localMeshsize_, parameters, parameters.timestep.dt);
    }
    if ((obstacle & OBSTACLE_TOP) == 0) {
      fgh.getVectorElement(i, j, k, 1) = computeTurbulentG3D(localVelocity_, localViscosity_, localMeshsize_, parameters, parameters.timestep.dt);
    }
    if ((obstacle & OBSTACLE_BACK) == 0) {
      fgh.getVectorElement(i, j, k, 2) = computeTurbulentH3D(localVelocity_, localViscosity_, localMeshsize_, parameters, parameters.timestep.dt);
    }
  }
}
