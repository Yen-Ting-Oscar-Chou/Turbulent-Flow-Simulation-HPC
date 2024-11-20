#include "StdAfx.hpp"
#include <vector>

#include "PressureBufferFillStencil.hpp"

Stencils::PressureBufferFillStencil::PressureBufferFillStencil(const Parameters& parameters):
  BoundaryStencil(parameters), pressureLeft_(selectPressureLeft(parameters)),
  pressureRight_(selectPressureRight(parameters)),
  pressureBottom_(selectPressureBottom(parameters)),
  pressureTop_(selectPressureTop(parameters)),
  pressureFront_(selectPressureFront(parameters)),
  pressureBack_(selectPressureBack(parameters)) {
  // Read the pressure values in each of the six (3D) boundary faces of a subdomain
  // and store them consecutively in one-dimensional buffer arrays.

  if (parameters.geometry.dim == 2) {
    // Initialize the vectors if 2D
    pressureLeft_.resize(parameters.parallel.localSize[1]);
    pressureRight_.resize(parameters.parallel.localSize[1]);
    pressureBottom_.resize(parameters.parallel.localSize[0]);
    pressureTop_.resize(parameters.parallel.localSize[0]);
  } else if (parameters.geometry.dim == 3) {
    // Initialize the vectors if 3D
    pressureLeft_.resize(parameters.parallel.localSize[1] * parameters.parallel.localSize[2]);
    pressureRight_.resize(parameters.parallel.localSize[1] * parameters.parallel.localSize[2]);
    pressureBottom_.resize(parameters.parallel.localSize[0] * parameters.parallel.localSize[2]);
    pressureTop_.resize(parameters.parallel.localSize[0] * parameters.parallel.localSize[2]);
    pressureFront_.resize(parameters.parallel.localSize[0] * parameters.parallel.localSize[1]);
    pressureBack_.resize(parameters.parallel.localSize[0] * parameters.parallel.localSize[1]);
  } else {
    throw std::invalid_argument("Unsupported dimensionality: must be 2 or 3.");
  }
}

RealType &Stencils::PressureBufferFillStencil::getScalar(FlowField &flowField, int i, int j) {
  return flowField.getPressure().getScalar(i, j);
}

RealType &Stencils::PressureBufferFillStencil::getScalar(FlowField &flowField, int i, int j, int k) {
  return flowField.getPressure().getScalar(i, j, k);
}

/* Methods for 2D case */
// iterate "- 2" due to ghost cells in each direction. 
void Stencils::PressureBufferFillStencil::applyLeftWall(FlowField& flowField, int i, int j) {
  pressureLeft_[j - 2] = flowField.getPressure().getScalar(i, j);
}

void Stencils::PressureBufferFillStencil::applyRightWall(FlowField& flowField, int i, int j) {
  pressureRight_[j - 2] = flowField.getPressure().getScalar(i, j);
}

void Stencils::PressureBufferFillStencil::applyBottomWall(FlowField& flowField, int i, int j) {
  pressureBottom_[i - 2] = flowField.getPressure().getScalar(i, j);
}

void Stencils::PressureBufferFillStencil::applyTopWall(FlowField& flowField, int i, int j) {
  pressureBottom_[i - 2] = flowField.getPressure().getScalar(i, j);
}

/* Methods for 3D case */

void Stencils::PressureBufferFillStencil::applyLeftWall(FlowField& flowField, int i, int j, int k) {
  pressureLeft_[j - 2 + flowField.getNy() * (k - 2)] = flowField.getPressure().getScalar(i, j, k);
}
void Stencils::PressureBufferFillStencil::applyRightWall(FlowField& flowField, int i, int j, int k) {
  pressureRight_[j - 2 + flowField.getNy() * (k - 2)] = flowField.getPressure().getScalar(i, j, k);
}
void Stencils::PressureBufferFillStencil::applyBottomWall(FlowField& flowField, int i, int j, int k) {
  pressureBottom_[i - 2 + flowField.getNx() * (k - 2)] = flowField.getPressure().getScalar(i, j, k);
}
void Stencils::PressureBufferFillStencil::applyTopWall(FlowField& flowField, int i, int j, int k) {
  pressureTop_[i - 2 + flowField.getNx() * (k - 2)] = flowField.getPressure().getScalar(i, j, k);
}
void Stencils::PressureBufferFillStencil::applyFrontWall(FlowField& flowField, int i, int j, int k) {
  pressureFront_[i - 2 + flowField.getNx() * (j - 2)] = flowField.getPressure().getScalar(i, j, k);
}
void Stencils::PressureBufferFillStencil::applyBackWall(FlowField& flowField, int i, int j, int k) {
  pressureBack_[i - 2 + flowField.getNx() * (j - 2)] = flowField.getPressure().getScalar(i, j, k);
}
