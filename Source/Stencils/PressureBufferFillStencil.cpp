#include "StdAfx.hpp"

#include "PressureBufferFillStencil.hpp"

Stencils::PressureBufferFillStencil::PressureBufferFillStencil(const Parameters& parameters):
  BoundaryStencil(parameters) {
  // Read the pressure values in each of the six (3D) boundary faces of a subdomain
  // and store them consecutively in one-dimensional buffer arrays.

  if (parameters.geometry.dim == 2) {
    // Initialize the vectors if 2D
    pressureLeft_.resize(parameters.parallel.localSize[1] + 2);
    pressureRight_.resize(parameters.parallel.localSize[1] + 2);
    pressureBottom_.resize(parameters.parallel.localSize[0] + 2);
    pressureTop_.resize(parameters.parallel.localSize[0] + 2);
  } else if (parameters.geometry.dim == 3) {
    // Initialize the vectors if 3D
    pressureLeft_.resize((parameters.parallel.localSize[1] + 2)* (parameters.parallel.localSize[2] + 2));
    pressureRight_.resize((parameters.parallel.localSize[1] + 2) * (parameters.parallel.localSize[2] + 2));
    pressureBottom_.resize((parameters.parallel.localSize[0] + 2) * (parameters.parallel.localSize[2] + 2));
    pressureTop_.resize((parameters.parallel.localSize[0] + 2) * (parameters.parallel.localSize[2] + 2));
    pressureFront_.resize((parameters.parallel.localSize[0] + 2) * (parameters.parallel.localSize[1] + 2));
    pressureBack_.resize((parameters.parallel.localSize[0] + 2) * (parameters.parallel.localSize[1] + 2));
  } else {
    throw std::invalid_argument("Unsupported dimensionality: must be 2 or 3.");
  }
}

/* Methods for 2D case */
// iterate "- 1" due to ghost cells in each direction. 
void Stencils::PressureBufferFillStencil::applyLeftWall(FlowField& flowField, int i, int j) {
  pressureLeft_[j - 1] = flowField.getPressure().getScalar(i + 1, j);
}

void Stencils::PressureBufferFillStencil::applyRightWall(FlowField& flowField, int i, int j) {
  pressureRight_[j - 1] = flowField.getPressure().getScalar(i - 1, j);
}


void Stencils::PressureBufferFillStencil::applyBottomWall(FlowField& flowField, int i, int j) {
  pressureBottom_[i - 1] = flowField.getPressure().getScalar(i, j + 1);
}

void Stencils::PressureBufferFillStencil::applyTopWall(FlowField& flowField, int i, int j) {
  pressureTop_[i - 1] = flowField.getPressure().getScalar(i, j - 1);
}

/* Methods for 3D case */

void Stencils::PressureBufferFillStencil::applyLeftWall(FlowField& flowField, int i, int j, int k) {
  pressureLeft_[j - 1 + (flowField.getNy() + 2) * (k - 1)] = flowField.getPressure().getScalar(i + 1, j, k);
}

void Stencils::PressureBufferFillStencil::applyRightWall(FlowField& flowField, int i, int j, int k) {
  pressureRight_[j - 1 + (flowField.getNy() + 2) * (k - 1)] = flowField.getPressure().getScalar(i - 1, j, k);
}

void Stencils::PressureBufferFillStencil::applyBottomWall(FlowField& flowField, int i, int j, int k) {
  pressureBottom_[i - 1 + (flowField.getNx() + 2) * (k - 1)] = flowField.getPressure().getScalar(i, j + 1, k);
}

void Stencils::PressureBufferFillStencil::applyTopWall(FlowField& flowField, int i, int j, int k) {
  pressureTop_[i - 1 + (flowField.getNx() + 2) * (k - 1)] = flowField.getPressure().getScalar(i, j - 1, k);
}

void Stencils::PressureBufferFillStencil::applyFrontWall(FlowField& flowField, int i, int j, int k) {
  pressureFront_[i - 1 + (flowField.getNx() + 2) * (j - 1)] = flowField.getPressure().getScalar(i, j, k + 1);
}

void Stencils::PressureBufferFillStencil::applyBackWall(FlowField& flowField, int i, int j, int k) {
  pressureBack_[i - 1 + (flowField.getNx() + 2) * (j - 1)] = flowField.getPressure().getScalar(i, j, k - 1);
}

std::vector<RealType>& Stencils::PressureBufferFillStencil::getpressureLeft(){
  return pressureLeft_;
}
std::vector<RealType>& Stencils::PressureBufferFillStencil::getpressureRight(){
  return pressureRight_;
}
std::vector<RealType>& Stencils::PressureBufferFillStencil::getpressureBottom(){
  return pressureBottom_;
}
std::vector<RealType>& Stencils::PressureBufferFillStencil::getpressureTop(){
  return pressureTop_;
}
std::vector<RealType>& Stencils::PressureBufferFillStencil::getpressureFront(){
  return pressureFront_;
}
std::vector<RealType>& Stencils::PressureBufferFillStencil::getpressureBack(){
  return pressureBack_;
}