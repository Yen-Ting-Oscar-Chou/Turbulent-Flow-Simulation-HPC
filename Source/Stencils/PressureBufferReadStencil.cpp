#include "StdAfx.hpp"

#include "PressureBufferReadStencil.hpp"

Stencils::PressureBufferReadStencil::PressureBufferReadStencil(const Parameters& parameters): BoundaryStencil(parameters){
  if (parameters.geometry.dim == 2) {
    // Initialize the vectors if 2D
    pressureLeft_.resize(parameters.parallel.localSize[1] + 2);
    pressureRight_.resize(parameters.parallel.localSize[1] + 2);
    pressureBottom_.resize(parameters.parallel.localSize[0] + 2);
    pressureTop_.resize(parameters.parallel.localSize[0] + 2);
    pressureFront_.resize(0);
    pressureBack_.resize(0);
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
// Write the pressure values from a 1-D array to the correct flowField boundaries.
// 2D case
void Stencils::PressureBufferReadStencil::applyLeftWall(FlowField& flowField, int i, int j) {
  flowField.getPressure().getScalar(i, j) = pressureLeft_[j - 1];
}

void Stencils::PressureBufferReadStencil::applyRightWall(FlowField& flowField, int i, int j) {
  flowField.getPressure().getScalar(i, j) = pressureRight_[j - 1];
}

void Stencils::PressureBufferReadStencil::applyBottomWall(FlowField& flowField, int i, int j) {
  flowField.getPressure().getScalar(i, j) = pressureBottom_[i - 1];
  
}

void Stencils::PressureBufferReadStencil::applyTopWall(FlowField& flowField, int i, int j) {
  flowField.getPressure().getScalar(i, j) = pressureTop_[i - 1];
}

// 3D case
void Stencils::PressureBufferReadStencil::applyLeftWall(FlowField& flowField, int i, int j, int k) {
  flowField.getPressure().getScalar(i, j, k) = pressureLeft_[j - 1 + (flowField.getNy() + 2) * (k - 1)];
}
void Stencils::PressureBufferReadStencil::applyRightWall(FlowField& flowField, int i, int j, int k) {
  flowField.getPressure().getScalar(i, j, k) = pressureRight_[j - 1 + (flowField.getNy() + 2) * (k - 1)];
}
void Stencils::PressureBufferReadStencil::applyBottomWall(FlowField& flowField, int i, int j, int k) {
  flowField.getPressure().getScalar(i, j, k) = pressureBottom_[i - 1 + (flowField.getNx() + 2) * (k - 1)];
}
void Stencils::PressureBufferReadStencil::applyTopWall(FlowField& flowField, int i, int j, int k) {
  flowField.getPressure().getScalar(i, j, k) = pressureTop_[i - 1 + (flowField.getNx() + 2) * (k - 1)];
}
void Stencils::PressureBufferReadStencil::applyFrontWall(FlowField& flowField, int i, int j, int k) {
  flowField.getPressure().getScalar(i, j, k) = pressureFront_[i - 1 + (flowField.getNx() + 2) * (j - 1)];
}
void Stencils::PressureBufferReadStencil::applyBackWall(FlowField& flowField, int i, int j, int k) {
  flowField.getPressure().getScalar(i, j, k) = pressureBack_[i - 1 + (flowField.getNx() + 2) * (j - 1)];
}

std::vector<RealType>& Stencils::PressureBufferReadStencil::getpressureLeft(){
  return pressureLeft_;
}
std::vector<RealType>& Stencils::PressureBufferReadStencil::getpressureRight(){
  return pressureRight_;
}
std::vector<RealType>& Stencils::PressureBufferReadStencil::getpressureBottom(){
  return pressureBottom_;
}
std::vector<RealType>& Stencils::PressureBufferReadStencil::getpressureTop(){
  return pressureTop_;
}
std::vector<RealType>& Stencils::PressureBufferReadStencil::getpressureFront(){
  return pressureFront_;
}
std::vector<RealType>& Stencils::PressureBufferReadStencil::getpressureBack(){
  return pressureBack_;
}