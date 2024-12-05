#include "StdAfx.hpp"

#include "ViscosityBufferReadStencil.hpp"

Stencils::ViscosityBufferReadStencil::ViscosityBufferReadStencil(const Parameters& parameters): BoundaryStencil(parameters){
  if (parameters.geometry.dim == 2) {
    // Initialize the vectors if 2D
    viscosityLeft_.resize(parameters.parallel.localSize[1] + 1);
    viscosityRight_.resize(parameters.parallel.localSize[1] + 1);
    viscosityBottom_.resize(parameters.parallel.localSize[0] + 1);
    viscosityTop_.resize(parameters.parallel.localSize[0] + 1);
    viscosityFront_.resize(0);
    viscosityBack_.resize(0);
  } else if (parameters.geometry.dim == 3) {
    // Initialize the vectors if 3D
    viscosityLeft_.resize((parameters.parallel.localSize[1] + 1)* (parameters.parallel.localSize[2] + 1));
    viscosityRight_.resize((parameters.parallel.localSize[1] + 1) * (parameters.parallel.localSize[2] + 1));
    viscosityBottom_.resize((parameters.parallel.localSize[0] + 1) * (parameters.parallel.localSize[2] + 1));
    viscosityTop_.resize((parameters.parallel.localSize[0] + 1) * (parameters.parallel.localSize[2] + 1));
    viscosityFront_.resize((parameters.parallel.localSize[0] + 1) * (parameters.parallel.localSize[1] + 1));
    viscosityBack_.resize((parameters.parallel.localSize[0] + 1) * (parameters.parallel.localSize[1] + 1));
  } else {
    throw std::invalid_argument("Unsupported dimensionality: must be 2 or 3.");
  }
}
// Write the viscosity values from a 1-D array to the correct flowField boundaries.
// 2D case
void Stencils::ViscosityBufferReadStencil::applyLeftWall(TurbulentFlowField& flowField, int i, int j) {
  flowField.getViscosity().getScalar(i, j) = viscosityLeft_[j - 1];
}

void Stencils::ViscosityBufferReadStencil::applyRightWall(TurbulentFlowField& flowField, int i, int j) {
  flowField.getViscosity().getScalar(i + 1, j) = viscosityRight_[j - 1];
}

void Stencils::ViscosityBufferReadStencil::applyBottomWall(TurbulentFlowField& flowField, int i, int j) {
  flowField.getViscosity().getScalar(i, j) = viscosityBottom_[i - 1];
  
}

void Stencils::ViscosityBufferReadStencil::applyTopWall(TurbulentFlowField& flowField, int i, int j) {
  flowField.getViscosity().getScalar(i, j + 1) = viscosityTop_[i - 1];
}

// 3D case
void Stencils::ViscosityBufferReadStencil::applyLeftWall(TurbulentFlowField& flowField, int i, int j, int k) {
  flowField.getViscosity().getScalar(i, j, k) = viscosityLeft_[j - 1 + (flowField.getNy() + 1) * (k - 1)];
}
void Stencils::ViscosityBufferReadStencil::applyRightWall(TurbulentFlowField& flowField, int i, int j, int k) {
  flowField.getViscosity().getScalar(i + 1, j, k) = viscosityRight_[j - 1 + (flowField.getNy() + 1) * (k - 1)];
}
void Stencils::ViscosityBufferReadStencil::applyBottomWall(TurbulentFlowField& flowField, int i, int j, int k) {
  flowField.getViscosity().getScalar(i, j, k) = viscosityBottom_[i - 1 + (flowField.getNx() + 1) * (k - 1)];
}
void Stencils::ViscosityBufferReadStencil::applyTopWall(TurbulentFlowField& flowField, int i, int j, int k) {
  flowField.getViscosity().getScalar(i, j + 1, k) = viscosityTop_[i - 1 + (flowField.getNx() + 1) * (k - 1)];
}
void Stencils::ViscosityBufferReadStencil::applyFrontWall(TurbulentFlowField& flowField, int i, int j, int k) {
  flowField.getViscosity().getScalar(i, j, k) = viscosityFront_[i - 1 + (flowField.getNx() + 1) * (j - 1)];
}
void Stencils::ViscosityBufferReadStencil::applyBackWall(TurbulentFlowField& flowField, int i, int j, int k) {
  flowField.getViscosity().getScalar(i, j, k + 1) = viscosityBack_[i - 1 + (flowField.getNx() + 1) * (j - 1)];
}

std::vector<RealType>& Stencils::ViscosityBufferReadStencil::getviscosityLeft(){
  return viscosityLeft_;
}
std::vector<RealType>& Stencils::ViscosityBufferReadStencil::getviscosityRight(){
  return viscosityRight_;
}
std::vector<RealType>& Stencils::ViscosityBufferReadStencil::getviscosityBottom(){
  return viscosityBottom_;
}
std::vector<RealType>& Stencils::ViscosityBufferReadStencil::getviscosityTop(){
  return viscosityTop_;
}
std::vector<RealType>& Stencils::ViscosityBufferReadStencil::getviscosityFront(){
  return viscosityFront_;
}
std::vector<RealType>& Stencils::ViscosityBufferReadStencil::getviscosityBack(){
  return viscosityBack_;
}