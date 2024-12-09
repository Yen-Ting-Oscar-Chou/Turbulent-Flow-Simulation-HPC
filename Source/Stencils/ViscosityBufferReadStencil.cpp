#include "StdAfx.hpp"

#include "ViscosityBufferReadStencil.hpp"

Stencils::ViscosityBufferReadStencil::ViscosityBufferReadStencil(const Parameters& parameters): BufferStencilBase(parameters){
  if (parameters.geometry.dim == 2) {
    // Initialize the vectors if 2D
    bufferLeft_.resize(parameters.parallel.localSize[1] + 2);
    bufferRight_.resize(parameters.parallel.localSize[1] + 2);
    bufferBottom_.resize(parameters.parallel.localSize[0] + 2);
    bufferTop_.resize(parameters.parallel.localSize[0] + 2);
    bufferFront_.resize(0);
    bufferBack_.resize(0);
  } else if (parameters.geometry.dim == 3) {
    // Initialize the vectors if 3D
    bufferLeft_.resize((parameters.parallel.localSize[1] + 2)* (parameters.parallel.localSize[2] + 2));
    bufferRight_.resize((parameters.parallel.localSize[1] + 2) * (parameters.parallel.localSize[2] + 2));
    bufferBottom_.resize((parameters.parallel.localSize[0] + 2) * (parameters.parallel.localSize[2] + 2));
    bufferTop_.resize((parameters.parallel.localSize[0] + 2) * (parameters.parallel.localSize[2] + 2));
    bufferFront_.resize((parameters.parallel.localSize[0] + 2) * (parameters.parallel.localSize[1] + 2));
    bufferBack_.resize((parameters.parallel.localSize[0] + 2) * (parameters.parallel.localSize[1] + 2));
  } else {
    throw std::invalid_argument("Unsupported dimensionality: must be 2 or 3.");
  }
}
// Write the viscosity values from a 1-D array to the correct flowField boundaries.
// 2D case
void Stencils::ViscosityBufferReadStencil::applyLeftWall(TurbulentFlowField& flowField, int i, int j) {
  flowField.getViscosity().getScalar(i, j) = bufferLeft_[j - 1];
}

void Stencils::ViscosityBufferReadStencil::applyRightWall(TurbulentFlowField& flowField, int i, int j) {
  flowField.getViscosity().getScalar(i, j) = bufferRight_[j - 1];
}

void Stencils::ViscosityBufferReadStencil::applyBottomWall(TurbulentFlowField& flowField, int i, int j) {
  flowField.getViscosity().getScalar(i, j) = bufferBottom_[i - 1];
  
}

void Stencils::ViscosityBufferReadStencil::applyTopWall(TurbulentFlowField& flowField, int i, int j) {
  flowField.getViscosity().getScalar(i, j) = bufferTop_[i - 1];
}

// 3D case
void Stencils::ViscosityBufferReadStencil::applyLeftWall(TurbulentFlowField& flowField, int i, int j, int k) {
  flowField.getViscosity().getScalar(i, j, k) = bufferLeft_[j - 1 + (flowField.getNy() + 2) * (k - 1)];
}
void Stencils::ViscosityBufferReadStencil::applyRightWall(TurbulentFlowField& flowField, int i, int j, int k) {
  flowField.getViscosity().getScalar(i, j, k) = bufferRight_[j - 1 + (flowField.getNy() + 2) * (k - 1)];
}
void Stencils::ViscosityBufferReadStencil::applyBottomWall(TurbulentFlowField& flowField, int i, int j, int k) {
  flowField.getViscosity().getScalar(i, j, k) = bufferBottom_[i - 1 + (flowField.getNx() + 2) * (k - 1)];
}
void Stencils::ViscosityBufferReadStencil::applyTopWall(TurbulentFlowField& flowField, int i, int j, int k) {
  flowField.getViscosity().getScalar(i, j, k) = bufferTop_[i - 1 + (flowField.getNx() + 2) * (k - 1)];
}
void Stencils::ViscosityBufferReadStencil::applyFrontWall(TurbulentFlowField& flowField, int i, int j, int k) {
  flowField.getViscosity().getScalar(i, j, k) = bufferFront_[i - 1 + (flowField.getNx() + 2) * (j - 1)];
}
void Stencils::ViscosityBufferReadStencil::applyBackWall(TurbulentFlowField& flowField, int i, int j, int k) {
  flowField.getViscosity().getScalar(i, j, k) = bufferBack_[i - 1 + (flowField.getNx() + 2) * (j - 1)];
}