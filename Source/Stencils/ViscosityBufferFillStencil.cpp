#include "StdAfx.hpp"

#include "ViscosityBufferFillStencil.hpp"

Stencils::ViscosityBufferFillStencil::ViscosityBufferFillStencil(const Parameters& parameters):
  BufferStencilBase(parameters) {
  // Read the viscosity values in each of the six (3D) boundary faces of a subdomain
  // and store them consecutively in one-dimensional buffer arrays.

  if (parameters.geometry.dim == 2) {
    // Initialize the vectors if 2D
    bufferLeft_.resize(parameters.parallel.localSize[1] + 2);
    bufferRight_.resize(parameters.parallel.localSize[1] + 2);
    bufferRight_.resize(parameters.parallel.localSize[0] + 2);
    bufferTop_.resize(parameters.parallel.localSize[0] + 2);
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

/* Methods for 2D case */
// iterate "- 1" due to ghost cells in each direction. 
void Stencils::ViscosityBufferFillStencil::applyLeftWall(TurbulentFlowField& flowField, int i, int j) {
  bufferLeft_[j - 1] = flowField.getViscosity().getScalar(i + 1, j);
}

void Stencils::ViscosityBufferFillStencil::applyRightWall(TurbulentFlowField& flowField, int i, int j) {
  bufferRight_[j - 1] = flowField.getViscosity().getScalar(i - 1, j);
}


void Stencils::ViscosityBufferFillStencil::applyBottomWall(TurbulentFlowField& flowField, int i, int j) {
  bufferBottom_[i - 1] = flowField.getViscosity().getScalar(i, j + 1);
}

void Stencils::ViscosityBufferFillStencil::applyTopWall(TurbulentFlowField& flowField, int i, int j) {
  bufferTop_[i - 1] = flowField.getViscosity().getScalar(i, j - 1);
}

/* Methods for 3D case */

void Stencils::ViscosityBufferFillStencil::applyLeftWall(TurbulentFlowField& flowField, int i, int j, int k) {
  bufferLeft_[j - 1 + (flowField.getNy() + 2) * (k - 1)] = flowField.getViscosity().getScalar(i + 1, j, k);
}

void Stencils::ViscosityBufferFillStencil::applyRightWall(TurbulentFlowField& flowField, int i, int j, int k) {
  bufferRight_[j - 1 + (flowField.getNy() + 2) * (k - 1)] = flowField.getViscosity().getScalar(i - 1, j, k);
}

void Stencils::ViscosityBufferFillStencil::applyBottomWall(TurbulentFlowField& flowField, int i, int j, int k) {
  bufferBottom_[i - 1 + (flowField.getNx() + 2) * (k - 1)] = flowField.getViscosity().getScalar(i, j + 1, k);
}

void Stencils::ViscosityBufferFillStencil::applyTopWall(TurbulentFlowField& flowField, int i, int j, int k) {
  bufferTop_[i - 1 + (flowField.getNx() + 2) * (k - 1)] = flowField.getViscosity().getScalar(i, j - 1, k);
}

void Stencils::ViscosityBufferFillStencil::applyFrontWall(TurbulentFlowField& flowField, int i, int j, int k) {
  bufferFront_[i - 1 + (flowField.getNx() + 2) * (j - 1)] = flowField.getViscosity().getScalar(i, j, k + 1);
}

void Stencils::ViscosityBufferFillStencil::applyBackWall(TurbulentFlowField& flowField, int i, int j, int k) {
  bufferBack_[i - 1 + (flowField.getNx() + 2) * (j - 1)] = flowField.getViscosity().getScalar(i, j, k - 1);
}