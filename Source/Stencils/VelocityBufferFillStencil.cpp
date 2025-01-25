#include "StdAfx.hpp"

#include "VelocityBufferFillStencil.hpp"

Stencils::VelocityBufferFillStencil::VelocityBufferFillStencil(const Parameters& parameters):
  BufferStencilBase(parameters) {

  // Read the pressure values in each of the six (3D) boundary faces of a subdomain
  // and store them consecutively in one-dimensional buffer arrays.

  if (parameters.geometry.dim == 2) {
    // Initialize the vectors if 2D
    bufferLeft_.resize((parameters.parallel.localSize[1] + 2) * 2);
    bufferRight_.resize((parameters.parallel.localSize[1] + 2) * 2);
    bufferBottom_.resize((parameters.parallel.localSize[0] + 2) * 2);
    bufferTop_.resize((parameters.parallel.localSize[0] + 2) * 2);
    bufferFront_.resize(0);
    bufferBack_.resize(0);
  } else if (parameters.geometry.dim == 3) {
    // Initialize the vectors if 3D
    bufferLeft_.resize((parameters.parallel.localSize[1] + 2) * (parameters.parallel.localSize[2] + 2) * 3);
    bufferRight_.resize((parameters.parallel.localSize[1] + 2) * (parameters.parallel.localSize[2] + 2) * 3);
    bufferBottom_.resize((parameters.parallel.localSize[0] + 2) * (parameters.parallel.localSize[2] + 2) * 3);
    bufferTop_.resize((parameters.parallel.localSize[0] + 2) * (parameters.parallel.localSize[2] + 2) * 3);
    bufferFront_.resize((parameters.parallel.localSize[0] + 2) * (parameters.parallel.localSize[1] + 2) * 3);
    bufferBack_.resize((parameters.parallel.localSize[0] + 2) * (parameters.parallel.localSize[1] + 2) * 3);
  } else {
    throw std::invalid_argument("Unsupported dimensionality. Must be 2 or 3.");
  }
}

/* Methods for 2D case */

void Stencils::VelocityBufferFillStencil::applyLeftWall(FlowField& flowField, int i, int j) {
  bufferLeft_[2 * (j - 1)] = flowField.getVelocity().getVector(i + 1, j)[0];
  bufferLeft_[2 * (j - 1) + 1] = flowField.getVelocity().getVector(i + 1, j)[1];
}
void Stencils::VelocityBufferFillStencil::applyRightWall(FlowField& flowField, int i, int j) {
  bufferRight_[2 * (j - 1)] = flowField.getVelocity().getVector(i - 2, j)[0];
  bufferRight_[2 * (j - 1) + 1] = flowField.getVelocity().getVector(i - 1, j)[1];
}

void Stencils::VelocityBufferFillStencil::applyBottomWall(FlowField& flowField, int i, int j) {
  bufferBottom_[2 * (i - 1)] = flowField.getVelocity().getVector(i, j + 1)[0];
  bufferBottom_[2 * (i - 1) + 1] = flowField.getVelocity().getVector(i, j + 1)[1];
}

void Stencils::VelocityBufferFillStencil::applyTopWall(FlowField& flowField, int i, int j) {
  bufferTop_[2 * (i - 1)] = flowField.getVelocity().getVector(i, j - 1)[0];
  bufferTop_[2 * (i - 1) + 1] = flowField.getVelocity().getVector(i, j - 2)[1];
}

/* Methods for 3D case */

void Stencils::VelocityBufferFillStencil::applyLeftWall(FlowField& flowField, int i, int j, int k) {
  bufferLeft_[3 * (j - 1 + (flowField.getCellsY() - 1) * (k - 1))] = flowField.getVelocity().getVector(i + 1, j, k)[0];
  bufferLeft_[3 * (j - 1 + (flowField.getCellsY() - 1) * (k - 1)) + 1] = flowField.getVelocity().getVector(i + 1, j, k)[1];
  bufferLeft_[3 * (j - 1 + (flowField.getCellsY() - 1) * (k - 1)) + 2] = flowField.getVelocity().getVector(i + 1, j, k)[2];
}

void Stencils::VelocityBufferFillStencil::applyRightWall(FlowField& flowField, int i, int j, int k) {
  bufferRight_[3 * (j - 1 + (flowField.getCellsY() - 1) * (k - 1))] = flowField.getVelocity().getVector(i - 2, j, k)[0];
  bufferRight_[3 * (j - 1 + (flowField.getCellsY() - 1) * (k - 1)) + 1] = flowField.getVelocity().getVector(i - 1, j, k)[1];
  bufferRight_[3 * (j - 1 + (flowField.getCellsY() - 1) * (k - 1)) + 2] = flowField.getVelocity().getVector(i - 1, j, k)[2];
}

void Stencils::VelocityBufferFillStencil::applyBottomWall(FlowField& flowField, int i, int j, int k) {
  bufferBottom_[3 * (i - 1 + (flowField.getCellsX() - 1) * (k - 1))] = flowField.getVelocity().getVector(i, j + 1, k)[0];
  bufferBottom_[3 * (i - 1 + (flowField.getCellsX() - 1) * (k - 1)) + 1] = flowField.getVelocity().getVector(i, j + 1, k)[1];
  bufferBottom_[3 * (i - 1 + (flowField.getCellsX() - 1) * (k - 1)) + 2] = flowField.getVelocity().getVector(i, j + 1, k)[2];
}

void Stencils::VelocityBufferFillStencil::applyTopWall(FlowField& flowField, int i, int j, int k) {
  bufferTop_[3 * (i - 1 + (flowField.getCellsX() - 1) * (k - 1))] = flowField.getVelocity().getVector(i, j - 1, k)[0];
  bufferTop_[3 * (i - 1 + (flowField.getCellsX() - 1) * (k - 1)) + 1] = flowField.getVelocity().getVector(i, j - 2, k)[1];
  bufferTop_[3 * (i - 1 + (flowField.getCellsX() - 1) * (k - 1)) + 2] = flowField.getVelocity().getVector(i, j - 1, k)[2];
}

void Stencils::VelocityBufferFillStencil::applyFrontWall(FlowField& flowField, int i, int j, int k) {
  bufferFront_[3 * (i - 1 + (flowField.getCellsX() - 1) * (j - 1))] = flowField.getVelocity().getVector(i, j, k + 1)[0];
  bufferFront_[3 * (i - 1 + (flowField.getCellsX() - 1) * (j - 1)) + 1] = flowField.getVelocity().getVector(i, j, k + 1)[1];
  bufferFront_[3 * (i - 1 + (flowField.getCellsX() - 1) * (j - 1)) + 2] = flowField.getVelocity().getVector(i, j, k + 1)[2];
}

void Stencils::VelocityBufferFillStencil::applyBackWall(FlowField& flowField, int i, int j, int k) {
  bufferBack_[3 * (i - 1 + (flowField.getCellsX() - 1) * (j - 1))] = flowField.getVelocity().getVector(i, j, k - 1)[0];
  bufferBack_[3 * (i - 1 + (flowField.getCellsX() - 1) * (j - 1)) + 1] = flowField.getVelocity().getVector(i, j, k - 1)[1];
  bufferBack_[3 * (i - 1 + (flowField.getCellsX() - 1) * (j - 1)) + 2] = flowField.getVelocity().getVector(i, j, k - 2)[2];
}
