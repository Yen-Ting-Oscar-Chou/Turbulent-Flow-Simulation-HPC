#include "StdAfx.hpp"

#include "VelocityBufferReadStencil.hpp"

Stencils::VelocityBufferReadStencil::VelocityBufferReadStencil(const Parameters& parameters) {
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
    throw std::invalid_argument("Unsupported dimensionality: must be 2 or 3.");
  }
}
// Write the velocity values from a 1-D array to the correct flowField boundaries.
// 2D
void Stencils::VelocityBufferReadStencil::applyLeftWall([[maybe_unused]] const Parameters& parameters, FlowField& flowField, int i, int j) {
  flowField.getVelocity().getVectorElement(i - 1, j, 0) = bufferLeft_[2 * (j - 1)];
  flowField.getVelocity().getVectorElement(i, j, 1)     = bufferLeft_[2 * (j - 1) + 1];
}

void Stencils::VelocityBufferReadStencil::applyRightWall([[maybe_unused]] const Parameters& parameters, FlowField& flowField, int i, int j) {
  flowField.getVelocity().getVectorElement(i, j, 0) = bufferRight_[2 * (j - 1)];
  flowField.getVelocity().getVectorElement(i, j, 1) = bufferRight_[2 * (j - 1) + 1];
}

void Stencils::VelocityBufferReadStencil::applyBottomWall([[maybe_unused]] const Parameters& parameters, FlowField& flowField, int i, int j) {
  flowField.getVelocity().getVectorElement(i, j, 0)     = bufferBottom_[2 * (i - 1)];
  flowField.getVelocity().getVectorElement(i, j - 1, 1) = bufferBottom_[2 * (i - 1) + 1];
}

void Stencils::VelocityBufferReadStencil::applyTopWall([[maybe_unused]] const Parameters& parameters, FlowField& flowField, int i, int j) {
  flowField.getVelocity().getVectorElement(i, j, 0) = bufferTop_[2 * (i - 1)];
  flowField.getVelocity().getVectorElement(i, j, 1) = bufferTop_[2 * (i - 1) + 1];
}

// 3D
void Stencils::VelocityBufferReadStencil::applyLeftWall([[maybe_unused]] const Parameters& parameters, FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVectorElement(i - 1, j, k, 0) = bufferLeft_[3 * (j - 1 + (flowField.getCellsY() - 1) * (k - 1))];
  flowField.getVelocity().getVectorElement(i, j, k, 1)     = bufferLeft_[3 * (j - 1 + (flowField.getCellsY() - 1) * (k - 1)) + 1];
  flowField.getVelocity().getVectorElement(i, j, k, 2)     = bufferLeft_[3 * (j - 1 + (flowField.getCellsY() - 1) * (k - 1)) + 2];
}
void Stencils::VelocityBufferReadStencil::applyRightWall([[maybe_unused]] const Parameters& parameters, FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVectorElement(i, j, k, 0) = bufferRight_[3 * (j - 1 + (flowField.getCellsY() - 1) * (k - 1))];
  flowField.getVelocity().getVectorElement(i, j, k, 1) = bufferRight_[3 * (j - 1 + (flowField.getCellsY() - 1) * (k - 1)) + 1];
  flowField.getVelocity().getVectorElement(i, j, k, 2) = bufferRight_[3 * (j - 1 + (flowField.getCellsY() - 1) * (k - 1)) + 2];
}
void Stencils::VelocityBufferReadStencil::applyBottomWall([[maybe_unused]] const Parameters& parameters, FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVectorElement(i, j, k, 0)     = bufferBottom_[3 * (i - 1 + (flowField.getCellsX() - 1) * (k - 1))];
  flowField.getVelocity().getVectorElement(i, j - 1, k, 1) = bufferBottom_[3 * (i - 1 + (flowField.getCellsX() - 1) * (k - 1)) + 1];
  flowField.getVelocity().getVectorElement(i, j, k, 2)     = bufferBottom_[3 * (i - 1 + (flowField.getCellsX() - 1) * (k - 1)) + 2];
}
void Stencils::VelocityBufferReadStencil::applyTopWall([[maybe_unused]] const Parameters& parameters, FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVectorElement(i, j, k, 0) = bufferTop_[3 * (i - 1 + (flowField.getCellsX() - 1) * (k - 1))];
  flowField.getVelocity().getVectorElement(i, j, k, 1) = bufferTop_[3 * (i - 1 + (flowField.getCellsX() - 1) * (k - 1)) + 1];
  flowField.getVelocity().getVectorElement(i, j, k, 2) = bufferTop_[3 * (i - 1 + (flowField.getCellsX() - 1) * (k - 1)) + 2];
}
void Stencils::VelocityBufferReadStencil::applyFrontWall([[maybe_unused]] const Parameters& parameters, FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVectorElement(i, j, k, 0)     = bufferFront_[3 * (i - 1 + (flowField.getCellsX() - 1) * (j - 1))];
  flowField.getVelocity().getVectorElement(i, j, k, 1)     = bufferFront_[3 * (i - 1 + (flowField.getCellsX() - 1) * (j - 1)) + 1];
  flowField.getVelocity().getVectorElement(i, j, k - 1, 2) = bufferFront_[3 * (i - 1 + (flowField.getCellsX() - 1) * (j - 1)) + 2];
}
void Stencils::VelocityBufferReadStencil::applyBackWall([[maybe_unused]] const Parameters& parameters, FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVectorElement(i, j, k, 0) = bufferBack_[3 * (i - 1 + (flowField.getCellsX() - 1) * (j - 1))];
  flowField.getVelocity().getVectorElement(i, j, k, 1) = bufferBack_[3 * (i - 1 + (flowField.getCellsX() - 1) * (j - 1)) + 1];
  flowField.getVelocity().getVectorElement(i, j, k, 2) = bufferBack_[3 * (i - 1 + (flowField.getCellsX() - 1) * (j - 1)) + 2];
}