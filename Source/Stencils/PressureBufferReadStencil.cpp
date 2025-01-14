#include "StdAfx.hpp"

#include "PressureBufferReadStencil.hpp"

Stencils::PressureBufferReadStencil::PressureBufferReadStencil(const Parameters& parameters) {
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
    bufferLeft_.resize((parameters.parallel.localSize[1] + 2) * (parameters.parallel.localSize[2] + 2));
    bufferRight_.resize((parameters.parallel.localSize[1] + 2) * (parameters.parallel.localSize[2] + 2));
    bufferBottom_.resize((parameters.parallel.localSize[0] + 2) * (parameters.parallel.localSize[2] + 2));
    bufferTop_.resize((parameters.parallel.localSize[0] + 2) * (parameters.parallel.localSize[2] + 2));
    bufferFront_.resize((parameters.parallel.localSize[0] + 2) * (parameters.parallel.localSize[1] + 2));
    bufferBack_.resize((parameters.parallel.localSize[0] + 2) * (parameters.parallel.localSize[1] + 2));
  } else {
    throw std::invalid_argument("Unsupported dimensionality: must be 2 or 3.");
  }
}
// Write the pressure values from a 1-D array to the correct flowField boundaries.
// 2D case
void Stencils::PressureBufferReadStencil::applyLeftWall([[maybe_unused]] const Parameters& parameters, FlowField& flowField, int i, int j) {
  flowField.getPressure().getScalar(i, j) = bufferLeft_[j - 1];
}

void Stencils::PressureBufferReadStencil::applyRightWall([[maybe_unused]] const Parameters& parameters, FlowField& flowField, int i, int j) {
  flowField.getPressure().getScalar(i, j) = bufferRight_[j - 1];
}

void Stencils::PressureBufferReadStencil::applyBottomWall([[maybe_unused]] const Parameters& parameters, FlowField& flowField, int i, int j) {
  flowField.getPressure().getScalar(i, j) = bufferBottom_[i - 1];
}

void Stencils::PressureBufferReadStencil::applyTopWall([[maybe_unused]] const Parameters& parameters, FlowField& flowField, int i, int j) {
  flowField.getPressure().getScalar(i, j) = bufferTop_[i - 1];
}

// 3D case
void Stencils::PressureBufferReadStencil::applyLeftWall([[maybe_unused]] const Parameters& parameters, FlowField& flowField, int i, int j, int k) {
  flowField.getPressure().getScalar(i, j, k) = bufferLeft_[j - 1 + (flowField.getNy() + 2) * (k - 1)];
}
void Stencils::PressureBufferReadStencil::applyRightWall([[maybe_unused]] const Parameters& parameters, FlowField& flowField, int i, int j, int k) {
  flowField.getPressure().getScalar(i, j, k) = bufferRight_[j - 1 + (flowField.getNy() + 2) * (k - 1)];
}
void Stencils::PressureBufferReadStencil::applyBottomWall([[maybe_unused]] const Parameters& parameters, FlowField& flowField, int i, int j, int k) {
  flowField.getPressure().getScalar(i, j, k) = bufferBottom_[i - 1 + (flowField.getNx() + 2) * (k - 1)];
}
void Stencils::PressureBufferReadStencil::applyTopWall([[maybe_unused]] const Parameters& parameters, FlowField& flowField, int i, int j, int k) {
  flowField.getPressure().getScalar(i, j, k) = bufferTop_[i - 1 + (flowField.getNx() + 2) * (k - 1)];
}
void Stencils::PressureBufferReadStencil::applyFrontWall([[maybe_unused]] const Parameters& parameters, FlowField& flowField, int i, int j, int k) {
  flowField.getPressure().getScalar(i, j, k) = bufferFront_[i - 1 + (flowField.getNx() + 2) * (j - 1)];
}
void Stencils::PressureBufferReadStencil::applyBackWall([[maybe_unused]] const Parameters& parameters, FlowField& flowField, int i, int j, int k) {
  flowField.getPressure().getScalar(i, j, k) = bufferBack_[i - 1 + (flowField.getNx() + 2) * (j - 1)];
}
