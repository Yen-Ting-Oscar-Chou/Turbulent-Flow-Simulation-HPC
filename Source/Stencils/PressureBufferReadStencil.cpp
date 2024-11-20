#include "StdAfx.hpp"
#include <vector>

#include "PressureBufferReadStencil.hpp"

Stencils::PressureBufferReadStencil::PressureBufferReadStencil(const Parameters& parameters, const std::vector<RealType>& pressure)
  : BoundaryStencil(parameters), pressure_(pressure) {}

// Write the pressure values from a 1-D array to the correct flowField boundaries.
// 2D case
void Stencils::PressureBufferReadStencil::applyLeftWall(FlowField& flowField, int i, int j) {
  flowField.getPressure().getScalar(i - 1, j) = pressure_[j - 2];
}

void Stencils::PressureBufferReadStencil::applyRightWall(FlowField& flowField, int i, int j) {
  flowField.getPressure().getScalar(i + 1, j) = pressure_[j - 2];
}

void Stencils::PressureBufferReadStencil::applyBottomWall(FlowField& flowField, int i, int j) {
  flowField.getPressure().getScalar(i, j - 1) = pressure_[i - 2];
}

void Stencils::PressureBufferReadStencil::applyTopWall(FlowField& flowField, int i, int j) {
  flowField.getPressure().getScalar(i, j + 1) = pressure_[i - 2];
}
// 3D case
void Stencils::PressureBufferReadStencil::applyLeftWall(FlowField& flowField, int i, int j, int k) {
  flowField.getPressure().getScalar(i - 1, j, k) = pressure_[j - 2 + flowField.getNy() * (k - 2)];
}
void Stencils::PressureBufferReadStencil::applyRightWall(FlowField& flowField, int i, int j, int k) {
  flowField.getPressure().getScalar(i + 1, j, k) = pressure_[j - 2 + flowField.getNy() * (k - 2)];
}
void Stencils::PressureBufferReadStencil::applyBottomWall(FlowField& flowField, int i, int j, int k) {
  flowField.getPressure().getScalar(i, j - 1, k) = pressure_[i - 2 + flowField.getNx() * (k - 2)];
}
void Stencils::PressureBufferReadStencil::applyTopWall(FlowField& flowField, int i, int j, int k) {
  flowField.getPressure().getScalar(i, j + 1, k) = pressure_[i - 2 + flowField.getNx() * (k - 2)];
}
void Stencils::PressureBufferReadStencil::applyFrontWall(FlowField& flowField, int i, int j, int k) {
  flowField.getPressure().getScalar(i, j, k - 1) = pressure_[i - 2 + flowField.getNx() * (j - 2)];
}
void Stencils::PressureBufferReadStencil::applyBackWall(FlowField& flowField, int i, int j, int k) {
  flowField.getPressure().getScalar(i, j, k + 1) = pressure_[i - 2 + flowField.getNx() * (j - 2)];
}
