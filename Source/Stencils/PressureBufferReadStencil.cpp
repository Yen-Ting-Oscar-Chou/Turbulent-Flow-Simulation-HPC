
#include "StdAfx.hpp"
#include <vector>
#include "PressureBufferReadStencil.hpp"

Stencils::PressureBufferReadStencil::PressureBufferReadStencil(const Parameters& parameters, const std::vector<RealType>& pressure)
  : BoundaryStencil(parameters), pressure_(pressure) {}

// Write the pressure values from a 1-D array to the correct flowField boundaries.

void Stencils::PressureBufferReadStencil::applyLeftWall(FlowField& flowField, int i, int j) {
  flowField.getPressure().getScalar(i, j) = pressure_[j];
}

void Stencils::PressureBufferReadStencil::applyRightWall(FlowField& flowField, int i, int j) {
  flowField.getPressure().getScalar(i, j) = pressure_[j];
}

void Stencils::PressureBufferReadStencil::applyBottomWall(FlowField& flowField, int i, int j) {
  flowField.getPressure().getScalar(i, j) = pressure_[i];
}

void Stencils::PressureBufferReadStencil::applyTopWall(FlowField& flowField, int i, int j) {
  flowField.getPressure().getScalar(i, j) = pressure_[i];
}
void Stencils::PressureBufferReadStencil::applyLeftWall(FlowField& flowField, int i, int j, int k) {
  flowField.getPressure().getScalar(i, j, k) = pressure_[j + BoundaryStencil<FlowField>::parameters_.parallel.localSize[1] * k];
}
void Stencils::PressureBufferReadStencil::applyRightWall(FlowField& flowField, int i, int j, int k) {
  flowField.getPressure().getScalar(i, j, k) = pressure_[j + BoundaryStencil<FlowField>::parameters_.parallel.localSize[1] * k];
}
void Stencils::PressureBufferReadStencil::applyBottomWall(FlowField& flowField, int i, int j, int k) {
  flowField.getPressure().getScalar(i, j, k) = pressure_[i + BoundaryStencil<FlowField>::parameters_.parallel.localSize[0] * k];
}
void Stencils::PressureBufferReadStencil::applyTopWall(FlowField& flowField, int i, int j, int k) {
  flowField.getPressure().getScalar(i, j, k) = pressure_[i + BoundaryStencil<FlowField>::parameters_.parallel.localSize[0] * k];
}
void Stencils::PressureBufferReadStencil::applyFrontWall(FlowField& flowField, int i, int j, int k) {
  flowField.getPressure().getScalar(i, j, k) = pressure_[i + BoundaryStencil<FlowField>::parameters_.parallel.localSize[0] * j];
}
void Stencils::PressureBufferReadStencil::applyBackWall(FlowField& flowField, int i, int j, int k) {
  flowField.getPressure().getScalar(i, j, k) = pressure_[i + BoundaryStencil<FlowField>::parameters_.parallel.localSize[0] * j];
}
