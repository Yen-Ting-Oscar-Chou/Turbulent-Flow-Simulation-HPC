#include "VelocityBufferReadStencil.hpp"

Stencils::VelocityBufferReadStencil::VelocityBufferReadStencil(const Parameters& parameters, const std::vector<RealType>& velocity)
  : BoundaryStencil(parameters), velocity_(velocity) {}

// Write the pressure values from a 1-D array to the correct flowField boundaries.

void Stencils::VelocityBufferReadStencil::applyLeftWall(FlowField& flowField, int i, int j) {
  flowField.getPressure().getScalar(i, j) = velocity_[j];
}

void Stencils::VelocityBufferReadStencil::applyRightWall(FlowField& flowField, int i, int j) {
  flowField.getPressure().getScalar(i, j) = velocity_[j];
}

void Stencils::VelocityBufferReadStencil::applyBottomWall(FlowField& flowField, int i, int j) {
  flowField.getPressure().getScalar(i, j) = velocity_[i];
}

void Stencils::VelocityBufferReadStencil::applyTopWall(FlowField& flowField, int i, int j) {
  flowField.getPressure().getScalar(i, j) = velocity_[i];
}
void Stencils::VelocityBufferReadStencil::applyLeftWall(FlowField& flowField, int i, int j, int k) {
  flowField.getPressure().getScalar(i, j, k) = velocity_[j + BoundaryStencil<FlowField>::parameters_.parallel.localSize[1] * k];
}
void Stencils::VelocityBufferReadStencil::applyRightWall(FlowField& flowField, int i, int j, int k) {
  flowField.getPressure().getScalar(i, j, k) = velocity_[j + BoundaryStencil<FlowField>::parameters_.parallel.localSize[1] * k];
}
void Stencils::VelocityBufferReadStencil::applyBottomWall(FlowField& flowField, int i, int j, int k) {
  flowField.getPressure().getScalar(i, j, k) = velocity_[i + BoundaryStencil<FlowField>::parameters_.parallel.localSize[0] * k];
}
void Stencils::VelocityBufferReadStencil::applyTopWall(FlowField& flowField, int i, int j, int k) {
  flowField.getPressure().getScalar(i, j, k) = velocity_[i + BoundaryStencil<FlowField>::parameters_.parallel.localSize[0] * k];
}
void Stencils::VelocityBufferReadStencil::applyFrontWall(FlowField& flowField, int i, int j, int k) {
  flowField.getPressure().getScalar(i, j, k) = velocity_[i + BoundaryStencil<FlowField>::parameters_.parallel.localSize[0] * j];
}
void Stencils::VelocityBufferReadStencil::applyBackWall(FlowField& flowField, int i, int j, int k) {
  flowField.getPressure().getScalar(i, j, k) = velocity_[i + BoundaryStencil<FlowField>::parameters_.parallel.localSize[0] * j];
}
