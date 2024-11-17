#include "StdAfx.hpp"

#include "PressureBufferReadStencil.hpp"
#include <vector>

Stencils::PressureBufferFillStencil::PressureBufferFillStencil(const Parameters& parameters, const std::vector<RealType>& pressure)
  : pressure_(pressure), BoundaryStencil(parameters) {}

// Write the pressure values from a 1-D array to the correct flowField boundaries.

void PressureBufferFillStencil::applyLeftWall(FlowField& flowField, int i, int j) {
  flowField.getPressure().getScalar(i, j) = pressure_[j];
}
void PressureBufferFillStencil::applyRightWall(FlowField& flowField, int i, int j) {
  flowField.getPressure().getScalar(i, j) = pressure_[j];
}
void PressureBufferFillStencil::applyBottomWall(FlowField& flowField, int i, int j) {
  flowField.getPressure().getScalar(i, j) = pressure_[i];
}
void PressureBufferFillStencil::applyTopWall(FlowField& flowField, int i, int j) {
  flowField.getPressure().getScalar(i, j) = pressure_[i];
}
void PressureBufferFillStencil::applyLeftWall(FlowField& flowField, int i, int j, int k) {
  flowField.getPressure().getScalar(i, j, k) = pressure_[j + parameters.geometry.localSize[1] * k];
}
void PressureBufferFillStencil::applyRightWall(FlowField& flowField, int i, int j, int k) {
  flowField.getPressure().getScalar(i, j, k) = pressure_[j + parameters.geometry.localSize[1] * k];
}
void PressureBufferFillStencil::applyBottomWall(FlowField& flowField, int i, int j, int k) {
  flowField.getPressure().getScalar(i, j, k) = pressure_[i + parameters.geometry.localSize[0] * k];
}
void PressureBufferFillStencil::applyTopWall(FlowField& flowField, int i, int j, int k) {
  flowField.getPressure().getScalar(i, j, k) = pressure_[i + parameters.geometry.localSize[0] * k];
}
void PressureBufferFillStencil::applyFrontWall(FlowField& flowField, int i, int j, int l) {
  flowField.getPressure().getScalar(i, j, k) = pressure_[i + parameters.geometry.localSize[0] * j];
}
void PressureBufferFillStencil::applyBackWall(FlowField& flowField, int i, int j, int l) {
  flowField.getPressure().getScalar(i, j, k) = pressure_[i + parameters.geometry.localSize[0] * j];
}
