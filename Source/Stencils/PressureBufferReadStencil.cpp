#include "StdAfx.hpp"
#include <vector>

#include "PressureBufferReadStencil.hpp"

Stencils::PressureBufferReadStencil::PressureBufferReadStencil(const Parameters& parameters): BoundaryStencil(parameters){}
// Write the pressure values from a 1-D array to the correct flowField boundaries.
// 2D case
void Stencils::PressureBufferReadStencil::applyLeftWall(FlowField& flowField, int i, int j) {
  flowField.getPressure().getScalar(i - 1, j) = pressureLeft_[j - 2];
}

void Stencils::PressureBufferReadStencil::applyRightWall(FlowField& flowField, int i, int j) {
  flowField.getPressure().getScalar(i + 1, j) = pressureRight_[j - 2];
}

void Stencils::PressureBufferReadStencil::applyBottomWall(FlowField& flowField, int i, int j) {
  flowField.getPressure().getScalar(i, j - 1) = pressureBottom_[i - 2];
}

void Stencils::PressureBufferReadStencil::applyTopWall(FlowField& flowField, int i, int j) {
  flowField.getPressure().getScalar(i, j + 1) = pressureTop_[i - 2];
}

// 3D case
void Stencils::PressureBufferReadStencil::applyLeftWall(FlowField& flowField, int i, int j, int k) {
  flowField.getPressure().getScalar(i - 1, j, k) = pressureLeft_[j - 2 + flowField.getCellsY() * (k - 2)];
}
void Stencils::PressureBufferReadStencil::applyRightWall(FlowField& flowField, int i, int j, int k) {
  flowField.getPressure().getScalar(i + 1, j, k) = pressureRight_[j - 2 + flowField.getCellsY() * (k - 2)];
}
void Stencils::PressureBufferReadStencil::applyBottomWall(FlowField& flowField, int i, int j, int k) {
  flowField.getPressure().getScalar(i, j - 1, k) = pressureBottom_[i - 2 + flowField.getCellsX() * (k - 2)];
}
void Stencils::PressureBufferReadStencil::applyTopWall(FlowField& flowField, int i, int j, int k) {
  flowField.getPressure().getScalar(i, j + 1, k) = pressureTop_[i - 2 + flowField.getCellsX() * (k - 2)];
}
void Stencils::PressureBufferReadStencil::applyFrontWall(FlowField& flowField, int i, int j, int k) {
  flowField.getPressure().getScalar(i, j, k - 1) = pressureFront_[i - 2 + flowField.getCellsX() * (j - 2)];
}
void Stencils::PressureBufferReadStencil::applyBackWall(FlowField& flowField, int i, int j, int k) {
  flowField.getPressure().getScalar(i, j, k + 1) = pressureBack_[i - 2 + flowField.getCellsX() * (j - 2)];
}

std::vector<RealType>& Stencils::PressureBufferReadStencil::getpressureLeft(){
  return pressureLeft_;
}
std::vector<RealType>& Stencils::PressureBufferReadStencil::getpressureRight(){
  return pressureRight_;
}
std::vector<RealType>& Stencils::PressureBufferReadStencil::getpressureBottom(){
  return pressureBottom_;
}
std::vector<RealType>& Stencils::PressureBufferReadStencil::getpressureTop(){
  return pressureTop_;
}
std::vector<RealType>& Stencils::PressureBufferReadStencil::getpressureFront(){
  return pressureFront_;
}
std::vector<RealType>& Stencils::PressureBufferReadStencil::getpressureBack(){
  return pressureBack_;
}