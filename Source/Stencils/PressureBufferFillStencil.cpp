#include "StdAfx.hpp"

#include "PressureBufferFillStencil.hpp"
#include <vector>


Stencils::PressureBufferFillStencil::PressureBufferFillStencil(const Parameters& parameters): 
  BoundaryStencil(parameters) {}
// Read the pressure values in each of the six (3D) boundary faces of a sub-domain
// and store them consecutively in one-dimensional buffer arrays.
if(parameters.geometry.dim == 2) {
  pressureLeft_ = std::vector<RealType>(parameters.parallel.localSize[1]);
  pressureRight_ = std::vector<RealType>(parameters.parallel.localSize[1]);
  pressureBottom_ = std::vector<RealType>(parameters.parallel.localSize[0]);
  pressureTop_ = std::vector<RealType>(parameters.parallel.localSize[0]);
}
else if(parameters.geometry.dim == 3) {
  pressureLeft_ = std::vector<RealType>(parameters.parallel.localSize[1] * parameters.parallel.localSize[2]);
  pressureRight_ = std::vector<RealType>(parameters.parallel.localSize[1] * parameters.parallel.localSize[2]);
  pressureBottom_ = std::vector<RealType>(parameters.parallel.localSize[0] * parameters.parallel.localSize[2]);
  pressureTop_ = std::vector<RealType>(parameters.parallel.localSize[0] * parameters.parallel.localSize[2]);
  pressureFront_ = std::vector<RealType>(parameters.parallel.localSize[0] * parameters.parallel.localSize[1]);
  pressureBack_ = std::vector<RealType>(parameters.parallel.localSize[0] * parameters.parallel.localSize[1]);
}

void PressureBufferFillStencil::applyLeftWall(FlowField& flowField, int i, int j) {
  pressureLeft_[j] = flowField.getPressure().getScalar(i, j);
}
void PressureBufferFillStencil::applyRightWall(FlowField& flowField, int i, int j) {
  pressureRight_[j] = flowField.getPressure().getScalar(i, j);
}
void PressureBufferFillStencil::applyBottomWall(FlowField& flowField, int i, int j) {
  pressureBottom_[i] = flowField.getPressure().getScalar(i, j);
}
void PressureBufferFillStencil::applyTopWall(FlowField& flowField, int i, int j) {
  pressureTop_[i] = flowField.getPressure().getScalar(i, j);
}
void PressureBufferFillStencil::applyLeftWall(FlowField& flowField, int i, int j, int k) {
  pressureLeft_[j + parameters.geometry.localSize[1] * k] = flowField.getPressure().getScalar(i, j, k);
}
void PressureBufferFillStencil::applyRightWall(FlowField& flowField, int i, int j, int k) {
  pressureRight_[j + parameters.geometry.localSize[1] * k] = flowField.getPressure().getScalar(i, j, k);
}
void PressureBufferFillStencil::applyBottomWall(FlowField& flowField, int i, int j, int k) {
  pressureBottom_[i + parameters.geometry.localSize[0] * k] = flowField.getPressure().getScalar(i, j, k);
}
void PressureBufferFillStencil::applyTopWall(FlowField& flowField, int i, int j, int k) {
  pressureTop_[i + parameters.geometry.localSize[0] * k] = flowField.getPressure().getScalar(i, j, k);
}
void PressureBufferFillStencil::applyFrontWall(FlowField& flowField, int i, int j, int l) {
  pressureFront_[i + parameters.geometry.localSize[0] * j] = flowField.getPressure().getScalar(i, j, k);
}
void PressureBufferFillStencil::applyBackWall(FlowField& flowField, int i, int j, int l) {
  pressureBack_[i + parameters.geometry.localSize[0] * j] = flowField.getPressure().getScalar(i, j, k);
}
