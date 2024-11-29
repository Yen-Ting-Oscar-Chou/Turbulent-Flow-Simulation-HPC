#include "StdAfx.hpp"

#include "VelocityBufferReadStencil.hpp"

Stencils::VelocityBufferReadStencil::VelocityBufferReadStencil(const Parameters& parameters):BoundaryStencil(parameters) {
  if (parameters.geometry.dim == 2) {
    // Initialize the vectors if 2D
    velocityLeft_.resize(parameters.parallel.localSize[1] * 2);
    velocityRight_.resize(parameters.parallel.localSize[1] * 2);
    velocityBottom_.resize(parameters.parallel.localSize[0] * 2);
    velocityTop_.resize(parameters.parallel.localSize[0] * 2);
    velocityFront_.resize(0);
    velocityBack_.resize(0);
  } else if (parameters.geometry.dim == 3) {
    // Initialize the vectors if 3D
    velocityLeft_.resize(parameters.parallel.localSize[1] * parameters.parallel.localSize[2] * 3);
    velocityRight_.resize(parameters.parallel.localSize[1] * parameters.parallel.localSize[2] * 3);
    velocityBottom_.resize(parameters.parallel.localSize[0] * parameters.parallel.localSize[2] * 3);
    velocityTop_.resize(parameters.parallel.localSize[0] * parameters.parallel.localSize[2] * 3);
    velocityFront_.resize(parameters.parallel.localSize[0] * parameters.parallel.localSize[1] * 3);
    velocityBack_.resize(parameters.parallel.localSize[0] * parameters.parallel.localSize[1] * 3);
  } else {
    throw std::invalid_argument("Unsupported dimensionality. Must be 2 or 3.");
  }
}
// Write the velocity values from a 1-D array to the correct flowField boundaries.
// 2D
void Stencils::VelocityBufferReadStencil::applyLeftWall(FlowField& flowField, int i, int j) {
  flowField.getVelocity().getVector(i - 2, j)[0] = velocityLeft_[2 * (j - 2)];
  flowField.getVelocity().getVector(i - 1, j)[1] = velocityLeft_[2 * (j - 2) + 1];
}

void Stencils::VelocityBufferReadStencil::applyRightWall(FlowField& flowField, int i, int j) {
  flowField.getVelocity().getVector(i + 1, j)[0] = velocityRight_[2 * (j - 2)];
  flowField.getVelocity().getVector(i + 1, j)[1] = velocityRight_[2 * (j - 2) + 1];
}

void Stencils::VelocityBufferReadStencil::applyBottomWall(FlowField& flowField, int i, int j) {
  flowField.getVelocity().getVector(i, j - 1)[0] = velocityBottom_[2 * (i - 2)];
  flowField.getVelocity().getVector(i, j - 2)[1] = velocityBottom_[2 * (i - 2) + 1];
}

void Stencils::VelocityBufferReadStencil::applyTopWall(FlowField& flowField, int i, int j) {
  flowField.getVelocity().getVector(i, j + 1)[0] = velocityTop_[2 * (i - 2)];
  flowField.getVelocity().getVector(i, j + 1)[1] = velocityTop_[2 * (i - 2) + 1];
}
// 3D
void Stencils::VelocityBufferReadStencil::applyLeftWall(FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVector(i - 2, j, k)[0] = velocityLeft_[3 * (j - 2 + flowField.getCellsY() * (k - 2))];
  flowField.getVelocity().getVector(i - 1, j, k)[1] = velocityLeft_[3 * (j - 2 + flowField.getCellsY() * (k - 2)) + 1];
  flowField.getVelocity().getVector(i - 1, j, k)[2] = velocityLeft_[3 * (j - 2 + flowField.getCellsY() * (k - 2)) + 2];
}
void Stencils::VelocityBufferReadStencil::applyRightWall(FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVector(i + 1, j, k)[0] = velocityRight_[3 * (j - 2 + flowField.getCellsY() * (k - 2))];
  flowField.getVelocity().getVector(i + 1, j, k)[1] = velocityRight_[3 * (j - 2 + flowField.getCellsY() * (k - 2)) + 1];
  flowField.getVelocity().getVector(i + 1, j, k)[2] = velocityRight_[3 * (j - 2 + flowField.getCellsY() * (k - 2)) + 2];
}
void Stencils::VelocityBufferReadStencil::applyBottomWall(FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVector(i, j - 1, k)[0] = velocityBottom_[3 * (i - 2 + flowField.getCellsX() * (k - 2))];
  flowField.getVelocity().getVector(i, j - 2, k)[1] = velocityBottom_[3 * (i - 2 + flowField.getCellsX() * (k - 2)) + 1];
  flowField.getVelocity().getVector(i, j - 1, k)[2] = velocityBottom_[3 * (i - 2 + flowField.getCellsX() * (k - 2)) + 2];
}
void Stencils::VelocityBufferReadStencil::applyTopWall(FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVector(i, j + 1, k)[0] = velocityTop_[3 * (i - 2 + flowField.getCellsX() * (k - 2))];
  flowField.getVelocity().getVector(i, j + 1, k)[1] = velocityTop_[3 * (i - 2 + flowField.getCellsX() * (k - 2)) + 1];
  flowField.getVelocity().getVector(i, j + 1, k)[2] = velocityTop_[3 * (i - 2 + flowField.getCellsX() * (k - 2)) + 2];
}
void Stencils::VelocityBufferReadStencil::applyFrontWall(FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVector(i, j, k - 1)[0] = velocityFront_[3 * (i - 2 + flowField.getCellsX() * (j - 2))];
  flowField.getVelocity().getVector(i, j, k - 1)[1] = velocityFront_[3 * (i - 2 + flowField.getCellsX() * (j - 2)) + 1];
  flowField.getVelocity().getVector(i, j, k - 2)[2] = velocityFront_[3 * (i - 2 + flowField.getCellsX() * (j - 2)) + 2];
}
void Stencils::VelocityBufferReadStencil::applyBackWall(FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVector(i, j, k + 1)[0] = velocityBack_[3 * (i - 2 + flowField.getCellsX() * (j - 2))];
  flowField.getVelocity().getVector(i, j, k + 1)[1] = velocityBack_[3 * (i - 2 + flowField.getCellsX() * (j - 2)) + 1];
  flowField.getVelocity().getVector(i, j, k + 1)[2] = velocityBack_[3 * (i - 2 + flowField.getCellsX() * (j - 2)) + 2];
}

std::vector<RealType>& Stencils::VelocityBufferReadStencil::getvelocityLeft(){
  return velocityLeft_;
}
std::vector<RealType>& Stencils::VelocityBufferReadStencil::getvelocityRight(){
  return velocityRight_;
}
std::vector<RealType>& Stencils::VelocityBufferReadStencil::getvelocityTop(){
  return velocityTop_;
}
std::vector<RealType>& Stencils::VelocityBufferReadStencil::getvelocityBottom(){
  return velocityBottom_;
}
std::vector<RealType>& Stencils::VelocityBufferReadStencil::getvelocityFront(){
  return velocityFront_;
}
std::vector<RealType>& Stencils::VelocityBufferReadStencil::getvelocityBack(){
  return velocityBack_;
}