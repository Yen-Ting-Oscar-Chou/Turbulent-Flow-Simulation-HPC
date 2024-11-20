#include "StdAfx.hpp"
#include <vector>

#include "VelocityBufferFillStencil.hpp"

Stencils::VelocityBufferFillStencil::VelocityBufferFillStencil(const Parameters& parameters):
  BoundaryStencil(parameters),
  velocityLeft_(selectVelocityLeft(parameters)),
  velocityRight_(selectVelocityRight(parameters)),
  velocityBottom_(selectVelocityBottom(parameters)),
  velocityTop_(selectVelocityTop(parameters)),
  velocityFront_(selectVelocityFront(parameters)),
  velocityBack_(selectVelocityBack(parameters)) {

  // Read the pressure values in each of the six (3D) boundary faces of a subdomain
  // and store them consecutively in one-dimensional buffer arrays.

  if (parameters.geometry.dim == 2) {
    // Initialize the vectors if 2D
    velocityLeft_.resize(parameters.parallel.localSize[1] * 2);
    velocityRight_.resize(parameters.parallel.localSize[1] * 2);
    velocityBottom_.resize(parameters.parallel.localSize[0] * 2);
    velocityTop_.resize(parameters.parallel.localSize[0] * 2);
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

RealType &Stencils::VelocityBufferFillStencil::getVector(FlowField &flowField, int i, int j) {
  return *flowField.getVelocity().getVector(i, j);
}

RealType &Stencils::VelocityBufferFillStencil::getVector(FlowField &flowField, int i, int j, int k) {
  return *flowField.getVelocity().getVector(i, j, k);
}

/* Methods for 2D case */

void Stencils::VelocityBufferFillStencil::applyLeftWall(FlowField& flowField, int i, int j) {
  velocityLeft_[2 * (j - 2)] = flowField.getVelocity().getVector(i, j)[0];
  velocityLeft_[2 * (j - 2) + 1] = flowField.getVelocity().getVector(i, j)[1];
}
void Stencils::VelocityBufferFillStencil::applyRightWall(FlowField& flowField, int i, int j) {
  velocityRight_[2 * (j - 2)] = flowField.getVelocity().getVector(i - 1, j)[0];
  velocityRight_[2 * (j - 2) + 1] = flowField.getVelocity().getVector(i, j)[1];
}

void Stencils::VelocityBufferFillStencil::applyBottomWall(FlowField& flowField, int i, int j) {
  velocityBottom_[2 * (i - 2)] = flowField.getVelocity().getVector(i, j)[0];
  velocityBottom_[2 * (i - 2) + 1] = flowField.getVelocity().getVector(i, j)[1];
}

void Stencils::VelocityBufferFillStencil::applyTopWall(FlowField& flowField, int i, int j) {
  velocityTop_[2 * (i - 2)] = flowField.getVelocity().getVector(i, j)[0];
  velocityTop_[2 * (i - 2) + 1] = flowField.getVelocity().getVector(i, j - 1)[1];
}

/* Methods for 3D case */

void Stencils::VelocityBufferFillStencil::applyLeftWall(FlowField& flowField, int i, int j, int k) {
  velocityLeft_[3 * (j - 2 + flowField.getNy() * (k - 2))] = flowField.getVelocity().getVector(i, j, k)[0];
  velocityLeft_[3 * (j - 2 + flowField.getNy() * (k - 2)) + 1] = flowField.getVelocity().getVector(i, j, k)[1];
  velocityLeft_[3 * (j - 2 + flowField.getNy() * (k - 2)) + 2] = flowField.getVelocity().getVector(i, j, k)[2];
}

void Stencils::VelocityBufferFillStencil::applyRightWall(FlowField& flowField, int i, int j, int k) {
  velocityRight_[3 * (j - 2 + flowField.getNy() * (k - 2))] = flowField.getVelocity().getVector(i - 1, j, k)[0];
  velocityRight_[3 * (j - 2 + flowField.getNy() * (k - 2)) + 1] = flowField.getVelocity().getVector(i, j, k)[1];
  velocityRight_[3 * (j - 2 + flowField.getNy() * (k - 2)) + 2] = flowField.getVelocity().getVector(i, j, k)[2];
}

void Stencils::VelocityBufferFillStencil::applyBottomWall(FlowField& flowField, int i, int j, int k) {
  velocityBottom_[3 * (i - 2 + flowField.getNx() * (k - 2))] = flowField.getVelocity().getVector(i, j, k)[0];
  velocityBottom_[3 * (i - 2 + flowField.getNx() * (k - 2)) + 1] = flowField.getVelocity().getVector(i, j, k)[1];
  velocityBottom_[3 * (i - 2 + flowField.getNx() * (k - 2)) + 2] = flowField.getVelocity().getVector(i, j, k)[2];
}

void Stencils::VelocityBufferFillStencil::applyTopWall(FlowField& flowField, int i, int j, int k) {
  velocityTop_[3 * (i - 2 + flowField.getNx() * (k - 2))] = flowField.getVelocity().getVector(i, j, k)[0];
  velocityTop_[3 * (i - 2 + flowField.getNx() * (k - 2)) + 1] = flowField.getVelocity().getVector(i, j - 1, k)[1];
  velocityTop_[3 * (i - 2 + flowField.getNx() * (k - 2)) + 2] = flowField.getVelocity().getVector(i, j, k)[2];
}

void Stencils::VelocityBufferFillStencil::applyFrontWall(FlowField& flowField, int i, int j, int k) {
  velocityFront_[3 * (i - 2 + flowField.getNx() * (j - 2))] = flowField.getVelocity().getVector(i, j, k)[0];
  velocityFront_[3 * (i - 2 + flowField.getNx() * (j - 2)) + 1] = flowField.getVelocity().getVector(i, j, k)[1];
  velocityFront_[3 * (i - 2 + flowField.getNx() * (j - 2)) + 2] = flowField.getVelocity().getVector(i, j, k)[2];
}

void Stencils::VelocityBufferFillStencil::applyBackWall(FlowField& flowField, int i, int j, int k) {
  velocityFront_[3 * (i - 2 + flowField.getNx() * (j - 2))] = flowField.getVelocity().getVector(i, j, k)[0];
  velocityFront_[3 * (i - 2 + flowField.getNx() * (j - 2)) + 1] = flowField.getVelocity().getVector(i, j, k)[1];
  velocityFront_[3 * (i - 2 + flowField.getNx() * (j - 2)) + 2] = flowField.getVelocity().getVector(i, j, k - 1)[2];
}
