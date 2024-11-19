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
    velocityLeft_.resize(parameters.parallel.localSize[1]);
    velocityRight_.resize(parameters.parallel.localSize[1]);
    velocityBottom_.resize(parameters.parallel.localSize[0]);
    velocityTop_.resize(parameters.parallel.localSize[0]);
  } else if (parameters.geometry.dim == 3) {
    // Initialize the vectors if 3D
    velocityLeft_.resize(parameters.parallel.localSize[1] * parameters.parallel.localSize[2]);
    velocityRight_.resize(parameters.parallel.localSize[1] * parameters.parallel.localSize[2]);
    velocityBottom_.resize(parameters.parallel.localSize[0] * parameters.parallel.localSize[2]);
    velocityTop_.resize(parameters.parallel.localSize[0] * parameters.parallel.localSize[2]);
    velocityFront_.resize(parameters.parallel.localSize[0] * parameters.parallel.localSize[1]);
    velocityBack_.resize(parameters.parallel.localSize[0] * parameters.parallel.localSize[1]);
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
  velocityLeft_[j] = flowField.getPressure().getScalar(i, j);
}

void Stencils::VelocityBufferFillStencil::applyRightWall(FlowField& flowField, int i, int j) {
  velocityRight_[j] = flowField.getPressure().getScalar(i, j);
}

void Stencils::VelocityBufferFillStencil::applyBottomWall(FlowField& flowField, int i, int j) {
  velocityBottom_[j] = flowField.getPressure().getScalar(i, j);
}

void Stencils::VelocityBufferFillStencil::applyTopWall(FlowField& flowField, int i, int j) {
  velocityBottom_[j] = flowField.getPressure().getScalar(i, j);
}

/* Methods for 3D case */

void Stencils::VelocityBufferFillStencil::applyLeftWall(FlowField& flowField, int i, int j, int k) {
  velocityLeft_[j + BoundaryStencil<FlowField>::parameters_.parallel.localSize[1] * k] = flowField.getPressure().getScalar(i, j, k);
}

void Stencils::VelocityBufferFillStencil::applyRightWall(FlowField& flowField, int i, int j, int k) {
  velocityRight_[j + BoundaryStencil<FlowField>::parameters_.parallel.localSize[1] * k] = flowField.getPressure().getScalar(i, j, k);
}

void Stencils::VelocityBufferFillStencil::applyBottomWall(FlowField& flowField, int i, int j, int k) {
  velocityBottom_[i + BoundaryStencil<FlowField>::parameters_.parallel.localSize[0] * k] = flowField.getPressure().getScalar(i, j, k);
}

void Stencils::VelocityBufferFillStencil::applyTopWall(FlowField& flowField, int i, int j, int k) {
  velocityTop_[i + BoundaryStencil<FlowField>::parameters_.parallel.localSize[0] * k] = flowField.getPressure().getScalar(i, j, k);
}

void Stencils::VelocityBufferFillStencil::applyFrontWall(FlowField& flowField, int i, int j, int k) {
  velocityFront_[i + BoundaryStencil<FlowField>::parameters_.parallel.localSize[0] * j] = flowField.getPressure().getScalar(i, j, k);
}

void Stencils::VelocityBufferFillStencil::applyBackWall(FlowField& flowField, int i, int j, int k) {
  velocityBack_[i +BoundaryStencil<FlowField>::parameters_.parallel.localSize[0] * j] = flowField.getPressure().getScalar(i, j, k);
}
