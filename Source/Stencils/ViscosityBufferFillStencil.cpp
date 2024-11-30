#include "StdAfx.hpp"

#include "ViscosityBufferFillStencil.hpp"

Stencils::ViscosityBufferFillStencil::ViscosityBufferFillStencil(const Parameters& parameters):
  BoundaryStencil(parameters) {
  // Read the viscosity values in each of the six (3D) boundary faces of a subdomain
  // and store them consecutively in one-dimensional buffer arrays.

  if (parameters.geometry.dim == 2) {
    // Initialize the vectors if 2D
    viscosityLeft_.resize(parameters.parallel.localSize[1]);
    viscosityRight_.resize(parameters.parallel.localSize[1]);
    viscosityBottom_.resize(parameters.parallel.localSize[0]);
    viscosityTop_.resize(parameters.parallel.localSize[0]);
  } else if (parameters.geometry.dim == 3) {
    // Initialize the vectors if 3D
    viscosityLeft_.resize(parameters.parallel.localSize[1] * parameters.parallel.localSize[2]);
    viscosityRight_.resize(parameters.parallel.localSize[1] * parameters.parallel.localSize[2]);
    viscosityBottom_.resize(parameters.parallel.localSize[0] * parameters.parallel.localSize[2]);
    viscosityTop_.resize(parameters.parallel.localSize[0] * parameters.parallel.localSize[2]);
    viscosityFront_.resize(parameters.parallel.localSize[0] * parameters.parallel.localSize[1]);
    viscosityBack_.resize(parameters.parallel.localSize[0] * parameters.parallel.localSize[1]);
  } else {
    throw std::invalid_argument("Unsupported dimensionality: must be 2 or 3.");
  }
}

/* Methods for 2D case */
// iterate "- 2" due to ghost cells in each direction. 
void Stencils::ViscosityBufferFillStencil::applyLeftWall(TurbulentFlowField& flowField, int i, int j) {
  viscosityLeft_[j - 2] = flowField.getViscosity().getScalar(i, j);
}

void Stencils::ViscosityBufferFillStencil::applyRightWall(TurbulentFlowField& flowField, int i, int j) {
  viscosityRight_[j - 2] = flowField.getViscosity().getScalar(i, j);
}


void Stencils::ViscosityBufferFillStencil::applyBottomWall(TurbulentFlowField& flowField, int i, int j) {
  viscosityBottom_[i - 2] = flowField.getViscosity().getScalar(i, j);
}

void Stencils::ViscosityBufferFillStencil::applyTopWall(TurbulentFlowField& flowField, int i, int j) {
  viscosityTop_[i - 2] = flowField.getViscosity().getScalar(i, j);
}

/* Methods for 3D case */

void Stencils::ViscosityBufferFillStencil::applyLeftWall(TurbulentFlowField& flowField, int i, int j, int k) {
  viscosityLeft_[j - 2 + flowField.getNy() * (k - 2)] = flowField.getViscosity().getScalar(i, j, k);
}

void Stencils::ViscosityBufferFillStencil::applyRightWall(TurbulentFlowField& flowField, int i, int j, int k) {
  viscosityRight_[j - 2 + flowField.getNy() * (k - 2)] = flowField.getViscosity().getScalar(i, j, k);
}

void Stencils::ViscosityBufferFillStencil::applyBottomWall(TurbulentFlowField& flowField, int i, int j, int k) {
  viscosityBottom_[i - 2 + flowField.getNx() * (k - 2)] = flowField.getViscosity().getScalar(i, j, k);
}

void Stencils::ViscosityBufferFillStencil::applyTopWall(TurbulentFlowField& flowField, int i, int j, int k) {
  viscosityTop_[i - 2 + flowField.getNx() * (k - 2)] = flowField.getViscosity().getScalar(i, j, k);
}

void Stencils::ViscosityBufferFillStencil::applyFrontWall(TurbulentFlowField& flowField, int i, int j, int k) {
  viscosityFront_[i - 2 + flowField.getNx() * (j - 2)] = flowField.getViscosity().getScalar(i, j, k);
}

void Stencils::ViscosityBufferFillStencil::applyBackWall(TurbulentFlowField& flowField, int i, int j, int k) {
  viscosityBack_[i - 2 + flowField.getNx() * (j - 2)] = flowField.getViscosity().getScalar(i, j, k);
}

std::vector<RealType>& Stencils::ViscosityBufferFillStencil::getviscosityLeft(){
  return viscosityLeft_;
}
std::vector<RealType>& Stencils::ViscosityBufferFillStencil::getviscosityRight(){
  return viscosityRight_;
}
std::vector<RealType>& Stencils::ViscosityBufferFillStencil::getviscosityBottom(){
  return viscosityBottom_;
}
std::vector<RealType>& Stencils::ViscosityBufferFillStencil::getviscosityTop(){
  return viscosityTop_;
}
std::vector<RealType>& Stencils::ViscosityBufferFillStencil::getviscosityFront(){
  return viscosityFront_;
}
std::vector<RealType>& Stencils::ViscosityBufferFillStencil::getviscosityBack(){
  return viscosityBack_;
}