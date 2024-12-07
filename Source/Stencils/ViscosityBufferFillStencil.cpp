#include "StdAfx.hpp"

#include "ViscosityBufferFillStencil.hpp"

Stencils::ViscosityBufferFillStencil::ViscosityBufferFillStencil(const Parameters& parameters):
  BoundaryStencil(parameters) {
  // Read the viscosity values in each of the six (3D) boundary faces of a subdomain
  // and store them consecutively in one-dimensional buffer arrays.

  if (parameters.geometry.dim == 2) {
    // Initialize the vectors if 2D
    viscosityLeft_.resize(parameters.parallel.localSize[1] + 2);
    viscosityRight_.resize(parameters.parallel.localSize[1] + 2);
    viscosityBottom_.resize(parameters.parallel.localSize[0] + 2);
    viscosityTop_.resize(parameters.parallel.localSize[0] + 2);
  } else if (parameters.geometry.dim == 3) {
    // Initialize the vectors if 3D
    viscosityLeft_.resize((parameters.parallel.localSize[1] + 2)* (parameters.parallel.localSize[2] + 2));
    viscosityRight_.resize((parameters.parallel.localSize[1] + 2) * (parameters.parallel.localSize[2] + 2));
    viscosityBottom_.resize((parameters.parallel.localSize[0] + 2) * (parameters.parallel.localSize[2] + 2));
    viscosityTop_.resize((parameters.parallel.localSize[0] + 2) * (parameters.parallel.localSize[2] + 2));
    viscosityFront_.resize((parameters.parallel.localSize[0] + 2) * (parameters.parallel.localSize[1] + 2));
    viscosityBack_.resize((parameters.parallel.localSize[0] + 2) * (parameters.parallel.localSize[1] + 2));
  } else {
    throw std::invalid_argument("Unsupported dimensionality: must be 2 or 3.");
  }
}

/* Methods for 2D case */
// iterate "- 1" due to ghost cells in each direction. 
void Stencils::ViscosityBufferFillStencil::applyLeftWall(TurbulentFlowField& flowField, int i, int j) {
  viscosityLeft_[j - 1] = flowField.getViscosity().getScalar(i + 1, j);
}

void Stencils::ViscosityBufferFillStencil::applyRightWall(TurbulentFlowField& flowField, int i, int j) {
  viscosityRight_[j - 1] = flowField.getViscosity().getScalar(i - 1, j);
}


void Stencils::ViscosityBufferFillStencil::applyBottomWall(TurbulentFlowField& flowField, int i, int j) {
  viscosityBottom_[i - 1] = flowField.getViscosity().getScalar(i, j + 1);
}

void Stencils::ViscosityBufferFillStencil::applyTopWall(TurbulentFlowField& flowField, int i, int j) {
  viscosityTop_[i - 1] = flowField.getViscosity().getScalar(i, j - 1);
}

/* Methods for 3D case */

void Stencils::ViscosityBufferFillStencil::applyLeftWall(TurbulentFlowField& flowField, int i, int j, int k) {
  viscosityLeft_[j - 1 + (flowField.getNy() + 2) * (k - 1)] = flowField.getViscosity().getScalar(i + 1, j, k);
}

void Stencils::ViscosityBufferFillStencil::applyRightWall(TurbulentFlowField& flowField, int i, int j, int k) {
  viscosityRight_[j - 1 + (flowField.getNy() + 2) * (k - 1)] = flowField.getViscosity().getScalar(i - 1, j, k);
}

void Stencils::ViscosityBufferFillStencil::applyBottomWall(TurbulentFlowField& flowField, int i, int j, int k) {
  viscosityBottom_[i - 1 + (flowField.getNx() + 2) * (k - 1)] = flowField.getViscosity().getScalar(i, j + 1, k);
}

void Stencils::ViscosityBufferFillStencil::applyTopWall(TurbulentFlowField& flowField, int i, int j, int k) {
  viscosityTop_[i - 1 + (flowField.getNx() + 2) * (k - 1)] = flowField.getViscosity().getScalar(i, j - 1, k);
}

void Stencils::ViscosityBufferFillStencil::applyFrontWall(TurbulentFlowField& flowField, int i, int j, int k) {
  viscosityFront_[i - 1 + (flowField.getNx() + 2) * (j - 1)] = flowField.getViscosity().getScalar(i, j, k + 1);
}

void Stencils::ViscosityBufferFillStencil::applyBackWall(TurbulentFlowField& flowField, int i, int j, int k) {
  viscosityBack_[i - 1 + (flowField.getNx() + 2) * (j - 1)] = flowField.getViscosity().getScalar(i, j, k - 1);
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