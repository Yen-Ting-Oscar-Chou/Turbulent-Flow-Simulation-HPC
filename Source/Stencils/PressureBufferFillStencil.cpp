#include "StdAfx.hpp"

#include "PressureBufferFillStencil.hpp"

Stencils::PressureBufferFillStencil::PressureBufferFillStencil(const Parameters& parameters):
  BoundaryStencil(parameters) {
  // Read the pressure values in each of the six (3D) boundary faces of a subdomain
  // and store them consecutively in one-dimensional buffer arrays.

  if (parameters.geometry.dim == 2) {
    // Initialize the vectors if 2D
    pressureLeft_.resize(parameters.parallel.localSize[1]);
    pressureRight_.resize(parameters.parallel.localSize[1]);
    pressureBottom_.resize(parameters.parallel.localSize[0]);
    pressureTop_.resize(parameters.parallel.localSize[0]);
  } else if (parameters.geometry.dim == 3) {
    // Initialize the vectors if 3D
    pressureLeft_.resize(parameters.parallel.localSize[1] * parameters.parallel.localSize[2]);
    pressureRight_.resize(parameters.parallel.localSize[1] * parameters.parallel.localSize[2]);
    pressureBottom_.resize(parameters.parallel.localSize[0] * parameters.parallel.localSize[2]);
    pressureTop_.resize(parameters.parallel.localSize[0] * parameters.parallel.localSize[2]);
    pressureFront_.resize(parameters.parallel.localSize[0] * parameters.parallel.localSize[1]);
    pressureBack_.resize(parameters.parallel.localSize[0] * parameters.parallel.localSize[1]);
  } else {
    throw std::invalid_argument("Unsupported dimensionality: must be 2 or 3.");
  }
}

/* Methods for 2D case */
// iterate "- 2" due to ghost cells in each direction. 
void Stencils::PressureBufferFillStencil::applyLeftWall(FlowField& flowField, int i, int j) {
  pressureLeft_[j - 2] = flowField.getPressure().getScalar(i, j);
  /*
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    std::cout << "i: " << i << " j: " << " " << j << std::endl;
    std::cout << "Filled pressureLeft_[" << j - 2 << "] with value " << pressureLeft_[j - 2] 
    << " == " << flowField.getPressure().getScalar(i, j) << std::endl;
  }
  */
}

void Stencils::PressureBufferFillStencil::applyRightWall(FlowField& flowField, int i, int j) {
  pressureRight_[j - 2] = flowField.getPressure().getScalar(i, j);
  /*
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    std::cout << "i: " << i << " j: " << " " << j << std::endl;
    std::cout << "Filled pressureRight_[" << j - 2 << "] with value " << pressureRight_[j - 2]
    << " == " << flowField.getPressure().getScalar(i, j) << std::endl;
  }
  */
}


void Stencils::PressureBufferFillStencil::applyBottomWall(FlowField& flowField, int i, int j) {
  pressureBottom_[i - 2] = flowField.getPressure().getScalar(i, j);
  /*
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    std::cout << "Filled pressureBottom_[" << i - 2 << "] with value " << pressureBottom_[i - 2]
    << " == " << flowField.getPressure().getScalar(i, j) << std::endl;
  }
  */
}

void Stencils::PressureBufferFillStencil::applyTopWall(FlowField& flowField, int i, int j) {
  pressureTop_[i - 2] = flowField.getPressure().getScalar(i, j);
  /*
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    std::cout << "Filled pressureTop_[" << i - 2 << "] with value " << pressureTop_[i - 2]
    << " == " << flowField.getPressure().getScalar(i, j) << std::endl;
  }
  */
}

/* Methods for 3D case */

void Stencils::PressureBufferFillStencil::applyLeftWall(FlowField& flowField, int i, int j, int k) {
  pressureLeft_[j - 2 + flowField.getCellsY() * (k - 2)] = flowField.getPressure().getScalar(i, j, k);
}

void Stencils::PressureBufferFillStencil::applyRightWall(FlowField& flowField, int i, int j, int k) {
  pressureRight_[j - 2 + flowField.getCellsY() * (k - 2)] = flowField.getPressure().getScalar(i, j, k);
}

void Stencils::PressureBufferFillStencil::applyBottomWall(FlowField& flowField, int i, int j, int k) {
  pressureBottom_[i - 2 + flowField.getCellsX() * (k - 2)] = flowField.getPressure().getScalar(i, j, k);
}

void Stencils::PressureBufferFillStencil::applyTopWall(FlowField& flowField, int i, int j, int k) {
  pressureTop_[i - 2 + flowField.getCellsX() * (k - 2)] = flowField.getPressure().getScalar(i, j, k);
}

void Stencils::PressureBufferFillStencil::applyFrontWall(FlowField& flowField, int i, int j, int k) {
  pressureFront_[i - 2 + flowField.getCellsX() * (j - 2)] = flowField.getPressure().getScalar(i, j, k);
}

void Stencils::PressureBufferFillStencil::applyBackWall(FlowField& flowField, int i, int j, int k) {
  pressureBack_[i - 2 + flowField.getCellsX() * (j - 2)] = flowField.getPressure().getScalar(i, j, k);
}

std::vector<RealType>& Stencils::PressureBufferFillStencil::getpressureLeft(){
  return pressureLeft_;
}
std::vector<RealType>& Stencils::PressureBufferFillStencil::getpressureRight(){
  return pressureRight_;
}
std::vector<RealType>& Stencils::PressureBufferFillStencil::getpressureBottom(){
  return pressureBottom_;
}
std::vector<RealType>& Stencils::PressureBufferFillStencil::getpressureTop(){
  return pressureTop_;
}
std::vector<RealType>& Stencils::PressureBufferFillStencil::getpressureFront(){
  return pressureFront_;
}
std::vector<RealType>& Stencils::PressureBufferFillStencil::getpressureBack(){
  return pressureBack_;
}