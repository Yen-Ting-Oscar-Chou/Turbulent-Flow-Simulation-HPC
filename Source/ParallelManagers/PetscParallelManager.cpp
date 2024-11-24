#include "PetscParallelManager.hpp"

#include <mpi.h>
#include <petsclog.h>

#include "FlowField.hpp"

PetscParallelManager::PetscParallelManager(const Parameters & parameters, FlowField & flowField): flowfield_(flowField),
  _parameters(parameters) {

  // Get the local sizes in the x, y, and z directions from the parameters.
  // These values define the size of the computational domain in each dimension for the local subdomain.
   int Cx = flowField.getCellsX();
   int Cy = flowField.getCellsY();
   int Cz = flowField.getCellsZ();

   int Nx = _parameters.parallel.localSize[0];
   int Ny = _parameters.parallel.localSize[1];
   int Nz = _parameters.parallel.localSize[2];


   // Handles the dimensions internally.
   _pressureBufferFillStencil = new Stencils::PressureBufferFillStencil(parameters);
   _velocityBufferFillStencil = new Stencils::VelocityBufferFillStencil(parameters);

   _pressureBufferReadStencil = new Stencils::PressureBufferReadStencil(parameters);
   _velocityBufferReadStencil = new Stencils::VelocityBufferReadStencil(parameters);

   if (parameters.geometry.dim == 2) {
     pressureBuffers.leftRecv = _pressureBufferReadStencil->getpressureLeft();
     pressureBuffers.leftRecv.resize(Ny);
     pressureBuffers.rightRecv = _pressureBufferReadStencil->getpressureRight();
     pressureBuffers.rightRecv.resize(Ny);
     pressureBuffers.topRecv = _pressureBufferReadStencil->getpressureTop();
     pressureBuffers.topRecv.resize(Nx);
     pressureBuffers.bottomRecv = _pressureBufferReadStencil->getpressureBottom();
     pressureBuffers.bottomRecv.resize(Nx);
     pressureBuffers.frontRecv = _pressureBufferReadStencil->getpressureFront();
     pressureBuffers.frontRecv.resize(0);
     pressureBuffers.backRecv = _pressureBufferReadStencil->getpressureBack();
     pressureBuffers.backRecv.resize(0);

     velocityBuffers.leftRecv = _velocityBufferReadStencil->getvelocityLeft();
     velocityBuffers.leftRecv.resize(2 * Ny);
     velocityBuffers.rightRecv = _velocityBufferReadStencil->getvelocityRight();
     velocityBuffers.rightRecv.resize(2 * Ny);
     velocityBuffers.topRecv = _velocityBufferReadStencil->getvelocityTop();
     velocityBuffers.topRecv.resize(2 * Nx);
     velocityBuffers.bottomRecv = _velocityBufferReadStencil->getvelocityBottom();
     velocityBuffers.bottomRecv.resize(2 * Nx);
     velocityBuffers.frontRecv = _velocityBufferReadStencil->getvelocityFront();
     velocityBuffers.frontRecv.resize(0);
     velocityBuffers.backRecv = _velocityBufferReadStencil->getvelocityBack();
     velocityBuffers.backRecv.resize(0);
   } else if (parameters.geometry.dim == 3) {
     pressureBuffers.leftRecv = _pressureBufferReadStencil->getpressureLeft();
     pressureBuffers.leftRecv.resize(Ny * Nz);
     pressureBuffers.rightRecv = _pressureBufferReadStencil->getpressureRight();
     pressureBuffers.rightRecv.resize(Ny * Nz);
     pressureBuffers.topRecv = _pressureBufferReadStencil->getpressureTop();
     pressureBuffers.topRecv.resize(Nx * Nz);
     pressureBuffers.bottomRecv = _pressureBufferReadStencil->getpressureBottom();
     pressureBuffers.bottomRecv.resize(Nx * Nz);
     pressureBuffers.frontRecv = _pressureBufferReadStencil->getpressureFront();
     pressureBuffers.frontRecv.resize(Nx * Ny);
     velocityBuffers.backRecv = _velocityBufferReadStencil->getvelocityBack();
     velocityBuffers.backRecv.resize(Nx * Ny);

     velocityBuffers.leftRecv = _velocityBufferReadStencil->getvelocityLeft();
     velocityBuffers.leftRecv.resize(3 * Ny * Nz);
     velocityBuffers.rightRecv = _velocityBufferReadStencil->getvelocityRight();
     velocityBuffers.rightRecv.resize(3 * Ny * Nz);
     velocityBuffers.topRecv = _velocityBufferReadStencil->getvelocityTop();
     velocityBuffers.topRecv.resize(3 * Nx * Nz);
     velocityBuffers.bottomRecv = _velocityBufferReadStencil->getvelocityBottom();
     velocityBuffers.bottomRecv.resize(3 * Nx * Nz);
     velocityBuffers.frontRecv = _velocityBufferReadStencil->getvelocityFront();
     velocityBuffers.frontRecv.resize(3 * Nx * Ny);
     velocityBuffers.backRecv = _velocityBufferReadStencil->getvelocityBack();
     velocityBuffers.backRecv.resize(3 * Nx * Ny);
   } else {
     throw std::invalid_argument("Unsupported dimensionality: must be 2 or 3.");
   }

   int lowOffset = 0;
   int highOffset = 0;

   // Construct iterators.
   _parallelBoundaryPressureFillIterator = new ParallelBoundaryIterator<FlowField>(flowField, parameters, *_pressureBufferFillStencil, lowOffset, highOffset);
   _parallelBoundaryPressureReadIterator = new ParallelBoundaryIterator<FlowField>(flowField, parameters, *_pressureBufferReadStencil, lowOffset, highOffset);

   _parallelBoundaryVelocityFillIterator = new ParallelBoundaryIterator<FlowField>(flowField, parameters, *_pressureBufferFillStencil, lowOffset, highOffset);
   _parallelBoundaryVelocityReadIterator = new ParallelBoundaryIterator<FlowField>(flowField, parameters, *_pressureBufferReadStencil, lowOffset, highOffset);
}

PetscParallelManager::~PetscParallelManager() {
}

int PetscParallelManager::computeRankFromIndices(int i, int j, int k) const {
  if (i < 0 || i >= _parameters.parallel.numProcessors[0] ||
      j < 0 || j >= _parameters.parallel.numProcessors[1] ||
      k < 0 || k >= _parameters.parallel.numProcessors[2]) {
    return MPI_PROC_NULL;
  }
  int nrank = i + j * _parameters.parallel.numProcessors[0];
  if (_parameters.geometry.dim == 3) {
    nrank += k * _parameters.parallel.numProcessors[0] * _parameters.parallel.numProcessors[1];
  }
  return nrank;
}

// MPI_Isend(send_buffer.data(), send_buffer.size(), MPI_DOUBLE, dest, dest_tag, _communicator, &requests.at(counter++));
void PetscParallelManager::communicateVelocities() {
  int counter = 0;
  std::array<MPI_Request, 12> requests{};

  int index0 = _parameters.parallel.indices[0];
  int index1 = _parameters.parallel.indices[1];
  int index2 = _parameters.parallel.indices[2];
  int destleft = computeRankFromIndices(index0 - 1, index1, index2);
  int destright = computeRankFromIndices(index0 + 1, index1, index2);
  int desttop = computeRankFromIndices(index0, index1 + 1, index2);
  int destbottom = computeRankFromIndices(index0, index1 - 1, index2);
  int destfront = 0;
  int destback = 0;
  int nproc = _parameters.parallel.numProcessors;
  if (destleft < 0 || destleft >= nproc||
      destright < 0 || destright >= nproc ||
      desttop < 0 || desttop >= nproc ||
      destbottom < 0 || destbottom >= nproc
      ) {
    throw std::logic_error("destination rank out of bounds");
  }

  _parallelBoundaryVelocityFillIterator->iterate();
  // Put recv before send to avoid exceptions
  MPI_Irecv(velocityBuffers.leftRecv.data(), velocityBuffers.leftRecv, MY_MPI_FLOAT, destleft, 1, PETSC_COMM_WORLD, &requests.at(counter++));
  MPI_Irecv(velocityBuffers.rightRecv.data(), velocityBuffers.rightRecv.size(), MY_MPI_FLOAT, destright, 1, PETSC_COMM_WORLD, &requests.at(counter++));
  MPI_Irecv(velocityBuffers.bottomRecv.data(), velocityBuffers.leftRecv.size(), MY_MPI_FLOAT, destbottom, 1, PETSC_COMM_WORLD, &requests.at(counter++));
  MPI_Irecv(velocityBuffers.topRecv.data(), velocityBuffers.leftRecv.size(), MY_MPI_FLOAT, desttop, 1, PETSC_COMM_WORLD, &requests.at(counter++));

  if (_parameters.geometry.dim == 3) {
    destfront = computeRankFromIndices(index0, index1, index2 - 1);
    destback = computeRankFromIndices(index0, index1, index2 + 1);
    MPI_Irecv(velocityBuffers.frontSend.data(), velocityBuffers.frontSend.size(), MY_MPI_FLOAT, destfront, 1, PETSC_COMM_WORLD,&requests.at(counter++));
    MPI_Irecv(velocityBuffers.backSend.data(), velocityBuffers.backSend.size(), MY_MPI_FLOAT, destback, 1, PETSC_COMM_WORLD,&requests.at(counter++));
    if (destfront < 0 || destleft >= nproc||
        destback < 0 || destright >= nproc
        ) {
      throw std::logic_error("destination rank out of bounds");
    }
  }


  MPI_Isend(velocityBuffers.rightSend.data(), velocityBuffers.rightSend.size(), MY_MPI_FLOAT, destleft, 1, PETSC_COMM_WORLD, &requests.at(counter++));
  MPI_Isend(velocityBuffers.leftSend.data(), velocityBuffers.leftSend.size(), MY_MPI_FLOAT, destright, 1, PETSC_COMM_WORLD, &requests.at(counter++));
  MPI_Isend(velocityBuffers.topSend.data(), velocityBuffers.topSend.size(), MY_MPI_FLOAT, destbottom, 1, PETSC_COMM_WORLD, &requests.at(counter++));
  MPI_Isend(velocityBuffers.bottomSend.data(), velocityBuffers.bottomSend.size(), MY_MPI_FLOAT, desttop, 1, PETSC_COMM_WORLD, &requests.at(counter++));

  if (_parameters.geometry.dim == 3) {
    destfront = computeRankFromIndices(index0, index1, index2 - 1);
    destback = computeRankFromIndices(index0, index1, index2 + 1);
    MPI_Isend(velocityBuffers.frontSend.data(), velocityBuffers.frontSend.size(), MY_MPI_FLOAT, destfront, 1, PETSC_COMM_WORLD, &requests.at(counter++));
    MPI_Isend(velocityBuffers.backSend.data(), velocityBuffers.backSend.size(), MY_MPI_FLOAT, destback, 1, PETSC_COMM_WORLD, &requests.at(counter++));
  }

  MPI_Waitall(counter, requests.data(), MPI_STATUSES_IGNORE);
  _parallelBoundaryVelocityReadIterator->iterate();
}

void PetscParallelManager::communicatePressure() {
  int counter = 0;
  std::array<MPI_Request, 12> requests{};

  int index0 = _parameters.parallel.indices[0];
  int index1 = _parameters.parallel.indices[1];
  int index2 = _parameters.parallel.indices[2];
  int destleft = computeRankFromIndices(index0 - 1, index1, index2);
  int destright = computeRankFromIndices(index0 + 1, index1, index2);
  int desttop = computeRankFromIndices(index0, index1 + 1, index2);
  int destbottom = computeRankFromIndices(index0, index1 - 1, index2);
  int destfront = 0;
  int destback = 0;
  int nproc = _parameters.parallel.numProcessors;
  if (destleft < 0 || destleft >= nproc||
      destright < 0 || destright >= nproc ||
      desttop < 0 || desttop >= nproc ||
      destbottom < 0 || destbottom >= nproc
  ) {
    throw std::logic_error("destination rank out of bounds");
  }
  _parallelBoundaryPressureFillIterator->iterate();
  // Put recv before send to avoid exceptions
  MPI_Irecv(pressureBuffers.leftRecv.data(), pressureBuffers.leftRecv, MY_MPI_FLOAT, destleft, 1, PETSC_COMM_WORLD, &requests.at(counter++));
  MPI_Irecv(pressureBuffers.rightRecv.data(), pressureBuffers.rightRecv.size(), MY_MPI_FLOAT, destright, 1, PETSC_COMM_WORLD, &requests.at(counter++));
  MPI_Irecv(pressureBuffers.bottomRecv.data(), pressureBuffers.leftRecv.size(), MY_MPI_FLOAT, destbottom, 1, PETSC_COMM_WORLD, &requests.at(counter++));
  MPI_Irecv(pressureBuffers.topRecv.data(), pressureBuffers.leftRecv.size(), MY_MPI_FLOAT, desttop, 1, PETSC_COMM_WORLD, &requests.at(counter++));

  if (_parameters.geometry.dim == 3) {
    destfront = computeRankFromIndices(index0, index1, index2 - 1);
    destback = computeRankFromIndices(index0, index1, index2 + 1);
    MPI_Irecv(pressureBuffers.frontSend.data(), pressureBuffers.frontSend.size(), MY_MPI_FLOAT, destfront, 1, PETSC_COMM_WORLD, &requests.at(counter++));
    MPI_Irecv(pressureBuffers.backSend.data(), pressureBuffers.backSend.size(), MY_MPI_FLOAT, destback, 1, PETSC_COMM_WORLD, &requests.at(counter++));
    if (destfront < 0 || destleft >= nproc||
        destback < 0 || destright >= nproc
    ) {
      throw std::logic_error("destination rank out of bounds");
    }
  }

  MPI_Isend(pressureBuffers.rightSend.data(), pressureBuffers.rightSend.size(), MY_MPI_FLOAT, destleft, 1, PETSC_COMM_WORLD, &requests.at(counter++));
  MPI_Isend(pressureBuffers.leftSend.data(), pressureBuffers.leftSend.size(), MY_MPI_FLOAT, destright, 1, PETSC_COMM_WORLD, &requests.at(counter++));
  MPI_Isend(pressureBuffers.topSend.data(), pressureBuffers.topSend.size(), MY_MPI_FLOAT, destbottom, 1, PETSC_COMM_WORLD, &requests.at(counter++));
  MPI_Isend(pressureBuffers.bottomSend.data(), pressureBuffers.bottomSend.size(), MY_MPI_FLOAT, desttop, 1, PETSC_COMM_WORLD, &requests.at(counter++));

  if (_parameters.geometry.dim == 3) {
    destfront = computeRankFromIndices(index0, index1, index2 - 1);
    destback = computeRankFromIndices(index0, index1, index2 + 1);
    MPI_Isend(pressureBuffers.frontSend.data(), pressureBuffers.frontSend.size(), MY_MPI_FLOAT, destfront, 1, PETSC_COMM_WORLD, &requests.at(counter++));
    MPI_Isend(pressureBuffers.backSend.data(), pressureBuffers.backSend.size(), MY_MPI_FLOAT, destback, 1, PETSC_COMM_WORLD, &requests.at(counter++));
  }

  MPI_Waitall(counter, requests.data(), MPI_STATUSES_IGNORE);
  _parallelBoundaryPressureReadIterator->iterate();
}

