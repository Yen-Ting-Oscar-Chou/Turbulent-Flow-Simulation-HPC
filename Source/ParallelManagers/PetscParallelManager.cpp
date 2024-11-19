#include "PetscParallelManager.hpp"

#include "FlowField.hpp"
#include <mpi.h>

PetscParallelManager::PetscParallelManager(const Parameters & parameters, FlowField & flowField): _parameters(parameters),
  flowfield_(flowField) {

  // Get the local sizes in the x, y, and z directions from the parameters.
  // These values define the size of the computational domain in each dimension for the local subdomain.
   int Nx = _parameters.parallel.localSize[0];
   int Ny = _parameters.parallel.localSize[1];
   int Nz = _parameters.parallel.localSize[2];

   constexpr int ghostLayers = 2;
   constexpr int n2D = 2;
   constexpr int n3D = 3;

   // Calculate buffer sizes for 2D and 3D.
   if (_parameters.geometry.dim == 2) {

     // +3 accounts for the ghost cells (boundary and halo regions)
     presBufSize.LR = Ny + ghostLayers;
     presBufSize.TB = Nx + ghostLayers;
     presBufSize.FB = 0;

     velBufSize.LR = n2D * (Ny + ghostLayers);
     velBufSize.TB = n2D * (Nx + ghostLayers);
     velBufSize.FB = 0;
   } else if (_parameters.geometry.dim == 3) {
     presBufSize.LR = (Ny + ghostLayers) * (Nz + ghostLayers);
     presBufSize.TB = (Nx + ghostLayers) * (Nz + ghostLayers);
     presBufSize.FB = (Nx + ghostLayers) * (Ny + ghostLayers);

     velBufSize.LR = n3D * (Ny + ghostLayers) * (Nz + ghostLayers);
     velBufSize.TB = n3D * (Nx + ghostLayers) * (Nz + ghostLayers);
     velBufSize.FB = n3D * (Nx + ghostLayers) * (Ny + ghostLayers);
   }

   pressureBuffers.leftSend = (RealType*) calloc(presBufSize.LR, sizeof(RealType));
   pressureBuffers.rightSend = (RealType*) calloc(presBufSize.LR, sizeof(RealType));
   pressureBuffers.bottomSend = (RealType*) calloc(presBufSize.TB, sizeof(RealType));
   pressureBuffers.topSend = (RealType*) calloc(presBufSize.TB, sizeof(RealType));

   if (_parameters.geometry.dim == 3) {
     pressureBuffers.frontSend = (RealType*) calloc(presBufSize.FB, sizeof(RealType));
     pressureBuffers.backSend = (RealType*) calloc(presBufSize.FB, sizeof(RealType));
   }

   pressureBuffers.leftRecv = (RealType*) calloc(presBufSize.LR, sizeof(RealType));
   pressureBuffers.rightRecv = (RealType*) calloc(presBufSize.LR, sizeof(RealType));
   pressureBuffers.bottomRecv = (RealType*) calloc(presBufSize.TB, sizeof(RealType));
   pressureBuffers.topRecv = (RealType*) calloc(presBufSize.TB, sizeof(RealType));

   if (_parameters.geometry.dim == 3) {
     pressureBuffers.frontRecv = (RealType*) calloc(presBufSize.FB, sizeof(RealType));
     pressureBuffers.backRecv = (RealType*) calloc(presBufSize.FB, sizeof(RealType));
   }

   velocityBuffers.leftSend = (RealType*) calloc(velBufSize.LR, sizeof(RealType));
   velocityBuffers.rightSend = (RealType*) calloc(2 * velBufSize.LR, sizeof(RealType));
   velocityBuffers.bottomSend = (RealType*) calloc(velBufSize.TB, sizeof(RealType));
   velocityBuffers.topSend = (RealType*) calloc(2 * velBufSize.TB, sizeof(RealType));

   if (_parameters.geometry.dim == 3) {
     velocityBuffers.frontSend = (RealType*) calloc(velBufSize.FB, sizeof(RealType));
     velocityBuffers.backSend = (RealType*) calloc(2 * velBufSize.FB, sizeof(RealType));
   }

   velocityBuffers.leftRecv = (RealType*) calloc(2 * velBufSize.LR, sizeof(RealType));
   velocityBuffers.rightRecv = (RealType*) calloc(velBufSize.LR, sizeof(RealType));
   velocityBuffers.bottomRecv = (RealType*) calloc(2 * velBufSize.TB, sizeof(RealType));
   velocityBuffers.topRecv = (RealType*) calloc(velBufSize.TB, sizeof(RealType));

   if (_parameters.geometry.dim == 3) {
     velocityBuffers.frontRecv = (RealType*) calloc(2 * velBufSize.FB, sizeof(RealType));
     velocityBuffers.backRecv = (RealType*) calloc(velBufSize.FB, sizeof(RealType));
   }

   if (_parameters.geometry.dim == 2) {
      // TODO: build stencils
   } else if (_parameters.geometry.dim == 3) {
     // TODO: build stencils
   }

   // TODO: construct iterators
}

PetscParallelManager::~PetscParallelManager() {
  free(pressureBuffers.leftSend);
  free(pressureBuffers.rightSend);
  free(pressureBuffers.bottomSend);
  free(pressureBuffers.topSend);

  free(pressureBuffers.leftRecv);
  free(pressureBuffers.rightRecv);
  free(pressureBuffers.bottomRecv);
  free(pressureBuffers.topRecv);

  if (_parameters.geometry.dim == 3) {
    free(pressureBuffers.frontSend);
    free(pressureBuffers.backSend);
    free(pressureBuffers.frontRecv);
    free(pressureBuffers.backRecv);
  }

  // Free the velocities buffers
  free(velocityBuffers.leftSend);
  free(velocityBuffers.rightSend);
  free(velocityBuffers.bottomSend);
  free(velocityBuffers.topSend);

  free(velocityBuffers.leftRecv);
  free(velocityBuffers.rightRecv);
  free(velocityBuffers.bottomRecv);
  free(velocityBuffers.topRecv);

  if (_parameters.geometry.dim == 3) {
    free(velocityBuffers.frontSend);
    free(velocityBuffers.backSend);
    free(velocityBuffers.frontRecv);
    free(velocityBuffers.backRecv);
  }
}

void PetscParallelManager::communicateVelocities() {
  MPI_Status status;
  MPI_Request send_requestLeft, send_requestRight, send_requestBottom, send_requestTop, send_requestFront, send_requestBack;
  MPI_Request recv_requestLeft, recv_requestRight, recv_requestBottom, recv_requestTop, recv_requestFront, recv_requestBack;

  // TODO: Implement communicateVelocities logic

  MPI_Wait(&send_requestLeft,  &status);
  MPI_Wait(&send_requestRight,  &status);
  MPI_Wait(&send_requestBottom, &status);
  MPI_Wait(&send_requestTop,  &status);
  MPI_Wait(&send_requestFront,  &status);
  MPI_Wait(&send_requestBack, &status);
}

void PetscParallelManager::communicatePressure() {
  MPI_Status status;
  MPI_Request send_requestLeft, send_requestRight, send_requestBottom, send_requestTop, send_requestFront, send_requestBack;
  MPI_Request recv_requestLeft, recv_requestRight, recv_requestBottom, recv_requestTop, recv_requestFront, recv_requestBack;

  // TODO: Implement communicatePressure logic

  MPI_Wait(&send_requestLeft,  &status);
  MPI_Wait(&send_requestRight,  &status);
  MPI_Wait(&send_requestBottom, &status);
  MPI_Wait(&send_requestTop,  &status);
  MPI_Wait(&send_requestFront,  &status);
  MPI_Wait(&send_requestBack, &status);
}
