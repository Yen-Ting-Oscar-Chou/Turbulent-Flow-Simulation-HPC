#include "StdAfx.hpp"

#include "TurbulentPetscParallelManager.hpp"

ParallelManagers::TurbulentPetscParallelManager::TurbulentPetscParallelManager(Parameters& parameters, TurbulentFlowField& flowField):
  PetscParallelManager(parameters, flowField),
  _turbulentFlowField(flowField),
  _viscosityBufferFillStencil(parameters),
  _viscosityBufferReadStencil(parameters),
  _parallelBoundaryViscosityFillIterator(_turbulentFlowField, parameters, _viscosityBufferFillStencil, 2, -1),
  _parallelBoundaryViscosityReadIterator(_turbulentFlowField, parameters, _viscosityBufferReadStencil, 2, -1) {
}

void ParallelManagers::TurbulentPetscParallelManager::communicateViscosity() {
  int                         counter = 0;
  std::array<MPI_Request, 12> requests{};

  _parallelBoundaryViscosityFillIterator.iterate();

  MPI_Irecv(
    _viscosityBufferReadStencil.getviscosityLeft().data(), _viscosityBufferReadStencil.getviscosityLeft().size(), MY_MPI_FLOAT, _left, 1, PETSC_COMM_WORLD, &requests.at(counter++)
  );
  MPI_Irecv(
    _viscosityBufferReadStencil.getviscosityRight().data(),
    _viscosityBufferReadStencil.getviscosityRight().size(),
    MY_MPI_FLOAT,
    _right,
    1,
    PETSC_COMM_WORLD,
    &requests.at(counter++)
  );
  MPI_Irecv(
    _viscosityBufferReadStencil.getviscosityBottom().data(),
    _viscosityBufferReadStencil.getviscosityBottom().size(),
    MY_MPI_FLOAT,
    _bottom,
    1,
    PETSC_COMM_WORLD,
    &requests.at(counter++)
  );
  MPI_Irecv(
    _viscosityBufferReadStencil.getviscosityTop().data(), _viscosityBufferReadStencil.getviscosityTop().size(), MY_MPI_FLOAT, _top, 1, PETSC_COMM_WORLD, &requests.at(counter++)
  );
  MPI_Irecv(
    _viscosityBufferReadStencil.getviscosityFront().data(),
    _viscosityBufferReadStencil.getviscosityFront().size(),
    MY_MPI_FLOAT,
    _front,
    1,
    PETSC_COMM_WORLD,
    &requests.at(counter++)
  );
  MPI_Irecv(
    _viscosityBufferReadStencil.getviscosityBack().data(), _viscosityBufferReadStencil.getviscosityBack().size(), MY_MPI_FLOAT, _back, 1, PETSC_COMM_WORLD, &requests.at(counter++)
  );

  MPI_Isend(
    _viscosityBufferFillStencil.getviscosityLeft().data(), _viscosityBufferFillStencil.getviscosityLeft().size(), MY_MPI_FLOAT, _left, 1, PETSC_COMM_WORLD, &requests.at(counter++)
  );
  MPI_Isend(
    _viscosityBufferFillStencil.getviscosityRight().data(),
    _viscosityBufferFillStencil.getviscosityRight().size(),
    MY_MPI_FLOAT,
    _right,
    1,
    PETSC_COMM_WORLD,
    &requests.at(counter++)
  );
  MPI_Isend(
    _viscosityBufferFillStencil.getviscosityBottom().data(),
    _viscosityBufferFillStencil.getviscosityBottom().size(),
    MY_MPI_FLOAT,
    _bottom,
    1,
    PETSC_COMM_WORLD,
    &requests.at(counter++)
  );
  MPI_Isend(
    _viscosityBufferFillStencil.getviscosityTop().data(), _viscosityBufferFillStencil.getviscosityTop().size(), MY_MPI_FLOAT, _top, 1, PETSC_COMM_WORLD, &requests.at(counter++)
  );
  MPI_Isend(
    _viscosityBufferFillStencil.getviscosityFront().data(),
    _viscosityBufferFillStencil.getviscosityFront().size(),
    MY_MPI_FLOAT,
    _front,
    1,
    PETSC_COMM_WORLD,
    &requests.at(counter++)
  );
  MPI_Isend(
    _viscosityBufferFillStencil.getviscosityBack().data(), _viscosityBufferFillStencil.getviscosityBack().size(), MY_MPI_FLOAT, _back, 1, PETSC_COMM_WORLD, &requests.at(counter++)
  );

  MPI_Waitall(counter, requests.data(), MPI_STATUSES_IGNORE);

  _parallelBoundaryViscosityReadIterator.iterate();
}
