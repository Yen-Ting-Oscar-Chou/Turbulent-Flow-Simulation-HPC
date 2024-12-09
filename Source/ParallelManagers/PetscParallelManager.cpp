#include "StdAfx.hpp"

#include "PetscParallelManager.hpp"

ParallelManagers::PetscParallelManager::PetscParallelManager(Parameters& parameters, FlowField& flowField):
  _flowField(flowField),
  _pressureBufferFillStencil(parameters),
  _pressureBufferReadStencil(parameters),
  _velocityBufferFillStencil(parameters),
  _velocityBufferReadStencil(parameters),
  _parallelBoundaryPressureFillIterator(_flowField, parameters, _pressureBufferFillStencil, 1, 0),
  _parallelBoundaryPressureReadIterator(_flowField, parameters, _pressureBufferReadStencil, 1, 0),
  _parallelBoundaryVelocityFillIterator(_flowField, parameters, _velocityBufferFillStencil, 1, 0),
  _parallelBoundaryVelocityReadIterator(_flowField, parameters, _velocityBufferReadStencil, 1, 0),
  _parameters(parameters) {

  _left   = _parameters.parallel.leftNb;
  _right  = _parameters.parallel.rightNb;
  _top    = _parameters.parallel.topNb;
  _bottom = _parameters.parallel.bottomNb;
  _front  = _parameters.parallel.frontNb;
  _back   = _parameters.parallel.backNb;
}

void ParallelManagers::PetscParallelManager::communicateVelocities() {
  int                         counter = 0;
  std::array<MPI_Request, 12> requests{};

  _parallelBoundaryVelocityFillIterator.iterate();
  // Put recv before send to avoid exceptions
  MPI_Irecv(
    _velocityBufferReadStencil.getBufferLeft().data(), _velocityBufferReadStencil.getBufferLeft().size(), MY_MPI_FLOAT, _left, 1, PETSC_COMM_WORLD, &requests.at(counter++)
  );
  MPI_Irecv(
    _velocityBufferReadStencil.getBufferRight().data(), _velocityBufferReadStencil.getBufferRight().size(), MY_MPI_FLOAT, _right, 1, PETSC_COMM_WORLD, &requests.at(counter++)
  );
  MPI_Irecv(
    _velocityBufferReadStencil.getBufferBottom().data(),
    _velocityBufferReadStencil.getBufferBottom().size(),
    MY_MPI_FLOAT,
    _bottom,
    1,
    PETSC_COMM_WORLD,
    &requests.at(counter++)
  );
  MPI_Irecv(
    _velocityBufferReadStencil.getBufferTop().data(), _velocityBufferReadStencil.getBufferTop().size(), MY_MPI_FLOAT, _top, 1, PETSC_COMM_WORLD, &requests.at(counter++)
  );
  MPI_Irecv(
    _velocityBufferReadStencil.getBufferFront().data(), _velocityBufferReadStencil.getBufferFront().size(), MY_MPI_FLOAT, _front, 1, PETSC_COMM_WORLD, &requests.at(counter++)
  );
  MPI_Irecv(
    _velocityBufferReadStencil.getBufferBack().data(), _velocityBufferReadStencil.getBufferBack().size(), MY_MPI_FLOAT, _back, 1, PETSC_COMM_WORLD, &requests.at(counter++)
  );

  MPI_Isend(
    _velocityBufferFillStencil.getBufferLeft().data(), _velocityBufferFillStencil.getBufferLeft().size(), MY_MPI_FLOAT, _left, 1, PETSC_COMM_WORLD, &requests.at(counter++)
  );
  MPI_Isend(
    _velocityBufferFillStencil.getBufferRight().data(), _velocityBufferFillStencil.getBufferRight().size(), MY_MPI_FLOAT, _right, 1, PETSC_COMM_WORLD, &requests.at(counter++)
  );
  MPI_Isend(
    _velocityBufferFillStencil.getBufferBottom().data(),
    _velocityBufferFillStencil.getBufferBottom().size(),
    MY_MPI_FLOAT,
    _bottom,
    1,
    PETSC_COMM_WORLD,
    &requests.at(counter++)
  );
  MPI_Isend(
    _velocityBufferFillStencil.getBufferTop().data(), _velocityBufferFillStencil.getBufferTop().size(), MY_MPI_FLOAT, _top, 1, PETSC_COMM_WORLD, &requests.at(counter++)
  );
  MPI_Isend(
    _velocityBufferFillStencil.getBufferFront().data(), _velocityBufferFillStencil.getBufferFront().size(), MY_MPI_FLOAT, _front, 1, PETSC_COMM_WORLD, &requests.at(counter++)
  );
  MPI_Isend(
    _velocityBufferFillStencil.getBufferBack().data(), _velocityBufferFillStencil.getBufferBack().size(), MY_MPI_FLOAT, _back, 1, PETSC_COMM_WORLD, &requests.at(counter++)
  );

  MPI_Waitall(counter, requests.data(), MPI_STATUSES_IGNORE);


  _parallelBoundaryVelocityReadIterator.iterate();
}

void ParallelManagers::PetscParallelManager::communicatePressure() {
  int                         counter = 0;
  std::array<MPI_Request, 12> requests{};

  _parallelBoundaryPressureFillIterator.iterate();
  // Put recv before send to avoid exceptions
  MPI_Irecv(
    _pressureBufferReadStencil.getBufferLeft().data(), _pressureBufferReadStencil.getBufferLeft().size(), MY_MPI_FLOAT, _left, 1, PETSC_COMM_WORLD, &requests.at(counter++)
  );
  MPI_Irecv(
    _pressureBufferReadStencil.getBufferRight().data(), _pressureBufferReadStencil.getBufferRight().size(), MY_MPI_FLOAT, _right, 1, PETSC_COMM_WORLD, &requests.at(counter++)
  );
  MPI_Irecv(
    _pressureBufferReadStencil.getBufferBottom().data(),
    _pressureBufferReadStencil.getBufferBottom().size(),
    MY_MPI_FLOAT,
    _bottom,
    1,
    PETSC_COMM_WORLD,
    &requests.at(counter++)
  );
  MPI_Irecv(
    _pressureBufferReadStencil.getBufferTop().data(), _pressureBufferReadStencil.getBufferTop().size(), MY_MPI_FLOAT, _top, 1, PETSC_COMM_WORLD, &requests.at(counter++)
  );
  MPI_Irecv(
    _pressureBufferReadStencil.getBufferFront().data(), _pressureBufferReadStencil.getBufferFront().size(), MY_MPI_FLOAT, _front, 1, PETSC_COMM_WORLD, &requests.at(counter++)
  );
  MPI_Irecv(
    _pressureBufferReadStencil.getBufferBack().data(), _pressureBufferReadStencil.getBufferBack().size(), MY_MPI_FLOAT, _back, 1, PETSC_COMM_WORLD, &requests.at(counter++)
  );

  MPI_Isend(
    _pressureBufferFillStencil.getBufferLeft().data(), _pressureBufferFillStencil.getBufferLeft().size(), MY_MPI_FLOAT, _left, 1, PETSC_COMM_WORLD, &requests.at(counter++)
  );
  MPI_Isend(
    _pressureBufferFillStencil.getBufferRight().data(), _pressureBufferFillStencil.getBufferRight().size(), MY_MPI_FLOAT, _right, 1, PETSC_COMM_WORLD, &requests.at(counter++)
  );
  MPI_Isend(
    _pressureBufferFillStencil.getBufferBottom().data(),
    _pressureBufferFillStencil.getBufferBottom().size(),
    MY_MPI_FLOAT,
    _bottom,
    1,
    PETSC_COMM_WORLD,
    &requests.at(counter++)
  );
  MPI_Isend(
    _pressureBufferFillStencil.getBufferTop().data(), _pressureBufferFillStencil.getBufferTop().size(), MY_MPI_FLOAT, _top, 1, PETSC_COMM_WORLD, &requests.at(counter++)
  );
  MPI_Isend(
    _pressureBufferFillStencil.getBufferFront().data(), _pressureBufferFillStencil.getBufferFront().size(), MY_MPI_FLOAT, _front, 1, PETSC_COMM_WORLD, &requests.at(counter++)
  );
  MPI_Isend(
    _pressureBufferFillStencil.getBufferBack().data(), _pressureBufferFillStencil.getBufferBack().size(), MY_MPI_FLOAT, _back, 1, PETSC_COMM_WORLD, &requests.at(counter++)
  );

  MPI_Waitall(counter, requests.data(), MPI_STATUSES_IGNORE);
  _parallelBoundaryPressureReadIterator.iterate();
}
