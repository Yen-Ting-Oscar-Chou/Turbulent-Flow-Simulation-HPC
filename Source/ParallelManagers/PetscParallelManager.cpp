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
    _velocityBufferReadStencil.getvelocityLeft().data(), _velocityBufferReadStencil.getvelocityLeft().size(), MY_MPI_FLOAT, _left, 1, PETSC_COMM_WORLD, &requests.at(counter++)
  );
  MPI_Irecv(
    _velocityBufferReadStencil.getvelocityRight().data(), _velocityBufferReadStencil.getvelocityRight().size(), MY_MPI_FLOAT, _right, 1, PETSC_COMM_WORLD, &requests.at(counter++)
  );
  MPI_Irecv(
    _velocityBufferReadStencil.getvelocityBottom().data(),
    _velocityBufferReadStencil.getvelocityBottom().size(),
    MY_MPI_FLOAT,
    _bottom,
    1,
    PETSC_COMM_WORLD,
    &requests.at(counter++)
  );
  MPI_Irecv(
    _velocityBufferReadStencil.getvelocityTop().data(), _velocityBufferReadStencil.getvelocityTop().size(), MY_MPI_FLOAT, _top, 1, PETSC_COMM_WORLD, &requests.at(counter++)
  );
  MPI_Irecv(
    _velocityBufferReadStencil.getvelocityFront().data(), _velocityBufferReadStencil.getvelocityFront().size(), MY_MPI_FLOAT, _front, 1, PETSC_COMM_WORLD, &requests.at(counter++)
  );
  MPI_Irecv(
    _velocityBufferReadStencil.getvelocityBack().data(), _velocityBufferReadStencil.getvelocityBack().size(), MY_MPI_FLOAT, _back, 1, PETSC_COMM_WORLD, &requests.at(counter++)
  );

  MPI_Isend(
    _velocityBufferFillStencil.getvelocityLeft().data(), _velocityBufferFillStencil.getvelocityLeft().size(), MY_MPI_FLOAT, _left, 1, PETSC_COMM_WORLD, &requests.at(counter++)
  );
  MPI_Isend(
    _velocityBufferFillStencil.getvelocityRight().data(), _velocityBufferFillStencil.getvelocityRight().size(), MY_MPI_FLOAT, _right, 1, PETSC_COMM_WORLD, &requests.at(counter++)
  );
  MPI_Isend(
    _velocityBufferFillStencil.getvelocityBottom().data(),
    _velocityBufferFillStencil.getvelocityBottom().size(),
    MY_MPI_FLOAT,
    _bottom,
    1,
    PETSC_COMM_WORLD,
    &requests.at(counter++)
  );
  MPI_Isend(
    _velocityBufferFillStencil.getvelocityTop().data(), _velocityBufferFillStencil.getvelocityTop().size(), MY_MPI_FLOAT, _top, 1, PETSC_COMM_WORLD, &requests.at(counter++)
  );
  MPI_Isend(
    _velocityBufferFillStencil.getvelocityFront().data(), _velocityBufferFillStencil.getvelocityFront().size(), MY_MPI_FLOAT, _front, 1, PETSC_COMM_WORLD, &requests.at(counter++)
  );
  MPI_Isend(
    _velocityBufferFillStencil.getvelocityBack().data(), _velocityBufferFillStencil.getvelocityBack().size(), MY_MPI_FLOAT, _back, 1, PETSC_COMM_WORLD, &requests.at(counter++)
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
    _pressureBufferReadStencil.getpressureLeft().data(), _pressureBufferReadStencil.getpressureLeft().size(), MY_MPI_FLOAT, _left, 1, PETSC_COMM_WORLD, &requests.at(counter++)
  );
  MPI_Irecv(
    _pressureBufferReadStencil.getpressureRight().data(), _pressureBufferReadStencil.getpressureRight().size(), MY_MPI_FLOAT, _right, 1, PETSC_COMM_WORLD, &requests.at(counter++)
  );
  MPI_Irecv(
    _pressureBufferReadStencil.getpressureBottom().data(),
    _pressureBufferReadStencil.getpressureBottom().size(),
    MY_MPI_FLOAT,
    _bottom,
    1,
    PETSC_COMM_WORLD,
    &requests.at(counter++)
  );
  MPI_Irecv(
    _pressureBufferReadStencil.getpressureTop().data(), _pressureBufferReadStencil.getpressureTop().size(), MY_MPI_FLOAT, _top, 1, PETSC_COMM_WORLD, &requests.at(counter++)
  );
  MPI_Irecv(
    _pressureBufferReadStencil.getpressureFront().data(), _pressureBufferReadStencil.getpressureFront().size(), MY_MPI_FLOAT, _front, 1, PETSC_COMM_WORLD, &requests.at(counter++)
  );
  MPI_Irecv(
    _pressureBufferReadStencil.getpressureBack().data(), _pressureBufferReadStencil.getpressureBack().size(), MY_MPI_FLOAT, _back, 1, PETSC_COMM_WORLD, &requests.at(counter++)
  );

  MPI_Isend(
    _pressureBufferFillStencil.getpressureLeft().data(), _pressureBufferFillStencil.getpressureLeft().size(), MY_MPI_FLOAT, _left, 1, PETSC_COMM_WORLD, &requests.at(counter++)
  );
  MPI_Isend(
    _pressureBufferFillStencil.getpressureRight().data(), _pressureBufferFillStencil.getpressureRight().size(), MY_MPI_FLOAT, _right, 1, PETSC_COMM_WORLD, &requests.at(counter++)
  );
  MPI_Isend(
    _pressureBufferFillStencil.getpressureBottom().data(),
    _pressureBufferFillStencil.getpressureBottom().size(),
    MY_MPI_FLOAT,
    _bottom,
    1,
    PETSC_COMM_WORLD,
    &requests.at(counter++)
  );
  MPI_Isend(
    _pressureBufferFillStencil.getpressureTop().data(), _pressureBufferFillStencil.getpressureTop().size(), MY_MPI_FLOAT, _top, 1, PETSC_COMM_WORLD, &requests.at(counter++)
  );
  MPI_Isend(
    _pressureBufferFillStencil.getpressureFront().data(), _pressureBufferFillStencil.getpressureFront().size(), MY_MPI_FLOAT, _front, 1, PETSC_COMM_WORLD, &requests.at(counter++)
  );
  MPI_Isend(
    _pressureBufferFillStencil.getpressureBack().data(), _pressureBufferFillStencil.getpressureBack().size(), MY_MPI_FLOAT, _back, 1, PETSC_COMM_WORLD, &requests.at(counter++)
  );

  MPI_Waitall(counter, requests.data(), MPI_STATUSES_IGNORE);
  _parallelBoundaryPressureReadIterator.iterate();
}
