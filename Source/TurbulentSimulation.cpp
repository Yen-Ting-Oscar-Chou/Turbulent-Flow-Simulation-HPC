#include "StdAfx.hpp"

#include "TurbulentSimulation.hpp"

TurbulentSimulation::TurbulentSimulation(Parameters& parameters, TurbulentFlowField& flowField):
  Simulation(parameters, flowField),
  turbulentField_(flowField),
  rhsStencil__(RHS),
  rhsIterator__(flowField_, parameters, rhsStencil__),
  velocityStencil__(VELOCITY),
  velocityIterator__(flowField_, parameters, velocityStencil__),
  viscosityStencil_(VISCOSITY),
  // viscosityStencil_(),
  viscosityIterator_(turbulentField_, parameters, viscosityStencil_),
  turbulentFGHStencil_(TURBFGH),
  // turbulentFGHStencil_(),
  turbulentFGHIterator_(turbulentField_, parameters, turbulentFGHStencil_),
  maxViscStencil_(parameters),
  maxViscFieldIterator_(turbulentField_, parameters, maxViscStencil_),
  maxViscBoundaryIterator_(turbulentField_, parameters, maxViscStencil_),
  turbulentPetscParallelManager_(parameters, turbulentField_) {}

void TurbulentSimulation::initializeFlowField() {
  Simulation::initializeFlowField();
  std::vector<std::tuple<RealType, RealType, RealType>> coordinateList3D;
  std::vector<std::tuple<RealType, RealType>>           coordinateList2D;
  Stencils::ObstacleCoordinatesStencil                  obstStencil(parameters_, coordinateList2D, coordinateList3D);
  FieldIterator<FlowField>                              obstIterator(flowField_, parameters_, obstStencil);
  GlobalBoundaryIterator<FlowField>                     obstBoundaryIterator(flowField_, parameters_, obstStencil, 1, 0);
  obstIterator.iterate();
  obstBoundaryIterator.iterate();
  turbulentPetscParallelManager_.communicateObstacleCoordinates(coordinateList2D, coordinateList3D);
  if (parameters_.geometry.dim == 2) {
    Stencils::DistanceStencil         distStencil(coordinateList2D);
    FieldIterator<TurbulentFlowField> distIterator(turbulentField_, parameters_, distStencil);
    distIterator.iterate();
  } else {
    Stencils::DistanceStencil         distStencil(coordinateList3D);
    FieldIterator<TurbulentFlowField> distIterator(turbulentField_, parameters_, distStencil);
    distIterator.iterate();
  }
}

void TurbulentSimulation::plotVTK(int timeStep, RealType simulationTime) {
  Stencils::TurbulentVTKStencil     vtkStencil(parameters_);
  FieldIterator<TurbulentFlowField> vtkIterator(turbulentField_, parameters_, vtkStencil, 1, 0);

  vtkIterator.iterate();
  vtkStencil.write(turbulentField_, timeStep, simulationTime);
}

void TurbulentSimulation::solveTimestep() {
  viscosityIterator_.iterate();
  // TODO communicate viscosity
  turbulentPetscParallelManager_.communicateViscosity();
  // Determine and set max. timestep which is allowed in this simulation
  setTimeStep();
  // Compute FGH
  turbulentFGHIterator_.iterate();
  
  // Set global boundary values
  Simulation::wallFGHIterator_.iterate();
  // Compute the right hand side (RHS)
  rhsIterator__.iterate();
  // Solve for pressure
  Simulation::solver_->solve();
  // TODO WS2: communicate pressure values
  Simulation::petscParallelManager_.communicatePressure();
  // Compute velocity
  velocityIterator__.iterate();
  Simulation::obstacleIterator_.iterate();
  // TODO WS2: communicate velocity values
  Simulation::petscParallelManager_.communicateVelocities();
  // Iterate for velocities on the boundary
  Simulation::wallVelocityIterator_.iterate();
}

void TurbulentSimulation::setTimeStep() {
  RealType localMin, globalMin;
  ASSERTION(parameters_.geometry.dim == 2 || parameters_.geometry.dim == 3);
  RealType factor = 1.0 / (parameters_.meshsize.getDxMin() * parameters_.meshsize.getDxMin()) + 1.0 / (parameters_.meshsize.getDyMin() * parameters_.meshsize.getDyMin());
  // Determine maximum velocity
  maxUStencil_.reset();
  maxViscStencil_.reset();
  maxUFieldIterator_.iterate();
  maxUBoundaryIterator_.iterate();
  maxViscFieldIterator_.iterate();
  maxViscBoundaryIterator_.iterate();
  if (parameters_.geometry.dim == 3) {
    factor += 1.0 / (parameters_.meshsize.getDzMin() * parameters_.meshsize.getDzMin());
    parameters_.timestep.dt = 1.0 / (maxUStencil_.getMaxValues()[2] + EPSILON);
  } else {
    parameters_.timestep.dt = 1.0 / (maxUStencil_.getMaxValues()[0] + EPSILON);
  }

  // localMin = std::min(parameters_.timestep.dt, std::min(std::min(parameters_.flow.Re/(2 * factor), 1.0 /
  // maxUStencil_.getMaxValues()[0]), 1.0 / maxUStencil_.getMaxValues()[1]));
  localMin = std::min(
    1 / (2 * (1 / parameters_.flow.Re + maxViscStencil_.getMaxValue()) * factor),
    std::min(parameters_.timestep.dt, std::min(1 / (maxUStencil_.getMaxValues()[0] + EPSILON), 1 / (maxUStencil_.getMaxValues()[1] + EPSILON)))
  );

  // Here, we select the type of operation before compiling. This allows to use the correct
  // data type for MPI. Not a concern for small simulations, but useful if using heterogeneous
  // machines.

  globalMin = MY_FLOAT_MAX;
  MPI_Allreduce(&localMin, &globalMin, 1, MY_MPI_FLOAT, MPI_MIN, PETSC_COMM_WORLD);

  parameters_.timestep.dt = globalMin;
  parameters_.timestep.dt *= parameters_.timestep.tau;
}