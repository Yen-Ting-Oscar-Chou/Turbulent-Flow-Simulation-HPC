#include "StdAfx.hpp"

#include "Simulation.hpp"

#include "Solvers/PetscSolver.hpp"
#include "Solvers/SORSolver.hpp"

Simulation::Simulation(Parameters& parameters, FlowField& flowField):
  parameters_(parameters),
  flowField_(flowField),
  vtkStencil(parameters),
  vtkIterator(flowField_, parameters_, vtkStencil, 1, 0),
  petscParallelManager_(parameters_, flowField_),
  stencil_(),
  fieldIterator_(flowField_, parameters_, stencil_),
  wallIterator_(flowField_, parameters_, stencil_, 1, 0),
  boundaryIterator_(flowField_, parameters_, stencil_)
#ifdef ENABLE_PETSC
  ,
  solver_(std::make_unique<Solvers::PetscSolver>(flowField_, parameters))
#else
  ,
  solver_(std::make_unique<Solvers::SORSolver>(flowField_, parameters_))
#endif
{
}

void Simulation::initializeFlowField() {
  if (parameters_.simulation.scenario == CAVITY) {
    parameters_.walls.typeLeft   = DIRICHLET;
    parameters_.walls.typeRight  = DIRICHLET;
    parameters_.walls.typeBottom = DIRICHLET;
    parameters_.walls.typeTop    = DIRICHLET;
    parameters_.walls.typeFront  = DIRICHLET;
    parameters_.walls.typeBack   = DIRICHLET;
  } else if (parameters_.simulation.scenario == CHANNEL) {
    parameters_.walls.typeLeft   = DIRICHLET;
    parameters_.walls.typeRight  = NEUMANN;
    parameters_.walls.typeBottom = DIRICHLET;
    parameters_.walls.typeTop    = DIRICHLET;
    parameters_.walls.typeFront  = DIRICHLET;
    parameters_.walls.typeBack   = DIRICHLET;

    Stencils::BFStepInitStencil bfStepInitStencil(parameters_);
    FieldIterator<FlowField>    bfStepIterator(flowField_, parameters_, bfStepInitStencil, 0, 1);
    bfStepIterator.iterate();
    // wallVelocityIterator_.iterate();
    boundaryIterator_.iterate(WALLVELOCITY);
  }

  solver_->reInitMatrix();
}

void Simulation::solveTimestep() {
  // Determine and set max. timestep which is allowed in this simulation
  setTimeStep();
  // Compute FGH
  fieldIterator_.iterate(FGH);
  solveTimestepHelper();
}


void Simulation::solveTimestepHelper() {
  // Set global boundary values
  wallIterator_.iterate(WALLFGH);
  // Compute the right hand side (RHS)
  fieldIterator_.iterate(RHS);
  // Solve for pressure
  solver_->solve();
  // TODO WS2: communicate pressure values
  petscParallelManager_.communicatePressure();
  // Compute velocity
  fieldIterator_.iterate(VELOCITY);
  // obstacleIterator_.iterate();
  fieldIterator_.iterate(OBSTACLE);
  // TODO WS2: communicate velocity values
  petscParallelManager_.communicateVelocities();
  // Iterate for velocities on the boundary
  wallIterator_.iterate(WALLVELOCITY);
}

void Simulation::plotVTK(int timeStep, RealType simulationTime) {
  vtkIterator.iterate();
  vtkStencil.write(flowField_, timeStep, simulationTime);
}

void Simulation::setTimeStep() {
  RealType localMin, globalMin;
  ASSERTION(parameters_.geometry.dim == 2 || parameters_.geometry.dim == 3);
  RealType factor = 1.0 / (parameters_.meshsize.getDxMin() * parameters_.meshsize.getDxMin()) + 1.0 / (parameters_.meshsize.getDyMin() * parameters_.meshsize.getDyMin());
  // Determine maximum velocity
  stencil_.maxUStencil_.reset();
  fieldIterator_.iterate(MAXU);
  boundaryIterator_.iterate(MAXU);
  if (parameters_.geometry.dim == 3) {
    factor += 1.0 / (parameters_.meshsize.getDzMin() * parameters_.meshsize.getDzMin());
  }

  parameters_.timestep.dt = 1.0 / (stencil_.maxUStencil_.getMaxValue() + EPSILON);

  localMin = std::min(parameters_.flow.Re / (2 * factor), parameters_.timestep.dt);

  // Here, we select the type of operation before compiling. This allows to use the correct
  // data type for MPI. Not a concern for small simulations, but useful if using heterogeneous
  // machines.

  globalMin = MY_FLOAT_MAX;
  MPI_Allreduce(&localMin, &globalMin, 1, MY_MPI_FLOAT, MPI_MIN, PETSC_COMM_WORLD);

  parameters_.timestep.dt = globalMin;
  parameters_.timestep.dt *= parameters_.timestep.tau;
}
