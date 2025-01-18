#include "StdAfx.hpp"

#include "Simulation.hpp"

// #include "Solvers/PetscSolver.hpp"
#include "Solvers/SORSolver.hpp"

Simulation::Simulation(Parameters* parameters, FlowField* flowField):
  parameters_(parameters),
  flowField_(flowField),
  vtkStencil(*parameters),
  vtkIterator(*flowField_, *parameters_, vtkStencil, 1, 0),
  stencil_(new StencilDelegate()),
  fieldIterator_(),
  wallIterator_(1, 0),
  boundaryIterator_()
  // #ifdef ENABLE_PETSC
  //   ,
  //   solver_(std::make_unique<Solvers::PetscSolver>(*flowField_, *parameters))
  // #else
  ,
  solver_()
// #endif
{}

void Simulation::initializeFlowField() {
  if (parameters_->simulation.scenario == CAVITY) {
    parameters_->walls.typeLeft   = DIRICHLET;
    parameters_->walls.typeRight  = DIRICHLET;
    parameters_->walls.typeBottom = DIRICHLET;
    parameters_->walls.typeTop    = DIRICHLET;
    parameters_->walls.typeFront  = DIRICHLET;
    parameters_->walls.typeBack   = DIRICHLET;
  } else if (parameters_->simulation.scenario == CHANNEL) {
    parameters_->walls.typeLeft   = DIRICHLET;
    parameters_->walls.typeRight  = NEUMANN;
    parameters_->walls.typeBottom = DIRICHLET;
    parameters_->walls.typeTop    = DIRICHLET;
    parameters_->walls.typeFront  = DIRICHLET;
    parameters_->walls.typeBack   = DIRICHLET;

    Stencils::BFStepInitStencil bfStepInitStencil(*parameters_);
    FieldIterator<FlowField>    bfStepIterator(*flowField_, *parameters_, bfStepInitStencil, 0, 1);
    bfStepIterator.iterate();
    boundaryIterator_.iterate(WALLVELOCITY, *parameters_, *flowField_, *stencil_);
  }

  solver_.reInitMatrix();
}

void Simulation::solveTimestep() {
  // Determine and set max. timestep which is allowed in this simulation
  setTimeStep();
  // Compute FGH
  fieldIterator_.iterate(FGH, *parameters_, *flowField_, *stencil_);
  solveTimestepHelper();
}

void Simulation::solveTimestepHelper() {
  // Set global boundary values
  wallIterator_.iterate(WALLFGH, *parameters_, *flowField_, *stencil_);
  // Compute the right hand side (RHS)
  fieldIterator_.iterate(RHS, *parameters_, *flowField_, *stencil_);
  // Solve for pressure
  solver_.solve(*flowField_, *parameters_);
  // TODO WS2: communicate pressure values
  // Compute velocity
  fieldIterator_.iterate(VELOCITY, *parameters_, *flowField_, *stencil_);
  fieldIterator_.iterate(OBSTACLE, *parameters_, *flowField_, *stencil_);
  // TODO WS2: communicate velocity values
  // Iterate for velocities on the boundary
  wallIterator_.iterate(WALLVELOCITY, *parameters_, *flowField_, *stencil_);
}

void Simulation::plotVTK(int timeStep, RealType simulationTime) {
  vtkIterator.iterate();
  vtkStencil.write(*flowField_, timeStep, simulationTime);
}

void Simulation::setTimeStep() {
  RealType localMin;
  ASSERTION(parameters_->geometry.dim == 2 || parameters_->geometry.dim == 3);
  RealType factor = 1.0 / (parameters_->meshsize->getDxMin() * parameters_->meshsize->getDxMin()) + 1.0 / (parameters_->meshsize->getDyMin() * parameters_->meshsize->getDyMin());
  // Determine maximum velocity
  stencil_->maxUStencil_.reset();
  fieldIterator_.iterate(MAXU, *parameters_, *flowField_, *stencil_);
  boundaryIterator_.iterate(MAXU, *parameters_, *flowField_, *stencil_);
  if (parameters_->geometry.dim == 3) {
    factor += 1.0 / (parameters_->meshsize->getDzMin() * parameters_->meshsize->getDzMin());
  }

  parameters_->timestep.dt = 1.0 / (stencil_->maxUStencil_.getMaxValue() + EPSILON);

  localMin = std::min(parameters_->flow.Re / (2 * factor), parameters_->timestep.dt);

  // Here, we select the type of operation before compiling. This allows to use the correct
  // data type for MPI. Not a concern for small simulations, but useful if using heterogeneous
  // machines.

  parameters_->timestep.dt = localMin;
  parameters_->timestep.dt *= parameters_->timestep.tau;
}
