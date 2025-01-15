#include "StdAfx.hpp"

#include "TurbulentSimulation.hpp"

TurbulentSimulation::TurbulentSimulation(Parameters* parameters, TurbulentFlowField* flowField):
  Simulation(parameters, flowField),
  turbulentField_(flowField),
  turbulentFieldIterator_(),
  turbulentBoundaryIterator_() {}

void TurbulentSimulation::initializeFlowField() {
  Simulation::initializeFlowField();
  std::vector<std::tuple<RealType, RealType, RealType>> coordinateList3D;
  std::vector<std::tuple<RealType, RealType>>           coordinateList2D;
  Stencils::ObstacleCoordinatesStencil                  obstStencil(coordinateList2D, coordinateList3D);
  FieldIterator<FlowField>                              obstIterator(*flowField_, *parameters_, obstStencil);
  GlobalBoundaryIterator<FlowField>                     obstBoundaryIterator(*flowField_, *parameters_, obstStencil, 1, 0);
  obstIterator.iterate();
  obstBoundaryIterator.iterate();
  if (parameters_->geometry.dim == 2) {
    Stencils::DistanceStencil         distStencil(coordinateList2D);
    FieldIterator<TurbulentFlowField> distIterator(*turbulentField_, *parameters_, distStencil);
    distIterator.iterate();
  } else {
    Stencils::DistanceStencil         distStencil(coordinateList3D);
    FieldIterator<TurbulentFlowField> distIterator(*turbulentField_, *parameters_, distStencil);
    distIterator.iterate();
  }
}

void TurbulentSimulation::plotVTK(int timeStep, RealType simulationTime) {
  Stencils::TurbulentVTKStencil     vtkStencil(*parameters_);
  FieldIterator<TurbulentFlowField> vtkIterator(*turbulentField_, *parameters_, vtkStencil, 1, 0);

  vtkIterator.iterate();
  vtkStencil.write(*turbulentField_, timeStep, simulationTime);
}

void TurbulentSimulation::solveTimestep() {
  turbulentFieldIterator_.iterate(VISCOSITY, *parameters_, *turbulentField_, *stencil_);
  // Determine and set max. timestep which is allowed in this simulation
  setTimeStep();
  // Compute FGH
  turbulentFieldIterator_.iterate(TURBFGH, *parameters_, *turbulentField_, *stencil_);

  Simulation::solveTimestepHelper();
}

void TurbulentSimulation::setTimeStep() {
  RealType localMin;
  ASSERTION(parameters_->geometry.dim == 2 || parameters_->geometry.dim == 3);
  RealType factor = 1.0 / (parameters_->meshsize->getDxMin() * parameters_->meshsize->getDxMin()) + 1.0 / (parameters_->meshsize->getDyMin() * parameters_->meshsize->getDyMin());
  // Determine maximum velocity
  stencil_->maxUStencil_.reset();
  stencil_->maxViscStencil_.reset();
  fieldIterator_.iterate(MAXU, *parameters_, *flowField_, *stencil_);
  boundaryIterator_.iterate(MAXU, *parameters_, *flowField_, *stencil_);
  turbulentFieldIterator_.iterate(MAXV, *parameters_, *turbulentField_, *stencil_);
  turbulentBoundaryIterator_.iterate(MAXV, *parameters_, *turbulentField_, *stencil_);
  if (parameters_->geometry.dim == 3) {
    factor += 1.0 / (parameters_->meshsize->getDzMin() * parameters_->meshsize->getDzMin());
  }

  parameters_->timestep.dt = 1.0 / (stencil_->maxUStencil_.getMaxValue() + EPSILON);

  localMin = std::min(1 / (2 * (1 / parameters_->flow.Re + stencil_->maxViscStencil_.getMaxValue()) * factor), parameters_->timestep.dt);

  // Here, we select the type of operation before compiling. This allows to use the correct
  // data type for MPI. Not a concern for small simulations, but useful if using heterogeneous
  // machines.

  parameters_->timestep.dt = localMin;
  parameters_->timestep.dt *= parameters_->timestep.tau;
}