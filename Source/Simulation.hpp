#pragma once

#include "Definitions.hpp"
#include "FlowField.hpp"
#include "GlobalBoundaryFactory.hpp"
#include "Iterators.hpp"
#include "IteratorsGPU.hpp"
#include "ParallelManagers/PetscParallelManager.hpp"
#include "Solvers/LinearSolver.hpp"
#include "Stencils/BFInputStencils.hpp"
#include "Stencils/BFStepInitStencil.hpp"
#include "Stencils/FGHStencil.hpp"
#include "Stencils/InitTaylorGreenFlowFieldStencil.hpp"
#include "Stencils/MaxUStencil.hpp"
#include "Stencils/MovingWallStencils.hpp"
#include "Stencils/NeumannBoundaryStencils.hpp"
#include "Stencils/ObstacleStencil.hpp"
#include "Stencils/PeriodicBoundaryStencils.hpp"
#include "Stencils/RHSStencil.hpp"
#include "Stencils/VelocityStencil.hpp"
#include "Stencils/VTKStencil.hpp"


class Simulation {
protected:
  Parameters& parameters_;

  FlowField& flowField_;

  Stencils::MaxUStencil             maxUStencil_;
  FieldIterator<FlowField>          maxUFieldIterator_;
  GlobalBoundaryIterator<FlowField> maxUBoundaryIterator_;

  // Set up the boundary conditions
  GlobalBoundaryFactory             globalBoundaryFactory_;
  GlobalBoundaryIterator<FlowField> wallVelocityIterator_;
  GlobalBoundaryIterator<FlowField> wallFGHIterator_;

  StencilDelegate  fghStencil_;
  GPUFieldIterator fghIterator_;

  StencilDelegate  rhsStencil_;
  GPUFieldIterator rhsIterator_;

  StencilDelegate           velocityStencil_;
  Stencils::ObstacleStencil obstacleStencil_;
  GPUFieldIterator          velocityIterator_;
  FieldIterator<FlowField>  obstacleIterator_;

  ParallelManagers::PetscParallelManager petscParallelManager_;
  std::unique_ptr<Solvers::LinearSolver> solver_;

  Stencils::VTKStencil<FlowField> vtkStencil;
  FieldIterator<FlowField>        vtkIterator;

  virtual void setTimeStep();

  void solveTimestepHelper();

public:
  Simulation(Parameters& parameters, FlowField& flowField);
  virtual ~Simulation() = default;

  /** Initialises the flow field according to the scenario */
  virtual void initializeFlowField();

  virtual void solveTimestep();

  /** Plots the flow field */
  virtual void plotVTK(int timeStep, RealType simulationTime);
};
