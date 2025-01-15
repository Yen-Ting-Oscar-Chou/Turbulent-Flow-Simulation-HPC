#pragma once

#include "Definitions.hpp"
#include "FlowField.hpp"
#include "GlobalBoundaryFactory.hpp"
#include "Iterators.hpp"
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

class Simulation;

struct SimulationPtrs {
  Parameters* parametersPtr_;
  FlowField* flowFieldPtr_;

  SimulationPtrs(Parameters* parametersPtr, FlowField* flowFieldPtr)
    : parametersPtr_(parametersPtr), flowFieldPtr_(flowFieldPtr) {}
};;

class Simulation {
protected:
  Parameters& parameters_;

  FlowField& flowField_;

  std::unique_ptr<Solvers::LinearSolver> solver_;

  Stencils::VTKStencil<FlowField> vtkStencil;
  FieldIterator<FlowField>        vtkIterator;

  StencilDelegate                      stencil_;
  FieldIteratorGPU<FlowField>          fieldIterator_;
  GlobalBoundaryIteratorGPU<FlowField> boundaryIterator_;
  GlobalBoundaryIteratorGPU<FlowField> wallIterator_;

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

  /* static SimulationPtrs mapToGPU(int hostDevice, int targetDevice, Simulation& simulation) {

  }

  static void mapToCPU(int hostDevice, int targetDevice, SimulationPtrs& simulationPtrs) {

  }

  static void mapToCPUAndFree(int hostDevice, int targetDevice, SimulationPtrs& simulationPtrs) {

  } */
};
