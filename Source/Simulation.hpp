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
  Simulation*         simulationGPU_;
  FlowFieldGPUPtrs    flowFieldGPUptrs_;
  StencilDelegatePtrs stencilDelegatePtrs_;

  SimulationPtrs(Simulation* simulationGPU, FlowFieldGPUPtrs flowFieldGPUPtrs, StencilDelegatePtrs stencilDelegatePtrs):
    simulationGPU_(simulationGPU),
    flowFieldGPUptrs_(flowFieldGPUPtrs),
    stencilDelegatePtrs_(stencilDelegatePtrs) {}
};

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

  static SimulationPtrs mapToGPU(int hostDevice, int targetDevice, Simulation& simulation) {
    size_t      simulationSize = sizeof(Simulation);
    Simulation* simulationGPU  = static_cast<Simulation*>(omp_target_alloc(simulationSize, targetDevice));
    if (!simulationGPU) {
      std::cout << "Failed to allocate memory for simulation on GPU" << std::endl;
    }

    bool associatedSimulation = omp_target_associate_ptr(&simulation, simulationGPU, simulationSize, 0, targetDevice) == 0;
    if (!associatedSimulation) {
      std::cout << "Failed to associate simulation with GPU memory" << std::endl;
    }

    FlowFieldGPUPtrs    flowFieldPtrs       = FlowField::mapToGPU(hostDevice, targetDevice, simulation.flowField_);
    StencilDelegatePtrs stencilDelegatePtrs = StencilDelegate::mapToGPU(hostDevice, targetDevice, simulation.stencil_);

    bool copiedFlowFieldPtr = omp_target_memcpy(&simulationGPU->flowField_, flowFieldPtrs.flowFieldPtr_, sizeof(FlowField*), 0, 0, targetDevice, hostDevice) == 0;
    if (!copiedFlowFieldPtr) {
      std::cout << "Failed to copy FlowField pointer to Simulation object" << std::endl;
    }
    bool copiedStencilDelegatePtr = omp_target_memcpy(&simulationGPU->stencil_, stencilDelegatePtrs.stencilDelegateGPU_, sizeof(StencilDelegate*), 0, 0, targetDevice, hostDevice)
                                    == 0;
    if (!copiedStencilDelegatePtr) {
      std::cout << "Failed to copy StencilDelegate pointer to Simulation object" << std::endl;
    }

    bool copiedIteratorStencilDelegatePtr
      = omp_target_memcpy(&simulationGPU->fieldIterator_.stencil_, stencilDelegatePtrs.stencilDelegateGPU_, sizeof(StencilDelegate*), 0, 0, targetDevice, hostDevice) == 0;
    if (!copiedIteratorStencilDelegatePtr) {
      std::cout << "Failed to copy StencilDelegate pointer to Iterator object" << std::endl;
    }

    return SimulationPtrs(simulationGPU, flowFieldPtrs, stencilDelegatePtrs);
  }

  static void mapToCPU(int hostDevice, int targetDevice, Simulation& simulation, SimulationPtrs& simulationPtrs) {
    FlowField::mapToCPU(hostDevice, targetDevice, simulation.flowField_, simulationPtrs.flowFieldGPUptrs_);
    //StencilDelegate::mapToCPU(hostDevice, targetDevice, simulation.stencil_, simulationPtrs.stencilDelegatePtrs_);
  }

  static void mapToCPUAndFree(int hostDevice, int targetDevice, Simulation& simulation, SimulationPtrs& simulationPtrs) {
    bool disassociatedSimulation = omp_target_disassociate_ptr(&simulation, targetDevice) == 0;
    if (!disassociatedSimulation) {
      std::cout << "Error: Could not disassociate Simulation pointer." << std::endl;
    }

    FlowField::mapToCPUAndFree(hostDevice, targetDevice, simulation.flowField_, simulationPtrs.flowFieldGPUptrs_);
    StencilDelegate::mapToCPUAndFree(hostDevice, targetDevice, simulation.stencil_, simulationPtrs.stencilDelegatePtrs_);
  }

  #pragma omp declare target
  void test();
  #pragma omp end declare target
};
