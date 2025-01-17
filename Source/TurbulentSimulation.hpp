#pragma once

#include "ParallelManagers/TurbulentPetscParallelManager.hpp"
#include "Simulation.hpp"
#include "Stencils/DistanceStencil.hpp"
#include "Stencils/MaxViscStencil.hpp"
#include "Stencils/ObstacleCoordinatesStencil.hpp"
#include "Stencils/TurbulentFGHStencil.hpp"
#include "Stencils/TurbulentVTKStencil.hpp"
#include "Stencils/ViscosityStencil.hpp"
#include "TurbulentFlowField.hpp"

class TurbulentSimulation;

struct TurbulentSimulationPtrs {
  TurbulentSimulation*   simulationGPU_;
  TurbulentFlowFieldPtrs flowFieldGPUptrs_;
  StencilDelegatePtrs    stencilDelegatePtrs_;
  ParametersGPUPtrs      parameterPtrs_;

  TurbulentSimulationPtrs(TurbulentSimulation* simulationGPU, TurbulentFlowFieldPtrs flowFieldGPUPtrs, StencilDelegatePtrs stencilDelegatePtrs, ParametersGPUPtrs parameterGPUPtrs):
    simulationGPU_(simulationGPU),
    flowFieldGPUptrs_(flowFieldGPUPtrs),
    stencilDelegatePtrs_(stencilDelegatePtrs),
    parameterPtrs_(parameterGPUPtrs) {}
};

class TurbulentSimulation: public Simulation {
public:
  TurbulentFlowField*                           turbulentField_;
  FieldIteratorGPU<TurbulentFlowField>          turbulentFieldIterator_;
  GlobalBoundaryIteratorGPU<TurbulentFlowField> turbulentBoundaryIterator_;

  TurbulentSimulation(Parameters* parameters, TurbulentFlowField* flowField);
  virtual ~TurbulentSimulation() = default;

  /** Initialises the flow field according to the scenario */
  virtual void initializeFlowField() override;

#pragma omp declare target
  void solveTimestep();
  void setTimeStep();
#pragma omp end declare target

  /** Plots the flow field */
  virtual void plotVTK(int timeStep, RealType simulationTime) override;

  static TurbulentSimulationPtrs mapToGPU(int hostDevice, int targetDevice, TurbulentSimulation& simulation) {
    size_t               simulationSize = sizeof(TurbulentSimulation);
    TurbulentSimulation* simulationGPU  = static_cast<TurbulentSimulation*>(omp_target_alloc(simulationSize, targetDevice));
    if (!simulationGPU) {
      std::cout << "Failed to allocate memory for turbulentSimulation on GPU" << std::endl;
    }

    bool associatedSimulation = omp_target_associate_ptr(&simulation, simulationGPU, simulationSize, 0, targetDevice) == 0;
    if (!associatedSimulation) {
      std::cout << "Failed to associate turbulentSimulation with GPU memory" << std::endl;
    }

    TurbulentFlowFieldPtrs flowFieldPtrs       = TurbulentFlowField::mapToGPU(hostDevice, targetDevice, *simulation.turbulentField_);
    StencilDelegatePtrs    stencilDelegatePtrs = StencilDelegate::mapToGPU(hostDevice, targetDevice, *simulation.stencil_);
    ParametersGPUPtrs      parametersPtrs      = Parameters::mapToGPU(hostDevice, targetDevice, *simulation.parameters_);

    bool copiedFlowFieldPtr = omp_target_memcpy(&simulationGPU->turbulentField_, &flowFieldPtrs.flowFieldPtr_, sizeof(TurbulentFlowField*), 0, 0, targetDevice, hostDevice) == 0;
    if (!copiedFlowFieldPtr) {
      std::cout << "Failed to copy TurbulentFlowField pointer to TurbulentSimulation object" << std::endl;
    }

    bool copiedStencilDelegatePtr = omp_target_memcpy(&simulationGPU->stencil_, &stencilDelegatePtrs.stencilDelegateGPU_, sizeof(StencilDelegate*), 0, 0, targetDevice, hostDevice)
                                    == 0;
    if (!copiedStencilDelegatePtr) {
      std::cout << "Failed to copy StencilDelegate pointer to TurbulentSimulation object" << std::endl;
    }
    bool copiedParameterPtr = omp_target_memcpy(&simulationGPU->parameters_, &parametersPtrs.parametersGPUPtrs_, sizeof(StencilDelegate*), 0, 0, targetDevice, hostDevice) == 0;
    if (!copiedParameterPtr) {
      std::cout << "Failed to copy Parameter pointer to TurbulentSimulation object" << std::endl;
    }

    omp_target_memcpy(&simulationGPU->fieldIterator_.highOffset_, &simulation.fieldIterator_.highOffset_, sizeof(int), 0, 0, targetDevice, hostDevice);
    omp_target_memcpy(&simulationGPU->fieldIterator_.lowOffset_, &simulation.fieldIterator_.lowOffset_, sizeof(int), 0, 0, targetDevice, hostDevice);
    omp_target_memcpy(&simulationGPU->boundaryIterator_.highOffset_, &simulation.boundaryIterator_.highOffset_, sizeof(int), 0, 0, targetDevice, hostDevice);
    omp_target_memcpy(&simulationGPU->boundaryIterator_.lowOffset_, &simulation.boundaryIterator_.lowOffset_, sizeof(int), 0, 0, targetDevice, hostDevice);
    omp_target_memcpy(&simulationGPU->wallIterator_.highOffset_, &simulation.wallIterator_.highOffset_, sizeof(int), 0, 0, targetDevice, hostDevice);
    omp_target_memcpy(&simulationGPU->wallIterator_.lowOffset_, &simulation.wallIterator_.lowOffset_, sizeof(int), 0, 0, targetDevice, hostDevice);
    omp_target_memcpy(&simulationGPU->turbulentFieldIterator_.highOffset_, &simulation.turbulentFieldIterator_.highOffset_, sizeof(int), 0, 0, targetDevice, hostDevice);
    omp_target_memcpy(&simulationGPU->turbulentFieldIterator_.lowOffset_, &simulation.turbulentFieldIterator_.lowOffset_, sizeof(int), 0, 0, targetDevice, hostDevice);
    omp_target_memcpy(&simulationGPU->turbulentBoundaryIterator_.highOffset_, &simulation.turbulentBoundaryIterator_.highOffset_, sizeof(int), 0, 0, targetDevice, hostDevice);
    omp_target_memcpy(&simulationGPU->turbulentBoundaryIterator_.lowOffset_, &simulation.turbulentBoundaryIterator_.lowOffset_, sizeof(int), 0, 0, targetDevice, hostDevice);

    return TurbulentSimulationPtrs(simulationGPU, flowFieldPtrs, stencilDelegatePtrs, parametersPtrs);
  }

  static void mapToCPU(int hostDevice, int targetDevice, TurbulentSimulation& simulation, TurbulentSimulationPtrs& simulationPtrs) {
    TurbulentFlowField::mapToCPU(hostDevice, targetDevice, *simulation.turbulentField_, simulationPtrs.flowFieldGPUptrs_);
  }

  static void mapToCPUAndFree(int hostDevice, int targetDevice, TurbulentSimulation& simulation, TurbulentSimulationPtrs& simulationPtrs) {
    bool disassociatedSimulation = omp_target_disassociate_ptr(&simulation, targetDevice) == 0;
    if (!disassociatedSimulation) {
      std::cout << "Error: Could not disassociate Simulation pointer." << std::endl;
    }

    TurbulentFlowField::mapToCPUAndFree(hostDevice, targetDevice, *simulation.turbulentField_, simulationPtrs.flowFieldGPUptrs_);
    StencilDelegate::freeGPU(hostDevice, targetDevice, *simulation.stencil_, simulationPtrs.stencilDelegatePtrs_);
    omp_target_memcpy(&simulation.parameters_->timestep.dt, &simulationPtrs.parameterPtrs_.parametersGPUPtrs_->timestep.dt, sizeof(RealType), 0, 0, hostDevice, targetDevice);
    Parameters::freeGPU(hostDevice, targetDevice, *simulation.parameters_, simulationPtrs.parameterPtrs_);
  }
};
