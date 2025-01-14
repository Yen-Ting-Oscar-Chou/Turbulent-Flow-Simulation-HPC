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

class TurbulentSimulation: public Simulation {
private:
  TurbulentFlowField&                             turbulentField_;
  Stencils::ViscosityStencil                      viscosityStencil_;
  FieldIterator<TurbulentFlowField>               viscosityIterator_;
  Stencils::TurbulentFGHStencil                   turbulentFGHStencil_;
  FieldIterator<TurbulentFlowField>               turbulentFGHIterator_;
  Stencils::MaxViscStencil                        maxViscStencil_;
  FieldIterator<TurbulentFlowField>               maxViscFieldIterator_;
  GlobalBoundaryIterator<TurbulentFlowField>      maxViscBoundaryIterator_;
  ParallelManagers::TurbulentPetscParallelManager turbulentPetscParallelManager_;

protected:
  void setTimeStep() override;

public:
  TurbulentSimulation(Parameters& parameters, TurbulentFlowField& flowField);
  virtual ~TurbulentSimulation() = default;

  /** Initialises the flow field according to the scenario */
  virtual void initializeFlowField() override;

  virtual void solveTimestep() override;

  /** Plots the flow field */
  virtual void plotVTK(int timeStep, RealType simulationTime) override;
};
