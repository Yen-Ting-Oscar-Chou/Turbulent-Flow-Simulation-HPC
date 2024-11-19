#pragma once

#include "Simulation.hpp"
#include "Stencils/DistanceStencil.hpp"
#include "Stencils/ObstacleCoordinatesStencil.hpp"
#include "Stencils/TurbulentVTKStencil.hpp"
#include "Stencils/ViscosityStencil.hpp"
#include "TurbulentFlowField.hpp"

class TurbulentSimulation: public Simulation {
private:
  TurbulentFlowField&               turbulentField_;
  Stencils::ViscosityStencil        viscosityStencil_;
  FieldIterator<TurbulentFlowField> viscosityIterator_;

public:
  TurbulentSimulation(Parameters& parameters, TurbulentFlowField& flowField);
  virtual ~TurbulentSimulation() = default;

  /** Initialises the flow field according to the scenario */
  virtual void initializeFlowField() override;

  virtual void solveTimestep() override;

  /** Plots the flow field */
  virtual void plotVTK(int timeStep, RealType simulationTime) override;
};
