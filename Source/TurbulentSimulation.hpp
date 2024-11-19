#include "Stencils/DistanceStencil.hpp"
#include "Stencils/ObstacleCoordinatesStencil.hpp"
#include "Simulation.hpp"
#include "TurbulentFlowField.hpp"
#include "Stencils/TurbulentVTKStencil.hpp"

class TurbulentSimulation : public Simulation {
  private:
    TurbulentFlowField& turbulentField_;

public:
  TurbulentSimulation(Parameters& parameters, TurbulentFlowField& flowField);
  virtual ~TurbulentSimulation() = default;

  /** Initialises the flow field according to the scenario */
  virtual void initializeFlowField() override;

  //virtual void solveTimestep() override;

  /** Plots the flow field */
  virtual void plotVTK(int timeStep, RealType simulationTime) override;
};
