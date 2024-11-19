#include "StdAfx.hpp"

#include "TurbulentSimulation.hpp"

TurbulentSimulation::TurbulentSimulation(Parameters& parameters, TurbulentFlowField& flowField):
  Simulation(parameters, flowField),
  turbulentField_(flowField) {};

void TurbulentSimulation::initializeFlowField() {
  Simulation::initializeFlowField();
  std::list<std::tuple<RealType, RealType, RealType>> coordinateList3D;
  std::list<std::tuple<RealType, RealType>>           coordinateList2D;
  Stencils::ObstacleCoordinatesStencil                obstStencil(parameters_, coordinateList2D, coordinateList3D);
  FieldIterator<FlowField>                            obstIterator(flowField_, parameters_, obstStencil);
  obstIterator.iterate();
  if (parameters_.geometry.dim == 2) {
    Stencils::DistanceStencil         distStencil(parameters_, coordinateList2D);
    FieldIterator<TurbulentFlowField> distIterator(turbulentField_, parameters_, distStencil);
    distIterator.iterate();
  } else {
    Stencils::DistanceStencil         distStencil(parameters_, coordinateList3D);
    FieldIterator<TurbulentFlowField> distIterator(turbulentField_, parameters_, distStencil);
    distIterator.iterate();
  }
}

void TurbulentSimulation::plotVTK(int timeStep, RealType simulationTime) {
  Stencils::TurbulentVTKStencil     vtkStencil(parameters_);
  FieldIterator<TurbulentFlowField> vtkIterator(turbulentField_, parameters_, vtkStencil, 1, 0);

  vtkIterator.iterate();
  vtkStencil.write(turbulentField_, timeStep, simulationTime);
}