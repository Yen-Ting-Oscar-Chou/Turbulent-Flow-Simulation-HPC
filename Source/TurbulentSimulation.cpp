#include "StdAfx.hpp"

#include "TurbulentSimulation.hpp"

TurbulentSimulation::TurbulentSimulation(Parameters& parameters, TurbulentFlowField& flowField):
  Simulation(parameters, flowField),
  turbulentField_(flowField) {};

void TurbulentSimulation::initializeFlowField() {
  Simulation::initializeFlowField();
  if (parameters_.geometry.dim == 2) {
    std::list<std::tuple<RealType, RealType>> coordinateList;
    Stencils::ObstacleCoordinatesStencil      obstStencil(parameters_, coordinateList);
    FieldIterator<FlowField>                  obstIterator(flowField_, parameters_, obstStencil);
    obstIterator.iterate();

    Stencils::DistanceStencil distStencil(parameters_, coordinateList);
    FieldIterator<TurbulentFlowField> distIterator(turbulentField_, parameters_, distStencil);
    distIterator.iterate();
  } else {
    std::list<std::tuple<RealType, RealType, RealType>> coordinateList;
    Stencils::ObstacleCoordinatesStencil      obstStencil(parameters_, coordinateList);
    FieldIterator<FlowField>                  obstIterator(flowField_, parameters_, obstStencil);
    obstIterator.iterate();

    Stencils::DistanceStencil distStencil(parameters_, coordinateList);
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