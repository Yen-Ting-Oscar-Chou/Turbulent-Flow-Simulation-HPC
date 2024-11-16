#include "StdAfx.hpp"

#include "DistanceStencil.hpp"

Stencils::DistanceStencil::DistanceStencil(const Parameters& parameters, std::list<std::tuple<RealType, RealType>>& coordinatesList):
  FieldStencil<FlowField>(parameters),
  coordinatesList2D(coordinatesList) {}

Stencils::DistanceStencil::DistanceStencil(const Parameters& parameters, std::list<std::tuple<RealType, RealType, RealType>>& coordinatesList):
  FieldStencil<FlowField>(parameters),
  coordinatesList3D(coordinatesList) {}

void Stencils::DistanceStencil::apply(TurbulentFlowField& turbulentField, int i, int j) {
  const int obstacle = turbulentField.getFlags().getValue(i, j);
  if ((obstacle & OBSTACLE_SELF) == 0) {
    RealType coords[2] = {0.0, 0.0};
    computeGlobalCoordinates(coords, parameters_, i, j);
    RealType minDistance = MY_FLOAT_MAX;
    for (auto const& coordsObst : coordinatesList2D) {
        RealType coordsObstArr[2] = {std::get<0>(coordsObst), std::get<1>(coordsObst)};
        minDistance = std::min(minDistance, computeDistance(coords, coordsObstArr, parameters_));
    }
  }
}

void Stencils::DistanceStencil::apply(TurbulentFlowField& turbulentField, int i, int j, int k) {
  const int obstacle = turbulentField.getFlags().getValue(i, j, k);
  if ((obstacle & OBSTACLE_SELF) == 0) {
    RealType coords[3] = {0.0, 0.0, 0.0};
    computeGlobalCoordinates(coords, parameters_, i, j, k);
    RealType minDistance = MY_FLOAT_MAX;
    for (auto const& coordsObst : coordinatesList3D) {
        RealType coordsObstArr[3] = {std::get<0>(coordsObst), std::get<1>(coordsObst), std::get<2>(coordsObst)};
        minDistance = std::min(minDistance, computeDistance(coords, coordsObstArr, parameters_));
    }
  }
}
