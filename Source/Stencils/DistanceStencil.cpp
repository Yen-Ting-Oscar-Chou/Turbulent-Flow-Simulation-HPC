#include "StdAfx.hpp"

#include "DistanceStencil.hpp"

Stencils::DistanceStencil::DistanceStencil(std::vector<std::tuple<RealType, RealType>>& coordinatesList):
  coordinatesList2D(coordinatesList) {}

Stencils::DistanceStencil::DistanceStencil(std::vector<std::tuple<RealType, RealType, RealType>>& coordinatesList):
  coordinatesList3D(coordinatesList) {}

void Stencils::DistanceStencil::apply(const Parameters& parameters, TurbulentFlowField& turbulentField, int i, int j) {
  const int obstacle = turbulentField.getFlags().getValue(i, j);
  if ((obstacle & OBSTACLE_SELF) == 1) {
    return;
  }
  RealType coords[2] = {0.0, 0.0};
  computeGlobalCoordinates(coords, parameters, i, j);
  RealType minDistance = MY_FLOAT_MAX;
  for (const auto& coordsObst : coordinatesList2D) {
    RealType coordsObstArr[2] = {std::get<0>(coordsObst), std::get<1>(coordsObst)};
    minDistance               = std::min(minDistance, computeDistance(coords, coordsObstArr, parameters));
  }
  RealType& dist = turbulentField.getDistance().getScalar(i, j);
  dist           = minDistance;
}

void Stencils::DistanceStencil::apply(const Parameters& parameters, TurbulentFlowField& turbulentField, int i, int j, int k) {
  const int obstacle = turbulentField.getFlags().getValue(i, j, k);
  if ((obstacle & OBSTACLE_SELF) == 1) {
    return;
  }
  RealType coords[3] = {0.0, 0.0, 0.0};
  computeGlobalCoordinates(coords, parameters, i, j, k);
  RealType minDistance = MY_FLOAT_MAX;
  for (const auto& coordsObst : coordinatesList3D) {
    RealType coordsObstArr[3] = {std::get<0>(coordsObst), std::get<1>(coordsObst), std::get<2>(coordsObst)};
    minDistance               = std::min(minDistance, computeDistance(coords, coordsObstArr, parameters));
  }
  RealType& dist = turbulentField.getDistance().getScalar(i, j, k);
  dist           = minDistance;
}