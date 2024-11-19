#include "StdAfx.hpp"

#include "ObstacleCoordinatesStencil.hpp"

Stencils::ObstacleCoordinatesStencil::ObstacleCoordinatesStencil(const Parameters& parameters, std::list<std::tuple<RealType, RealType>>& coordinatesList2D, std::list<std::tuple<RealType, RealType, RealType>>& coordinatesList3D):
  FieldStencil<FlowField>(parameters),
  coordinatesList2D(coordinatesList2D),
  coordinatesList3D(coordinatesList3D) {}

void Stencils::ObstacleCoordinatesStencil::apply(FlowField& flowField, int i, int j) {
  const int obstacle = flowField.getFlags().getValue(i, j);
  if ((obstacle & OBSTACLE_SELF) == 1) {
    //TODO only use outer obstacle cells
    RealType coords[2] = {0.0, 0.0};
    computeGlobalCoordinates(coords, parameters_, i, j);
    coordinatesList2D.push_back(std::make_tuple(coords[0], coords[1]));
  }
}

void Stencils::ObstacleCoordinatesStencil::apply(FlowField& flowField, int i, int j, int k) {
  const int obstacle = flowField.getFlags().getValue(i, j, k);
  if ((obstacle & OBSTACLE_SELF) == 1) {
    //TODO only use outer obstacle cells
    RealType coords[3] = {0.0, 0.0, 0.0};
    computeGlobalCoordinates(coords, parameters_, i, j, k);
    coordinatesList3D.push_back(std::make_tuple(coords[0], coords[1], coords[2]));
  }
}
