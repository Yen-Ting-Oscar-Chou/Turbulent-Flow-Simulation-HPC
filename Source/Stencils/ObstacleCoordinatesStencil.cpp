#include "StdAfx.hpp"

#include "ObstacleCoordinatesStencil.hpp"

Stencils::ObstacleCoordinatesStencil::ObstacleCoordinatesStencil(
  const Parameters& parameters, std::vector<std::tuple<RealType, RealType>>& coordinatesList2D, std::vector<std::tuple<RealType, RealType, RealType>>& coordinatesList3D
):
  FieldStencil<FlowField>(parameters),
  BoundaryStencil<FlowField>(parameters),
  coordinatesList2D(coordinatesList2D),
  coordinatesList3D(coordinatesList3D) {}

void Stencils::ObstacleCoordinatesStencil::apply(FlowField& flowField, int i, int j) {
  const int obstacle = flowField.getFlags().getValue(i, j);
  if ((obstacle & OBSTACLE_SELF) == 1) {
    addToList(i, j);
  }
}

void Stencils::ObstacleCoordinatesStencil::addToList(int i, int j) {
  RealType coords[2] = {0.0, 0.0};
  computeGlobalCoordinates(coords, FieldStencil::parameters_, i, j);
  coordinatesList2D.push_back(std::make_tuple(coords[0], coords[1]));
}

void Stencils::ObstacleCoordinatesStencil::apply(FlowField& flowField, int i, int j, int k) {
  const int obstacle = flowField.getFlags().getValue(i, j, k);
  if ((obstacle & OBSTACLE_SELF) == 1) {
    // TODO only use outer obstacle cells
    addToList(i, j, k);
  }
}

void Stencils::ObstacleCoordinatesStencil::addToList(int i, int j, int k) {
  RealType coords[3] = {0.0, 0.0, 0.0};
  computeGlobalCoordinates(coords, FieldStencil::parameters_, i, j, k);
  coordinatesList3D.push_back(std::make_tuple(coords[0], coords[1], coords[2]));
}

void Stencils::ObstacleCoordinatesStencil::applyLeftWall([[maybe_unused]] FlowField& flowField, int i, int j) {
  std::string scenario = FieldStencil::parameters_.simulation.scenario;
  if (scenario == "cavity") {
    addToList(i, j);
  } else if (scenario == "channel") {
    addToList(i, j);
  }
}
void Stencils::ObstacleCoordinatesStencil::applyRightWall([[maybe_unused]] FlowField& flowField, int i, int j) {
  std::string scenario = FieldStencil::parameters_.simulation.scenario;
  if (scenario == "cavity") {
    addToList(i, j);
  }
};
void Stencils::ObstacleCoordinatesStencil::applyBottomWall([[maybe_unused]] FlowField& flowField, int i, int j) {
  std::string scenario = FieldStencil::parameters_.simulation.scenario;
  if (scenario == "cavity") {
    addToList(i, j);
  } else if (scenario == "channel") {
    addToList(i, j);
  }
};
void Stencils::ObstacleCoordinatesStencil::applyTopWall([[maybe_unused]] FlowField& flowField, int i, int j) {
  std::string scenario = FieldStencil::parameters_.simulation.scenario;
  if (scenario == "cavity") {
    addToList(i, j);
  } else if (scenario == "channel") {
    addToList(i, j);
  }
};

void Stencils::ObstacleCoordinatesStencil::applyLeftWall([[maybe_unused]] FlowField& flowField, int i, int j, int k) {
  std::string scenario = FieldStencil::parameters_.simulation.scenario;
  if (scenario == "cavity") {
    addToList(i, j, k);
  } else if (scenario == "channel") {
    addToList(i, j, k);
  }
}
void Stencils::ObstacleCoordinatesStencil::applyRightWall([[maybe_unused]] FlowField& flowField, int i, int j, int k) {
  std::string scenario = FieldStencil::parameters_.simulation.scenario;
  if (scenario == "cavity") {
    addToList(i, j, k);
  }
};
void Stencils::ObstacleCoordinatesStencil::applyBottomWall([[maybe_unused]] FlowField& flowField, int i, int j, int k) {
  std::string scenario = FieldStencil::parameters_.simulation.scenario;
  if (scenario == "cavity") {
    addToList(i, j, k);
  } else if (scenario == "channel") {
    addToList(i, j);
  }
};
void Stencils::ObstacleCoordinatesStencil::applyTopWall([[maybe_unused]] FlowField& flowField, int i, int j, int k) {
  std::string scenario = FieldStencil::parameters_.simulation.scenario;
  if (scenario == "cavity") {
    addToList(i, j, k);
  } else if (scenario == "channel") {
    addToList(i, j, k);
  }
};
void Stencils::ObstacleCoordinatesStencil::applyFrontWall([[maybe_unused]] FlowField& flowField, int i, int j, int k) {
  std::string scenario = FieldStencil::parameters_.simulation.scenario;
  if (scenario == "cavity") {
    addToList(i, j, k);
  } else if (scenario == "channel") {
    addToList(i, j, k);
  }
};
void Stencils::ObstacleCoordinatesStencil::applyBackWall([[maybe_unused]] FlowField& flowField, int i, int j, int k) {
  std::string scenario = FieldStencil::parameters_.simulation.scenario;
  if (scenario == "cavity") {
    addToList(i, j, k);
  } else if (scenario == "channel") {
    addToList(i, j, k);
  }
};