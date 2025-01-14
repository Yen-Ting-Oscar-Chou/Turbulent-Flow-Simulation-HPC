#include "StdAfx.hpp"

#include "ObstacleCoordinatesStencil.hpp"

Stencils::ObstacleCoordinatesStencil::ObstacleCoordinatesStencil(
  std::vector<std::tuple<RealType, RealType>>& coordinatesList2D, std::vector<std::tuple<RealType, RealType, RealType>>& coordinatesList3D
):
  coordinatesList2D(coordinatesList2D),
  coordinatesList3D(coordinatesList3D) {}

void Stencils::ObstacleCoordinatesStencil::apply(const Parameters& parameters, FlowField& flowField, int i, int j) {
  const int obstacle = flowField.getFlags().getValue(i, j);
  if ((obstacle & OBSTACLE_SELF) == 1) {
    addToList(parameters, i, j);
  }
}

void Stencils::ObstacleCoordinatesStencil::addToList(const Parameters& parameters, int i, int j) {
  RealType coords[2] = {0.0, 0.0};
  computeGlobalCoordinates(coords, parameters, i, j);
  coordinatesList2D.push_back(std::make_tuple(coords[0], coords[1]));
}

void Stencils::ObstacleCoordinatesStencil::apply(const Parameters& parameters, FlowField& flowField, int i, int j, int k) {
  const int obstacle = flowField.getFlags().getValue(i, j, k);
  if ((obstacle & OBSTACLE_SELF) == 1) {
    // TODO only use outer obstacle cells
    addToList(parameters, i, j, k);
  }
}

void Stencils::ObstacleCoordinatesStencil::addToList(const Parameters& parameters, int i, int j, int k) {
  RealType coords[3] = {0.0, 0.0, 0.0};
  computeGlobalCoordinates(coords, parameters, i, j, k);
  coordinatesList3D.push_back(std::make_tuple(coords[0], coords[1], coords[2]));
}

void Stencils::ObstacleCoordinatesStencil::applyLeftWall(const Parameters& parameters, [[maybe_unused]] FlowField& flowField, int i, int j) {
  std::string scenario = parameters.simulation.scenario;
  if (scenario == "cavity") {
    addToList(parameters, i, j);
  } else if (scenario == "channel") {
    addToList(parameters, i, j);
  }
}
void Stencils::ObstacleCoordinatesStencil::applyRightWall(const Parameters& parameters, [[maybe_unused]] FlowField& flowField, int i, int j) {
  std::string scenario = parameters.simulation.scenario;
  if (scenario == "cavity") {
    addToList(parameters, i, j);
  }
}
void Stencils::ObstacleCoordinatesStencil::applyBottomWall(const Parameters& parameters, [[maybe_unused]] FlowField& flowField, int i, int j) {
  std::string scenario = parameters.simulation.scenario;
  if (scenario == "cavity") {
    addToList(parameters, i, j);
  } else if (scenario == "channel") {
    addToList(parameters, i, j);
  }
}
void Stencils::ObstacleCoordinatesStencil::applyTopWall(const Parameters& parameters, [[maybe_unused]] FlowField& flowField, int i, int j) {
  std::string scenario = parameters.simulation.scenario;
  if (scenario == "cavity") {
    addToList(parameters, i, j);
  } else if (scenario == "channel") {
    addToList(parameters, i, j);
  }
}

void Stencils::ObstacleCoordinatesStencil::applyLeftWall(const Parameters& parameters, [[maybe_unused]] FlowField& flowField, int i, int j, int k) {
  std::string scenario = parameters.simulation.scenario;
  if (scenario == "cavity") {
    addToList(parameters, i, j, k);
  } else if (scenario == "channel") {
    addToList(parameters, i, j, k);
  }
}
void Stencils::ObstacleCoordinatesStencil::applyRightWall(const Parameters& parameters, [[maybe_unused]] FlowField& flowField, int i, int j, int k) {
  std::string scenario = parameters.simulation.scenario;
  if (scenario == "cavity") {
    addToList(parameters, i, j, k);
  }
}
void Stencils::ObstacleCoordinatesStencil::applyBottomWall(const Parameters& parameters, [[maybe_unused]] FlowField& flowField, int i, int j, int k) {
  std::string scenario = parameters.simulation.scenario;
  if (scenario == "cavity") {
    addToList(parameters, i, j, k);
  } else if (scenario == "channel") {
    addToList(parameters, i, j, k);
  }
}
void Stencils::ObstacleCoordinatesStencil::applyTopWall(const Parameters& parameters, [[maybe_unused]] FlowField& flowField, int i, int j, int k) {
  std::string scenario = parameters.simulation.scenario;
  if (scenario == "cavity") {
    addToList(parameters, i, j, k);
  } else if (scenario == "channel") {
    addToList(parameters, i, j, k);
  }
}
void Stencils::ObstacleCoordinatesStencil::applyFrontWall(const Parameters& parameters, [[maybe_unused]] FlowField& flowField, int i, int j, int k) {
  std::string scenario = parameters.simulation.scenario;
  if (scenario == "cavity") {
    addToList(parameters, i, j, k);
  } else if (scenario == "channel") {
    addToList(parameters, i, j, k);
  }
}
void Stencils::ObstacleCoordinatesStencil::applyBackWall(const Parameters& parameters, [[maybe_unused]] FlowField& flowField, int i, int j, int k) {
  std::string scenario = parameters.simulation.scenario;
  if (scenario == "cavity") {
    addToList(parameters, i, j, k);
  } else if (scenario == "channel") {
    addToList(parameters, i, j, k);
  }
}