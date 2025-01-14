#include "StdAfx.hpp"

#include "PeriodicBoundaryStencils.hpp"

void Stencils::PeriodicBoundaryVelocityStencil::applyLeftWall([[maybe_unused]] const Parameters& parameters, FlowField& flowField, [[maybe_unused]] int i, int j) {
  flowField.getVelocity().getVectorElement(0, j, 0) = flowField.getVelocity().getVectorElement(flowField.getNx(), j, 0);
  flowField.getVelocity().getVectorElement(1, j, 1) = flowField.getVelocity().getVectorElement(flowField.getNx() + 1, j, 1);
}

void Stencils::PeriodicBoundaryVelocityStencil::applyRightWall([[maybe_unused]] const Parameters& parameters, FlowField& flowField, [[maybe_unused]] int i, int j) {
  flowField.getVelocity().getVectorElement(flowField.getNx() + 2, j, 0) = flowField.getVelocity().getVectorElement(2, j, 0);
  flowField.getVelocity().getVectorElement(flowField.getNx() + 2, j, 1) = flowField.getVelocity().getVectorElement(2, j, 1);
}

void Stencils::PeriodicBoundaryVelocityStencil::applyBottomWall([[maybe_unused]] const Parameters& parameters, FlowField& flowField, int i, [[maybe_unused]] int j) {
  flowField.getVelocity().getVectorElement(i, 0, 1) = flowField.getVelocity().getVectorElement(i, flowField.getNy(), 1);
  flowField.getVelocity().getVectorElement(i, 1, 0) = flowField.getVelocity().getVectorElement(i, flowField.getNy() + 1, 0);
}

void Stencils::PeriodicBoundaryVelocityStencil::applyTopWall([[maybe_unused]] const Parameters& parameters, FlowField& flowField, int i, [[maybe_unused]] int j) {
  flowField.getVelocity().getVectorElement(i, flowField.getNy() + 2, 1) = flowField.getVelocity().getVectorElement(i, 2, 1);
  flowField.getVelocity().getVectorElement(i, flowField.getNy() + 2, 0) = flowField.getVelocity().getVectorElement(i, 2, 0);
}

void Stencils::PeriodicBoundaryVelocityStencil::applyLeftWall([[maybe_unused]] const Parameters& parameters, FlowField& flowField, [[maybe_unused]] int i, int j, int k) {
  flowField.getVelocity().getVectorElement(0, j, k, 0) = flowField.getVelocity().getVectorElement(flowField.getNx(), j, k, 0);
  flowField.getVelocity().getVectorElement(1, j, k, 1) = flowField.getVelocity().getVectorElement(flowField.getNx() + 1, j, k, 1);
  flowField.getVelocity().getVectorElement(1, j, k, 2) = flowField.getVelocity().getVectorElement(flowField.getNx() + 1, j, k, 2);
}

void Stencils::PeriodicBoundaryVelocityStencil::applyRightWall([[maybe_unused]] const Parameters& parameters, FlowField& flowField, [[maybe_unused]] int i, int j, int k) {
  flowField.getVelocity().getVectorElement(flowField.getNx() + 2, j, k, 0) = flowField.getVelocity().getVectorElement(2, j, k, 0);
  flowField.getVelocity().getVectorElement(flowField.getNx() + 2, j, k, 1) = flowField.getVelocity().getVectorElement(2, j, k, 1);
  flowField.getVelocity().getVectorElement(flowField.getNx() + 2, j, k, 2) = flowField.getVelocity().getVectorElement(2, j, k, 2);
}

void Stencils::PeriodicBoundaryVelocityStencil::applyBottomWall([[maybe_unused]] const Parameters& parameters, FlowField& flowField, int i, [[maybe_unused]] int j, int k) {
  flowField.getVelocity().getVectorElement(i, 1, k, 0) = flowField.getVelocity().getVectorElement(i, flowField.getNy() + 1, k, 0);
  flowField.getVelocity().getVectorElement(i, 0, k, 1) = flowField.getVelocity().getVectorElement(i, flowField.getNy(), k, 1);
  flowField.getVelocity().getVectorElement(i, 1, k, 2) = flowField.getVelocity().getVectorElement(i, flowField.getNy() + 1, k, 2);
}

void Stencils::PeriodicBoundaryVelocityStencil::applyTopWall([[maybe_unused]] const Parameters& parameters, FlowField& flowField, int i, [[maybe_unused]] int j, int k) {
  flowField.getVelocity().getVectorElement(i, flowField.getNy() + 2, k, 0) = flowField.getVelocity().getVectorElement(i, 2, k, 0);
  flowField.getVelocity().getVectorElement(i, flowField.getNy() + 2, k, 1) = flowField.getVelocity().getVectorElement(i, 2, k, 1);
  flowField.getVelocity().getVectorElement(i, flowField.getNy() + 2, k, 2) = flowField.getVelocity().getVectorElement(i, 2, k, 2);
}

void Stencils::PeriodicBoundaryVelocityStencil::applyFrontWall([[maybe_unused]] const Parameters& parameters, FlowField& flowField, int i, int j, [[maybe_unused]] int k) {
  flowField.getVelocity().getVectorElement(i, j, 1, 0) = flowField.getVelocity().getVectorElement(i, j, flowField.getNz() + 1, 0);
  flowField.getVelocity().getVectorElement(i, j, 1, 1) = flowField.getVelocity().getVectorElement(i, j, flowField.getNz() + 1, 1);
  flowField.getVelocity().getVectorElement(i, j, 0, 2) = flowField.getVelocity().getVectorElement(i, j, flowField.getNz(), 2);
}

void Stencils::PeriodicBoundaryVelocityStencil::applyBackWall([[maybe_unused]] const Parameters& parameters, FlowField& flowField, int i, int j, [[maybe_unused]] int k) {
  flowField.getVelocity().getVectorElement(i, j, flowField.getNz() + 2, 0) = flowField.getVelocity().getVectorElement(i, j, 2, 0);
  flowField.getVelocity().getVectorElement(i, j, flowField.getNz() + 2, 1) = flowField.getVelocity().getVectorElement(i, j, 2, 1);
  flowField.getVelocity().getVectorElement(i, j, flowField.getNz() + 2, 2) = flowField.getVelocity().getVectorElement(i, j, 2, 2);
}

void Stencils::PeriodicBoundaryFGHStencil::applyLeftWall(
  [[maybe_unused]] const Parameters& parameters, [[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j
) {}
void Stencils::PeriodicBoundaryFGHStencil::applyRightWall(
  [[maybe_unused]] const Parameters& parameters, [[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j
) {}
void Stencils::PeriodicBoundaryFGHStencil::applyBottomWall(
  [[maybe_unused]] const Parameters& parameters, [[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j
) {}
void Stencils::PeriodicBoundaryFGHStencil::applyTopWall(
  [[maybe_unused]] const Parameters& parameters, [[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j
) {}

void Stencils::PeriodicBoundaryFGHStencil::applyLeftWall(
  [[maybe_unused]] const Parameters& parameters, [[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k
) {}
void Stencils::PeriodicBoundaryFGHStencil::applyRightWall(
  [[maybe_unused]] const Parameters& parameters, [[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k
) {}
void Stencils::PeriodicBoundaryFGHStencil::applyBottomWall(
  [[maybe_unused]] const Parameters& parameters, [[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k
) {}
void Stencils::PeriodicBoundaryFGHStencil::applyTopWall(
  [[maybe_unused]] const Parameters& parameters, [[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k
) {}
void Stencils::PeriodicBoundaryFGHStencil::applyFrontWall(
  [[maybe_unused]] const Parameters& parameters, [[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k
) {}
void Stencils::PeriodicBoundaryFGHStencil::applyBackWall(
  [[maybe_unused]] const Parameters& parameters, [[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k
) {}
