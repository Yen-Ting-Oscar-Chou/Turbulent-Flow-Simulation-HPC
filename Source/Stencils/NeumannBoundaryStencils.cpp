#include "StdAfx.hpp"

#include "NeumannBoundaryStencils.hpp"

void Stencils::NeumannVelocityBoundaryStencil::applyLeftWall(const Parameters& parameters, FlowField& flowField, int i, int j) {
  flowField.getVelocity().getVectorElement(i - 1, j, 0) = flowField.getVelocity().getVectorElement(i, j, 0);
  flowField.getVelocity().getVectorElement(i, j, 1)     = flowField.getVelocity().getVectorElement(i + 1, j, 1);
}

void Stencils::NeumannVelocityBoundaryStencil::applyRightWall(const Parameters& parameters, FlowField& flowField, int i, int j) {
  flowField.getVelocity().getVectorElement(i, j, 0) = flowField.getVelocity().getVectorElement(i - 1, j, 0);
  flowField.getVelocity().getVectorElement(i, j, 1) = flowField.getVelocity().getVectorElement(i - 1, j, 1);
}

void Stencils::NeumannVelocityBoundaryStencil::applyBottomWall(const Parameters& parameters, FlowField& flowField, int i, int j) {
  flowField.getVelocity().getVectorElement(i, j, 0)     = flowField.getVelocity().getVectorElement(i, j + 1, 0);
  flowField.getVelocity().getVectorElement(i, j - 1, 1) = flowField.getVelocity().getVectorElement(i, j, 1);
}

void Stencils::NeumannVelocityBoundaryStencil::applyTopWall(const Parameters& parameters, FlowField& flowField, int i, int j) {
  flowField.getVelocity().getVectorElement(i, j, 0) = flowField.getVelocity().getVectorElement(i, j - 1, 0);
  flowField.getVelocity().getVectorElement(i, j, 1) = flowField.getVelocity().getVectorElement(i, j - 1, 1);
}

void Stencils::NeumannVelocityBoundaryStencil::applyLeftWall(const Parameters& parameters, FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVectorElement(i - 1, j, k, 0) = flowField.getVelocity().getVectorElement(i, j, k, 0);
  flowField.getVelocity().getVectorElement(i, j, k, 1)     = flowField.getVelocity().getVectorElement(i + 1, j, k, 1);
  flowField.getVelocity().getVectorElement(i, j, k, 2)     = flowField.getVelocity().getVectorElement(i + 1, j, k, 2);
}

void Stencils::NeumannVelocityBoundaryStencil::applyRightWall(const Parameters& parameters, FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVectorElement(i, j, k, 0) = flowField.getVelocity().getVectorElement(i - 1, j, k, 0);
  flowField.getVelocity().getVectorElement(i, j, k, 1) = flowField.getVelocity().getVectorElement(i - 1, j, k, 1);
  flowField.getVelocity().getVectorElement(i, j, k, 2) = flowField.getVelocity().getVectorElement(i - 1, j, k, 2);
}

void Stencils::NeumannVelocityBoundaryStencil::applyBottomWall(const Parameters& parameters, FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVectorElement(i, j, k, 0)     = flowField.getVelocity().getVectorElement(i, j + 1, k, 0);
  flowField.getVelocity().getVectorElement(i, j - 1, k, 1) = flowField.getVelocity().getVectorElement(i, j, k, 1);
  flowField.getVelocity().getVectorElement(i, j, k, 2)     = flowField.getVelocity().getVectorElement(i, j + 1, k, 2);
}

void Stencils::NeumannVelocityBoundaryStencil::applyTopWall(const Parameters& parameters, FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVectorElement(i, j, k, 0) = flowField.getVelocity().getVectorElement(i, j - 1, k, 0);
  flowField.getVelocity().getVectorElement(i, j, k, 1) = flowField.getVelocity().getVectorElement(i, j - 1, k, 1);
  flowField.getVelocity().getVectorElement(i, j, k, 2) = flowField.getVelocity().getVectorElement(i, j - 1, k, 2);
}

void Stencils::NeumannVelocityBoundaryStencil::applyFrontWall(const Parameters& parameters, FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVectorElement(i, j, k, 0)     = flowField.getVelocity().getVectorElement(i, j, k + 1, 0);
  flowField.getVelocity().getVectorElement(i, j, k, 1)     = flowField.getVelocity().getVectorElement(i, j, k + 1, 1);
  flowField.getVelocity().getVectorElement(i, j, k - 1, 2) = flowField.getVelocity().getVectorElement(i, j, k, 2);
}

void Stencils::NeumannVelocityBoundaryStencil::applyBackWall(const Parameters& parameters, FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVectorElement(i, j, k, 0) = flowField.getVelocity().getVectorElement(i, j, k - 1, 0);
  flowField.getVelocity().getVectorElement(i, j, k, 1) = flowField.getVelocity().getVectorElement(i, j, k - 1, 1);
  flowField.getVelocity().getVectorElement(i, j, k, 2) = flowField.getVelocity().getVectorElement(i, j, k - 1, 2);
}

// These are left empty. The right values should be computed by the FGH body stencil.
void Stencils::NeumannFGHBoundaryStencil::applyLeftWall(
  [[maybe_unused]] const Parameters& parameters, [[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j
) {}
void Stencils::NeumannFGHBoundaryStencil::applyRightWall(
  [[maybe_unused]] const Parameters& parameters, [[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j
) {}
void Stencils::NeumannFGHBoundaryStencil::applyBottomWall(
  [[maybe_unused]] const Parameters& parameters, [[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j
) {}
void Stencils::NeumannFGHBoundaryStencil::applyTopWall(
  [[maybe_unused]] const Parameters& parameters, [[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j
) {}

void Stencils::NeumannFGHBoundaryStencil::applyLeftWall(
  [[maybe_unused]] const Parameters& parameters, [[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k
) {}
void Stencils::NeumannFGHBoundaryStencil::applyRightWall(
  [[maybe_unused]] const Parameters& parameters, [[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k
) {}
void Stencils::NeumannFGHBoundaryStencil::applyBottomWall(
  [[maybe_unused]] const Parameters& parameters, [[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k
) {}
void Stencils::NeumannFGHBoundaryStencil::applyTopWall(
  [[maybe_unused]] const Parameters& parameters, [[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k
) {}
void Stencils::NeumannFGHBoundaryStencil::applyFrontWall(
  [[maybe_unused]] const Parameters& parameters, [[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k
) {}
void Stencils::NeumannFGHBoundaryStencil::applyBackWall(
  [[maybe_unused]] const Parameters& parameters, [[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k
) {}
