#include "StdAfx.hpp"

#include "MovingWallStencils.hpp"

void Stencils::MovingWallVelocityStencil::applyLeftWall(const Parameters& parameters, FlowField& flowField, int i, int j) {
  flowField.getVelocity().getVectorElement(i, j, 0) = parameters.walls.vectorLeft[0];
  flowField.getVelocity().getVectorElement(i, j, 1) = 2 * parameters.walls.vectorLeft[1] - flowField.getVelocity().getVectorElement(i + 1, j, 1);
}

void Stencils::MovingWallVelocityStencil::applyRightWall(const Parameters& parameters, FlowField& flowField, int i, int j) {
  flowField.getVelocity().getVectorElement(i - 1, j, 0) = parameters.walls.vectorRight[0];
  flowField.getVelocity().getVectorElement(i, j, 1)     = 2 * parameters.walls.vectorRight[1] - flowField.getVelocity().getVectorElement(i - 1, j, 1);
}

void Stencils::MovingWallVelocityStencil::applyBottomWall(const Parameters& parameters, FlowField& flowField, int i, int j) {
  flowField.getVelocity().getVectorElement(i, j, 0) = 2 * parameters.walls.vectorBottom[0] - flowField.getVelocity().getVectorElement(i, j + 1, 0);
  flowField.getVelocity().getVectorElement(i, j, 1) = parameters.walls.vectorBottom[1];
}

void Stencils::MovingWallVelocityStencil::applyTopWall(const Parameters& parameters, FlowField& flowField, int i, int j) {
  flowField.getVelocity().getVectorElement(i, j, 0)      = 2 * parameters.walls.vectorTop[0] - flowField.getVelocity().getVectorElement(i, j - 1, 0);
  flowField.getVelocity().getVectorElement(i, j - 1, 1) = parameters.walls.vectorTop[1];
}

void Stencils::MovingWallVelocityStencil::applyLeftWall(const Parameters& parameters, FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVectorElement(i, j, k, 0) = parameters.walls.vectorLeft[0];
  flowField.getVelocity().getVectorElement(i, j, k, 1) = 2 * parameters.walls.vectorLeft[1] - flowField.getVelocity().getVectorElement(i + 1, j, k, 1);
  flowField.getVelocity().getVectorElement(i, j, k, 2) = 2 * parameters.walls.vectorLeft[2] - flowField.getVelocity().getVectorElement(i + 1, j, k, 2);
}

void Stencils::MovingWallVelocityStencil::applyRightWall(const Parameters& parameters, FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVectorElement(i - 1, j, k, 0) = parameters.walls.vectorRight[0];
  flowField.getVelocity().getVectorElement(i, j, k, 1)     = 2 * parameters.walls.vectorRight[1] - flowField.getVelocity().getVectorElement(i - 1, j, k, 1);
  flowField.getVelocity().getVectorElement(i, j, k, 2)     = 2 * parameters.walls.vectorRight[2] - flowField.getVelocity().getVectorElement(i - 1, j, k, 2);
}

void Stencils::MovingWallVelocityStencil::applyBottomWall(const Parameters& parameters, FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVectorElement(i, j, k, 0) = 2 * parameters.walls.vectorBottom[0] - flowField.getVelocity().getVectorElement(i, j + 1, k, 0);
  flowField.getVelocity().getVectorElement(i, j, k, 1) = parameters.walls.vectorBottom[1];
  flowField.getVelocity().getVectorElement(i, j, k, 2) = 2 * parameters.walls.vectorBottom[2] - flowField.getVelocity().getVectorElement(i, j + 1, k, 2);
}

void Stencils::MovingWallVelocityStencil::applyTopWall(const Parameters& parameters, FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVectorElement(i, j, k, 0)      = 2 * parameters.walls.vectorTop[0] - flowField.getVelocity().getVectorElement(i, j - 1, k, 0);
  flowField.getVelocity().getVectorElement(i, j - 1, k, 1) = parameters.walls.vectorTop[1];
  flowField.getVelocity().getVectorElement(i, j, k, 2)      = 2 * parameters.walls.vectorTop[2] - flowField.getVelocity().getVectorElement(i, j - 1, k, 2);
}

void Stencils::MovingWallVelocityStencil::applyFrontWall(const Parameters& parameters, FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVectorElement(i, j, k, 0) = 2 * parameters.walls.vectorFront[0] - flowField.getVelocity().getVectorElement(i, j, k + 1, 0);
  flowField.getVelocity().getVectorElement(i, j, k, 1) = 2 * parameters.walls.vectorFront[1] - flowField.getVelocity().getVectorElement(i, j, k + 1, 1);
  flowField.getVelocity().getVectorElement(i, j, k, 2) = parameters.walls.vectorFront[2];
}

void Stencils::MovingWallVelocityStencil::applyBackWall(const Parameters& parameters, FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVectorElement(i, j, k, 0)      = 2 * parameters.walls.vectorBack[0] - flowField.getVelocity().getVectorElement(i, j, k - 1, 0);
  flowField.getVelocity().getVectorElement(i, j, k, 1)      = 2 * parameters.walls.vectorBack[1] - flowField.getVelocity().getVectorElement(i, j, k - 1, 1);
  flowField.getVelocity().getVectorElement(i, j, k - 1, 2) = parameters.walls.vectorBack[2];
}

void Stencils::MovingWallFGHStencil::applyLeftWall(const Parameters& parameters, FlowField& flowField, int i, int j) {
  flowField.getFGH().getVectorElement(i, j, 0) = parameters.walls.vectorLeft[0];
}

void Stencils::MovingWallFGHStencil::applyRightWall(const Parameters& parameters, FlowField& flowField, int i, int j) {
  flowField.getFGH().getVectorElement(i - 1, j, 0) = parameters.walls.vectorRight[0];
}

void Stencils::MovingWallFGHStencil::applyBottomWall(const Parameters& parameters, FlowField& flowField, int i, int j) {
  flowField.getFGH().getVectorElement(i, j, 1) = parameters.walls.vectorBottom[1];
}

void Stencils::MovingWallFGHStencil::applyTopWall(const Parameters& parameters, FlowField& flowField, int i, int j) {
  flowField.getFGH().getVectorElement(i, j - 1, 1) = parameters.walls.vectorTop[1];
}

void Stencils::MovingWallFGHStencil::applyLeftWall(const Parameters& parameters, FlowField& flowField, int i, int j, int k) {
  flowField.getFGH().getVectorElement(i, j, k, 0) = parameters.walls.vectorLeft[0];
}

void Stencils::MovingWallFGHStencil::applyRightWall(const Parameters& parameters, FlowField& flowField, int i, int j, int k) {
  flowField.getFGH().getVectorElement(i - 1, j, k, 0) = parameters.walls.vectorRight[0];
}

void Stencils::MovingWallFGHStencil::applyBottomWall(const Parameters& parameters, FlowField& flowField, int i, int j, int k) {
  flowField.getFGH().getVectorElement(i, j, k, 1) = parameters.walls.vectorBottom[1];
}

void Stencils::MovingWallFGHStencil::applyTopWall(const Parameters& parameters, FlowField& flowField, int i, int j, int k) {
  flowField.getFGH().getVectorElement(i, j - 1, k, 1) = parameters.walls.vectorTop[1];
}

void Stencils::MovingWallFGHStencil::applyFrontWall(const Parameters& parameters, FlowField& flowField, int i, int j, int k) {
  flowField.getFGH().getVectorElement(i, j, k, 2) = parameters.walls.vectorFront[2];
}

void Stencils::MovingWallFGHStencil::applyBackWall(const Parameters& parameters, FlowField& flowField, int i, int j, int k) {
  flowField.getFGH().getVectorElement(i, j, k - 1, 2) = parameters.walls.vectorBack[2];
}
