#include "StdAfx.hpp"

#include "ObstacleStencil.hpp"

void Stencils::ObstacleStencil::apply(const Parameters& parameters, FlowField& flowField, int i, int j) {
  const int    obstacle = flowField.getFlags().getValue(i, j);
  VectorField& velocity = flowField.getVelocity();

  // Check if current cell is obstacle cell
  if ((obstacle & OBSTACLE_SELF) == 1) {
    // If top cell is fluid, then the no-slip boundary has to be enforced
    if ((obstacle & OBSTACLE_TOP) == 0) {
      const RealType dy_t                = parameters.meshsize->getDy(i, j + 1);
      const RealType dy                  = parameters.meshsize->getDy(i, j);
      velocity.getVectorElement(i, j, 0) = -dy / dy_t * velocity.getVectorElement(i, j + 1, 0);
    }
    // Same for bottom
    if ((obstacle & OBSTACLE_BOTTOM) == 0) {
      const RealType dy_b                = parameters.meshsize->getDy(i, j - 1);
      const RealType dy                  = parameters.meshsize->getDy(i, j);
      velocity.getVectorElement(i, j, 0) = -dy / dy_b * velocity.getVectorElement(i, j - 1, 0);
    }
    // If right cell is fluid, then the no-slip boundary has to be enforced
    if ((obstacle & OBSTACLE_RIGHT) == 0) {
      const RealType dx_r                = parameters.meshsize->getDx(i + 1, j);
      const RealType dx                  = parameters.meshsize->getDx(i, j);
      velocity.getVectorElement(i, j, 1) = -dx / dx_r * velocity.getVectorElement(i + 1, j, 1);
    }
    // Same for left
    if ((obstacle & OBSTACLE_LEFT) == 0) {
      const RealType dx_l                = parameters.meshsize->getDx(i - 1, j);
      const RealType dx                  = parameters.meshsize->getDx(i, j);
      velocity.getVectorElement(i, j, 1) = -dx / dx_l * velocity.getVectorElement(i - 1, j, 1);
    }

    // Set normal velocity to zero if right neighbour is not obstacle
    if ((obstacle & OBSTACLE_RIGHT) == 0) {
      velocity.getVectorElement(i, j, 0) = 0.0;
    }

    // Set normal velocity to zero if top neighbour is not obstacle
    if ((obstacle & OBSTACLE_TOP) == 0) {
      velocity.getVectorElement(i, j, 1) = 0.0;
    }
  }
}

void Stencils::ObstacleStencil::apply(const Parameters& parameters, FlowField& flowField, int i, int j, int k) {
  const int    obstacle = flowField.getFlags().getValue(i, j);
  VectorField& velocity = flowField.getVelocity();

  // Check if current cell is obstacle cell
  if ((obstacle & OBSTACLE_SELF) == 1) {
    // If top cell is fluid: two velocities have to be set: direction 0 and 2.
    if ((obstacle & OBSTACLE_TOP) == 0) {
      const RealType dy_t                   = parameters.meshsize->getDy(i, j + 1, k);
      const RealType dy                     = parameters.meshsize->getDy(i, j, k);
      velocity.getVectorElement(i, j, k, 0) = -dy / dy_t * velocity.getVectorElement(i, j + 1, k, 0);
      velocity.getVectorElement(i, j, k, 2) = -dy / dy_t * velocity.getVectorElement(i, j + 1, k, 2);
    }
    if ((obstacle & OBSTACLE_BOTTOM) == 0) {
      const RealType dy_b                   = parameters.meshsize->getDy(i, j - 1, k);
      const RealType dy                     = parameters.meshsize->getDy(i, j, k);
      velocity.getVectorElement(i, j, k, 0) = -dy / dy_b * velocity.getVectorElement(i, j - 1, k, 0);
      velocity.getVectorElement(i, j, k, 2) = -dy / dy_b * velocity.getVectorElement(i, j - 1, k, 2);
    }

    // If right cell is fluid: two velocities have to be set: direction 1 and 2.
    if ((obstacle & OBSTACLE_RIGHT) == 0) {
      const RealType dx_r                   = parameters.meshsize->getDx(i + 1, j, k);
      const RealType dx                     = parameters.meshsize->getDx(i, j, k);
      velocity.getVectorElement(i, j, k, 1) = -dx / dx_r * velocity.getVectorElement(i + 1, j, k, 1);
      velocity.getVectorElement(i, j, k, 2) = -dx / dx_r * velocity.getVectorElement(i + 1, j, k, 2);
    }
    if ((obstacle & OBSTACLE_LEFT) == 0) {
      const RealType dx_l                   = parameters.meshsize->getDx(i - 1, j, k);
      const RealType dx                     = parameters.meshsize->getDx(i, j, k);
      velocity.getVectorElement(i, j, k, 1) = -dx / dx_l * velocity.getVectorElement(i - 1, j, k, 1);
      velocity.getVectorElement(i, j, k, 2) = -dx / dx_l * velocity.getVectorElement(i - 1, j, k, 2);
    }

    // Same for fluid cell in front
    if ((obstacle & OBSTACLE_BACK) == 0) {
      const RealType dz_f                   = parameters.meshsize->getDx(i, j, k + 1);
      const RealType dz                     = parameters.meshsize->getDx(i, j, k);
      velocity.getVectorElement(i, j, k, 1) = -dz / dz_f * velocity.getVectorElement(i, j, k + 1, 1);
      velocity.getVectorElement(i, j, k, 0) = -dz / dz_f * velocity.getVectorElement(i, j, k + 1, 0);
    }
    if ((obstacle & OBSTACLE_FRONT) == 0) {
      const RealType dz_b                   = parameters.meshsize->getDx(i, j, k - 1);
      const RealType dz                     = parameters.meshsize->getDx(i, j, k);
      velocity.getVectorElement(i, j, k, 1) = -dz / dz_b * velocity.getVectorElement(i, j, k - 1, 1);
      velocity.getVectorElement(i, j, k, 0) = -dz / dz_b * velocity.getVectorElement(i, j, k - 1, 0);
    }

    // Now the normal velocities need to be set to zero to ensure no flow at interfaces between solid and fluid.
    if ((obstacle & OBSTACLE_RIGHT) == 0) {
      velocity.getVectorElement(i, j, k, 0) = 0.0;
    }
    if ((obstacle & OBSTACLE_TOP) == 0) {
      velocity.getVectorElement(i, j, k, 1) = 0.0;
    }
    if ((obstacle & OBSTACLE_BACK) == 0) {
      velocity.getVectorElement(i, j, k, 2) = 0.0;
    }
  }
}
