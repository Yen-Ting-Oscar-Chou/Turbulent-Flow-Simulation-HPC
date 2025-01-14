#include "StdAfx.hpp"

#include "InitTaylorGreenFlowFieldStencil.hpp"

Stencils::InitTaylorGreenFlowFieldStencil::InitTaylorGreenFlowFieldStencil():
  pi2_(2.0 * 3.141592653589793238) {}

void Stencils::InitTaylorGreenFlowFieldStencil::apply(const Parameters& parameters, FlowField& flowField, int i, int j) {
  const RealType coordX = parameters.meshsize.getPosX(i, j) + 0.5 * parameters.meshsize.getDx(i, j);
  const RealType coordY = parameters.meshsize.getPosY(i, j) + 0.5 * parameters.meshsize.getDy(i, j);

  // Initialize velocities
  flowField.getVelocity().getVectorElement(
    i, j, 0
  ) = sin(pi2_ * (coordX + 0.5 * parameters.meshsize.getDx(i, j)) / parameters.geometry.lengthX) * sin(pi2_ * coordY / parameters.geometry.lengthY);
  flowField.getVelocity().getVectorElement(
    i, j, 1
  ) = cos(pi2_ * coordX / parameters.geometry.lengthX) * cos(pi2_ * (coordY + 0.5 * parameters.meshsize.getDy(i, j)) / parameters.geometry.lengthY);
}

void Stencils::InitTaylorGreenFlowFieldStencil::apply(const Parameters& parameters, FlowField& flowField, int i, int j, int k) {
  const RealType coordX = parameters.meshsize.getPosX(i, j, k) + 0.5 * parameters.meshsize.getDx(i, j, k);
  const RealType coordY = parameters.meshsize.getPosY(i, j, k) + 0.5 * parameters.meshsize.getDy(i, j, k);
  const RealType coordZ = parameters.meshsize.getPosZ(i, j, k) + 0.5 * parameters.meshsize.getDz(i, j, k);
  
  // Initialize velocities
  flowField.getVelocity().getVectorElement(
    i, j, k, 0
  ) = cos(pi2_ * (coordX + 0.5 * parameters.meshsize.getDx(i, j, k)) / parameters.geometry.lengthX) * sin(pi2_ * coordY / parameters.geometry.lengthY)
      * sin(pi2_ * coordZ / parameters.geometry.lengthZ);
  flowField.getVelocity().getVectorElement(
    i, j, k, 1
  ) = sin(pi2_ * coordX / parameters.geometry.lengthX) * cos(pi2_ * (coordY + 0.5 * parameters.meshsize.getDy(i, j, k)) / parameters.geometry.lengthY)
      * sin(pi2_ * coordZ / parameters.geometry.lengthZ);
  flowField.getVelocity().getVectorElement(
    i, j, k, 2
  ) = sin(pi2_ * coordX / parameters.geometry.lengthX) * sin(pi2_ * coordY / parameters.geometry.lengthY)
      * cos(pi2_ * (coordZ + 0.5 * parameters.meshsize.getDz(i, j, k)) / parameters.geometry.lengthZ);
}
