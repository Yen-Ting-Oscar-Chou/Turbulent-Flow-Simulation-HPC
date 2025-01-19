#include "StdAfx.hpp"

#include "MaxUStencil.hpp"

Stencils::MaxUStencil::MaxUStencil() { reset(); }

void Stencils::MaxUStencil::apply(const Parameters& parameters, FlowField& flowField, int i, int j) { cellMaxValue(parameters, flowField, i, j); }

void Stencils::MaxUStencil::apply(const Parameters& parameters, FlowField& flowField, int i, int j, int k) { cellMaxValue(parameters, flowField, i, j, k); }

void Stencils::MaxUStencil::applyLeftWall(const Parameters& parameters, FlowField& flowField, int i, int j) { cellMaxValue(parameters, flowField, i, j); }

void Stencils::MaxUStencil::applyRightWall(const Parameters& parameters, FlowField& flowField, int i, int j) { cellMaxValue(parameters, flowField, i, j); }

void Stencils::MaxUStencil::applyBottomWall(const Parameters& parameters, FlowField& flowField, int i, int j) { cellMaxValue(parameters, flowField, i, j); }

void Stencils::MaxUStencil::applyTopWall(const Parameters& parameters, FlowField& flowField, int i, int j) { cellMaxValue(parameters, flowField, i, j); }

void Stencils::MaxUStencil::applyLeftWall(const Parameters& parameters, FlowField& flowField, int i, int j, int k) { cellMaxValue(parameters, flowField, i, j, k); }

void Stencils::MaxUStencil::applyRightWall(const Parameters& parameters, FlowField& flowField, int i, int j, int k) { cellMaxValue(parameters, flowField, i, j, k); }

void Stencils::MaxUStencil::applyBottomWall(const Parameters& parameters, FlowField& flowField, int i, int j, int k) { cellMaxValue(parameters, flowField, i, j, k); }

void Stencils::MaxUStencil::applyTopWall(const Parameters& parameters, FlowField& flowField, int i, int j, int k) { cellMaxValue(parameters, flowField, i, j, k); }

void Stencils::MaxUStencil::applyFrontWall(const Parameters& parameters, FlowField& flowField, int i, int j, int k) { cellMaxValue(parameters, flowField, i, j, k); }

void Stencils::MaxUStencil::applyBackWall(const Parameters& parameters, FlowField& flowField, int i, int j, int k) { cellMaxValue(parameters, flowField, i, j, k); }

void Stencils::MaxUStencil::cellMaxValue(const Parameters& parameters, FlowField& flowField, int i, int j) {
  const RealType velocityX = flowField.getVelocity().getVectorElement(i, j, 0);
  const RealType velocityY = flowField.getVelocity().getVectorElement(i, j, 1);
  const RealType dx        = parameters.meshsize->getDx(i, j);
  const RealType dy        = parameters.meshsize->getDy(i, j);

  const RealType scaledVelocityX = fabs(velocityX) / dx;
  const RealType scaledVelocityY = fabs(velocityY) / dy;

  if (scaledVelocityX > maxValue_) {
    maxValue_ = scaledVelocityX;
  }
  if (scaledVelocityY > maxValue_) {
    maxValue_ = scaledVelocityY;
  }
}

void Stencils::MaxUStencil::cellMaxValue(const Parameters& parameters, FlowField& flowField, int i, int j, int k) {
  const RealType velocityX = flowField.getVelocity().getVectorElement(i, j, k, 0);
  const RealType velocityY = flowField.getVelocity().getVectorElement(i, j, k, 1);
  const RealType velocityZ = flowField.getVelocity().getVectorElement(i, j, k, 2);
  const RealType dx        = parameters.meshsize->getDx(i, j, k);
  const RealType dy        = parameters.meshsize->getDy(i, j, k);
  const RealType dz        = parameters.meshsize->getDz(i, j, k);

  const RealType scaledVelocityX = fabs(velocityX) / dx;
  const RealType scaledVelocityY = fabs(velocityY) / dy;
  const RealType scaledVelocityZ = fabs(velocityZ) / dz;
  if (scaledVelocityX > maxValue_ || scaledVelocityY > maxValue_ || scaledVelocityZ > maxValue_) {
#pragma omp critical(maxUVelX)
    if (scaledVelocityX > maxValue_) {
      maxValue_ = scaledVelocityX;
    }
#pragma omp critical(maxUVelY)
    if (scaledVelocityY > maxValue_) {
      maxValue_ = scaledVelocityY;
    }
#pragma omp critical(maxUVelZ)
    if (scaledVelocityZ > maxValue_) {
      maxValue_ = scaledVelocityZ;
    }
  }
}

void Stencils::MaxUStencil::reset() { maxValue_ = 0.0; }

const RealType Stencils::MaxUStencil::getMaxValue() const { return maxValue_; }
