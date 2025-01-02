#include "StdAfx.hpp"

#include "MaxViscStencil.hpp"

Stencils::MaxViscStencil::MaxViscStencil(const Parameters& parameters):
  BoundaryStencil<TurbulentFlowField>(parameters) {
  reset();
}

void Stencils::MaxViscStencil::apply(const Parameters& parameters, TurbulentFlowField& flowField, int i, int j) {
  const RealType visc = flowField.getViscosity().getScalar(i, j);
  maxValue_           = std::max(maxValue_, visc);
}

void Stencils::MaxViscStencil::apply(const Parameters& parameters, TurbulentFlowField& flowField, int i, int j, int k) {
  const RealType visc = flowField.getViscosity().getScalar(i, j, k);
  maxValue_           = std::max(maxValue_, visc);
}

void Stencils::MaxViscStencil::applyLeftWall(TurbulentFlowField& flowField, int i, int j) {
  const RealType visc = flowField.getViscosity().getScalar(i, j);
  maxValue_           = std::max(maxValue_, visc);
}
void Stencils::MaxViscStencil::applyRightWall(TurbulentFlowField& flowField, int i, int j) {
  const RealType visc = flowField.getViscosity().getScalar(i, j);
  maxValue_           = std::max(maxValue_, visc);
}
void Stencils::MaxViscStencil::applyBottomWall(TurbulentFlowField& flowField, int i, int j) {
  const RealType visc = flowField.getViscosity().getScalar(i, j);
  maxValue_           = std::max(maxValue_, visc);
}
void Stencils::MaxViscStencil::applyTopWall(TurbulentFlowField& flowField, int i, int j) {
  const RealType visc = flowField.getViscosity().getScalar(i, j);
  maxValue_           = std::max(maxValue_, visc);
}

void Stencils::MaxViscStencil::applyLeftWall(TurbulentFlowField& flowField, int i, int j, int k) {
  const RealType visc = flowField.getViscosity().getScalar(i, j, k);
  maxValue_           = std::max(maxValue_, visc);
}
void Stencils::MaxViscStencil::applyRightWall(TurbulentFlowField& flowField, int i, int j, int k) {
  const RealType visc = flowField.getViscosity().getScalar(i, j, k);
  maxValue_           = std::max(maxValue_, visc);
}
void Stencils::MaxViscStencil::applyBottomWall(TurbulentFlowField& flowField, int i, int j, int k) {
  const RealType visc = flowField.getViscosity().getScalar(i, j, k);
  maxValue_           = std::max(maxValue_, visc);
}
void Stencils::MaxViscStencil::applyTopWall(TurbulentFlowField& flowField, int i, int j, int k) {
  const RealType visc = flowField.getViscosity().getScalar(i, j, k);
  maxValue_           = std::max(maxValue_, visc);
}
void Stencils::MaxViscStencil::applyFrontWall(TurbulentFlowField& flowField, int i, int j, int k) {
  const RealType visc = flowField.getViscosity().getScalar(i, j, k);
  maxValue_           = std::max(maxValue_, visc);
}
void Stencils::MaxViscStencil::applyBackWall(TurbulentFlowField& flowField, int i, int j, int k) {
  const RealType visc = flowField.getViscosity().getScalar(i, j, k);
  maxValue_           = std::max(maxValue_, visc);
}

void Stencils::MaxViscStencil::reset() { maxValue_ = 0; }

RealType Stencils::MaxViscStencil::getMaxValue() const { return maxValue_; }