#include "StdAfx.hpp"

#include "MaxViscStencil.hpp"

Stencils::MaxViscStencil::MaxViscStencil() { reset(); }

void Stencils::MaxViscStencil::apply([[maybe_unused]] [[maybe_unused]] const Parameters& parameters, TurbulentFlowField& flowField, int i, int j) {
  const RealType visc = flowField.getViscosity().getScalar(i, j);
  maxValue_           = fmax(maxValue_, visc);
}

void Stencils::MaxViscStencil::apply([[maybe_unused]] const Parameters& parameters, TurbulentFlowField& flowField, int i, int j, int k) {
  const RealType visc = flowField.getViscosity().getScalar(i, j, k);
  maxValue_           = fmax(maxValue_, visc);
}

void Stencils::MaxViscStencil::applyLeftWall([[maybe_unused]] const Parameters& parameters, TurbulentFlowField& flowField, int i, int j) {
  const RealType visc = flowField.getViscosity().getScalar(i, j);
  maxValue_           = fmax(maxValue_, visc);
}
void Stencils::MaxViscStencil::applyRightWall([[maybe_unused]] const Parameters& parameters, TurbulentFlowField& flowField, int i, int j) {
  const RealType visc = flowField.getViscosity().getScalar(i, j);
  maxValue_           = fmax(maxValue_, visc);
}
void Stencils::MaxViscStencil::applyBottomWall([[maybe_unused]] const Parameters& parameters, TurbulentFlowField& flowField, int i, int j) {
  const RealType visc = flowField.getViscosity().getScalar(i, j);
  maxValue_           = fmax(maxValue_, visc);
}
void Stencils::MaxViscStencil::applyTopWall([[maybe_unused]] const Parameters& parameters, TurbulentFlowField& flowField, int i, int j) {
  const RealType visc = flowField.getViscosity().getScalar(i, j);
  maxValue_           = fmax(maxValue_, visc);
}

void Stencils::MaxViscStencil::applyLeftWall([[maybe_unused]] const Parameters& parameters, TurbulentFlowField& flowField, int i, int j, int k) {
  const RealType visc = flowField.getViscosity().getScalar(i, j, k);
  maxValue_           = fmax(maxValue_, visc);
}
void Stencils::MaxViscStencil::applyRightWall([[maybe_unused]] const Parameters& parameters, TurbulentFlowField& flowField, int i, int j, int k) {
  const RealType visc = flowField.getViscosity().getScalar(i, j, k);
  maxValue_           = fmax(maxValue_, visc);
}
void Stencils::MaxViscStencil::applyBottomWall([[maybe_unused]] const Parameters& parameters, TurbulentFlowField& flowField, int i, int j, int k) {
  const RealType visc = flowField.getViscosity().getScalar(i, j, k);
  maxValue_           = fmax(maxValue_, visc);
}
void Stencils::MaxViscStencil::applyTopWall([[maybe_unused]] const Parameters& parameters, TurbulentFlowField& flowField, int i, int j, int k) {
  const RealType visc = flowField.getViscosity().getScalar(i, j, k);
  maxValue_           = fmax(maxValue_, visc);
}
void Stencils::MaxViscStencil::applyFrontWall([[maybe_unused]] const Parameters& parameters, TurbulentFlowField& flowField, int i, int j, int k) {
  const RealType visc = flowField.getViscosity().getScalar(i, j, k);
  maxValue_           = fmax(maxValue_, visc);
}
void Stencils::MaxViscStencil::applyBackWall([[maybe_unused]] const Parameters& parameters, TurbulentFlowField& flowField, int i, int j, int k) {
  const RealType visc = flowField.getViscosity().getScalar(i, j, k);
  maxValue_           = fmax(maxValue_, visc);
}

void Stencils::MaxViscStencil::reset() { maxValue_ = 0; }

RealType Stencils::MaxViscStencil::getMaxValue() const { return maxValue_; }