#include "StdAfx.hpp"

#include "MaxViscStencil.hpp"

Stencils::MaxViscStencil::MaxViscStencil() { reset(); }

void Stencils::MaxViscStencil::apply([[maybe_unused]] [[maybe_unused]] const Parameters& parameters, TurbulentFlowField& flowField, int i, int j) {
  cellMaxValue(parameters, flowField, i, j);
}

void Stencils::MaxViscStencil::apply([[maybe_unused]] const Parameters& parameters, TurbulentFlowField& flowField, int i, int j, int k) {
  cellMaxValue(parameters, flowField, i, j);
}

void Stencils::MaxViscStencil::applyLeftWall([[maybe_unused]] const Parameters& parameters, TurbulentFlowField& flowField, int i, int j) {
  cellMaxValue(parameters, flowField, i, j);
}
void Stencils::MaxViscStencil::applyRightWall([[maybe_unused]] const Parameters& parameters, TurbulentFlowField& flowField, int i, int j) {
  cellMaxValue(parameters, flowField, i, j);
}
void Stencils::MaxViscStencil::applyBottomWall([[maybe_unused]] const Parameters& parameters, TurbulentFlowField& flowField, int i, int j) {
  cellMaxValue(parameters, flowField, i, j);
}
void Stencils::MaxViscStencil::applyTopWall([[maybe_unused]] const Parameters& parameters, TurbulentFlowField& flowField, int i, int j) {
  cellMaxValue(parameters, flowField, i, j);
}

void Stencils::MaxViscStencil::applyLeftWall([[maybe_unused]] const Parameters& parameters, TurbulentFlowField& flowField, int i, int j, int k) {
  cellMaxValue(parameters, flowField, i, j, k);
}
void Stencils::MaxViscStencil::applyRightWall([[maybe_unused]] const Parameters& parameters, TurbulentFlowField& flowField, int i, int j, int k) {
  cellMaxValue(parameters, flowField, i, j, k);
}
void Stencils::MaxViscStencil::applyBottomWall([[maybe_unused]] const Parameters& parameters, TurbulentFlowField& flowField, int i, int j, int k) {
  cellMaxValue(parameters, flowField, i, j, k);
}
void Stencils::MaxViscStencil::applyTopWall([[maybe_unused]] const Parameters& parameters, TurbulentFlowField& flowField, int i, int j, int k) {
  cellMaxValue(parameters, flowField, i, j, k);
}
void Stencils::MaxViscStencil::applyFrontWall([[maybe_unused]] const Parameters& parameters, TurbulentFlowField& flowField, int i, int j, int k) {
  cellMaxValue(parameters, flowField, i, j, k);
}
void Stencils::MaxViscStencil::applyBackWall([[maybe_unused]] const Parameters& parameters, TurbulentFlowField& flowField, int i, int j, int k) {
  cellMaxValue(parameters, flowField, i, j, k);
}

void Stencils::MaxViscStencil::cellMaxValue([[maybe_unused]] const Parameters& parameters, TurbulentFlowField& flowField, int i, int j, int k) {
  const RealType visc = fabs(flowField.getViscosity().getScalar(i, j, k));
  if (visc > maxValue_) {
#pragma omp critical(maxVisc)
    if (visc > maxValue_) {
      maxValue_ = visc;
    }
  }
}

void Stencils::MaxViscStencil::reset() { maxValue_ = 0; }

RealType Stencils::MaxViscStencil::getMaxValue() const { return maxValue_; }