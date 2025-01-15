#include "StdAfx.hpp"

#include "RHSStencil.hpp"

void Stencils::RHSStencil::apply(const Parameters& parameters, FlowField& flowField, int i, int j) {
  flowField.getRHS().getScalar(i, j) = 1.0 / parameters.timestep.dt *
        ((flowField.getFGH().getVectorElement(i, j, 0) - flowField.getFGH().getVectorElement(i - 1, j, 0)) / parameters.meshsize->getDx(i, j) +
         (flowField.getFGH().getVectorElement(i, j, 1) - flowField.getFGH().getVectorElement(i, j - 1, 1)) / parameters.meshsize->getDy(i, j));
}

void Stencils::RHSStencil::apply(const Parameters& parameters, FlowField& flowField, int i, int j, int k) {
  flowField.getRHS().getScalar(i, j, k) = 1.0 / parameters.timestep.dt *
        ((flowField.getFGH().getVectorElement(i, j, k, 0) - flowField.getFGH().getVectorElement(i - 1, j, k, 0)) / parameters.meshsize->getDx(i, j, k) +
         (flowField.getFGH().getVectorElement(i, j, k, 1) - flowField.getFGH().getVectorElement(i, j - 1, k, 1)) / parameters.meshsize->getDy(i, j, k) +
         (flowField.getFGH().getVectorElement(i, j, k, 2) - flowField.getFGH().getVectorElement(i, j, k - 1, 2)) / parameters.meshsize->getDz(i, j, k));
}
