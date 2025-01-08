#include "StdAfx.hpp"

#include "RHSStencil.hpp"

void Stencils::RHSStencil::apply(const Parameters& parameters, FlowField& flowField, int i, int j) {
  flowField.getRHS().getScalar(i, j) = 1.0 / parameters.timestep.dt *
        ((flowField.getFGH().getVector(i, j)[0] - flowField.getFGH().getVector(i - 1, j)[0]) / parameters.meshsize.getDx(i, j) +
         (flowField.getFGH().getVector(i, j)[1] - flowField.getFGH().getVector(i, j - 1)[1]) / parameters.meshsize.getDy(i, j));
}

void Stencils::RHSStencil::apply(const Parameters& parameters, FlowField& flowField, int i, int j, int k) {
  flowField.getRHS().getScalar(i, j, k) = 1.0 / parameters.timestep.dt *
        ((flowField.getFGH().getVector(i, j, k)[0] - flowField.getFGH().getVector(i - 1, j, k)[0]) / parameters.meshsize.getDx(i, j, k) +
         (flowField.getFGH().getVector(i, j, k)[1] - flowField.getFGH().getVector(i, j - 1, k)[1]) / parameters.meshsize.getDy(i, j, k) +
         (flowField.getFGH().getVector(i, j, k)[2] - flowField.getFGH().getVector(i, j, k - 1)[2]) / parameters.meshsize.getDz(i, j, k));
}

void Stencils::RHSStencil::applyGPU(const Parameters& parameters, ScalarField& rhs, VectorField& fgh, int i, int j) {
  rhs.getScalar(i, j) = 1.0 / parameters.timestep.dt *
          ((fgh.getVectorElement(i, j, 0) - fgh.getVectorElement(i - 1, j, 0)) / parameters.meshsize.getDx(i, j) +
          (fgh.getVectorElement(i, j, 1) - fgh.getVectorElement(i, j - 1, 1)) / parameters.meshsize.getDy(i, j));
}

void Stencils::RHSStencil::applyGPU(const Parameters& parameters, ScalarField& rhs, VectorField& fgh, int i, int j, int k) {
  rhs.getScalar(i, j, k) = 1.0 / parameters.timestep.dt *
          ((fgh.getVectorElement(i, j, k, 0) - fgh.getVectorElement(i - 1, j, k, 0)) / parameters.meshsize.getDx(i, j, k) +
          (fgh.getVectorElement(i, j, k, 1) - fgh.getVectorElement(i, j - 1, k, 1)) / parameters.meshsize.getDy(i, j, k) + 
          (fgh.getVectorElement(i, j, k, 2) - fgh.getVectorElement(i, j, k - 1, 2)) / parameters.meshsize.getDz(i, j, k));
}
