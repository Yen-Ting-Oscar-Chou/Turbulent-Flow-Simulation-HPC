#pragma once

#include "FieldStencil.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"

namespace Stencils {

  /** Stencil to compute the velocity once the pressure has been found.
   */
  class VelocityStencil: public FieldStencil<FlowField> {
  public:
    VelocityStencil() = default;
    ~VelocityStencil() override = default;

#pragma omp declare target
    void apply(const Parameters& parameters, FlowField& flowField, int i, int j) override;
    void apply(const Parameters& parameters, FlowField& flowField, int i, int j, int k) override;
    void applyGPU(const Parameters& parameters, VectorField& velocity, VectorField& fgh, ScalarField& pressure, IntScalarField& flags, int i, int j);
    void applyGPU(const Parameters& parameters, VectorField& velocity, VectorField& fgh, ScalarField& pressure, IntScalarField& flags, int i, int j, int k);
#pragma omp end declare target
  };

} // namespace Stencils
