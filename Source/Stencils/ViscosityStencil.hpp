#pragma once

#include "FieldStencil.hpp"
#include "Parameters.hpp"
#include "TurbulentFlowField.hpp"

namespace Stencils {

  class ViscosityStencil: public FieldStencil<TurbulentFlowField> {
  public:
    RealType localVelocity_[27 * 3];
    RealType localMeshsize_[27 * 3];
    RealType coords[3];
    ViscosityStencil()           = default;
    ~ViscosityStencil() override = default;
    // using FieldStencil::apply;

#pragma omp declare target
    void apply(const Parameters& parameters, TurbulentFlowField& turbulentField, int i, int j) override;
    void apply(const Parameters& parameters, TurbulentFlowField& turbulentField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k) override;
    void applyGPU(const Parameters& parameters, VectorField& velocity, ScalarField& viscosity, ScalarField& distance, IntScalarField& flags, int i, int j);
    void applyGPU(const Parameters& parameters, VectorField& velocity, ScalarField& viscosity, ScalarField& distance, IntScalarField& flags, int i, int j, int k);
#pragma omp end declare target
  };

} // namespace Stencils
