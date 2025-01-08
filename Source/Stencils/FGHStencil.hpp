#pragma once

#include "FieldStencil.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"

namespace Stencils {

  class FGHStencil: public FieldStencil<FlowField> {

  public:
    // A local velocity variable that will be used to approximate derivatives. Size matches 3D
    // case, but can be used for 2D as well.
    RealType localVelocity_[27 * 3];
    RealType localMeshsize_[27 * 3];
    FGHStencil() = default;
    ~FGHStencil() override = default;

#pragma omp declare target
    void apply(const Parameters& parameters, FlowField& flowField, int i, int j) override;
    void apply(const Parameters& parameters, FlowField& flowField, int i, int j, int k) override;
    void applyGPU(const Parameters& parameters, VectorField& velocity, VectorField& fgh, [[maybe_unused]] IntScalarField& flags, int i, int j);
    void applyGPU(const Parameters& parameters, VectorField& velocity, VectorField& fgh, [[maybe_unused]] IntScalarField& flags, int i, int j, int k);
#pragma omp end declare target
  };

} // namespace Stencils
