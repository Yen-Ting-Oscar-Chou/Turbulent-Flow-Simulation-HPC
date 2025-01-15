#pragma once

#include "FieldStencil.hpp"
#include "TurbulentFlowField.hpp"
#include "Parameters.hpp"

namespace Stencils {

  class TurbulentFGHStencil: public FieldStencil<TurbulentFlowField> {

  public:
    // A local velocity variable that will be used to approximate derivatives. Size matches 3D
    // case, but can be used for 2D as well.
    RealType localVelocity_[27 * 3];
    RealType localMeshsize_[27 * 3];
    RealType localViscosity_[27 * 3];
    TurbulentFGHStencil() = default;
    ~TurbulentFGHStencil() override = default;

#pragma omp declare target
    void apply(const Parameters& parameters, TurbulentFlowField& flowField, int i, int j) override;
    void apply(const Parameters& parameters, TurbulentFlowField& flowField, int i, int j, int k) override;
#pragma omp end declare target
  };

} // namespace Stencils
