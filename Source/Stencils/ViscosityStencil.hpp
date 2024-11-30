#pragma once

#include "FieldStencil.hpp"
#include "Parameters.hpp"
#include "TurbulentFlowField.hpp"

namespace Stencils {

  class ViscosityStencil: public FieldStencil<TurbulentFlowField> {
  private:
    RealType localVelocity_[27 * 3];
    RealType localMeshsize_[27 * 3];
  public:
    ViscosityStencil(const Parameters& parameters);
    ~ViscosityStencil() override = default;
    using FieldStencil::apply;
    void apply(TurbulentFlowField& turbulentField, int i, int j);
    void apply(TurbulentFlowField& turbulentField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k);
  };

} // namespace Stencils
