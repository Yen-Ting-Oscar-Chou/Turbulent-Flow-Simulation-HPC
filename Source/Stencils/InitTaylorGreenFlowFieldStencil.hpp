#pragma once

#include "Definitions.hpp"
#include "FieldStencil.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"

namespace Stencils {

  /** Stencil for initialising the taylor-green vortex flow in a periodic domain.
   */
  class InitTaylorGreenFlowFieldStencil: public FieldStencil<FlowField> {
  private:
    const RealType pi2_;

  public:
    InitTaylorGreenFlowFieldStencil();

    void apply(const Parameters& parameters, FlowField& flowField, int i, int j) override;
    void apply(const Parameters& parameters, FlowField& flowField, int i, int j, int k) override;
  };

} // namespace Stencils
