#pragma once

#include "FieldStencil.hpp"
#include "Parameters.hpp"
#include "TurbulentFlowField.hpp"

namespace Stencils {

  /** Compute all velocities on obstacle cells
   */
  class MaxViscStencil: public FieldStencil<TurbulentFlowField> {
  private:
    RealType maxValue_; //! Stores the maximum module of every component

  public:
    MaxViscStencil(const Parameters& parameters);
    ~MaxViscStencil() override = default;

    void apply(TurbulentFlowField& flowField, int i, int j) override;
    void apply(TurbulentFlowField& flowField, int i, int j, int k) override;

    RealType getMaxValue() const;
    void           reset();
  };

} // namespace Stencils
