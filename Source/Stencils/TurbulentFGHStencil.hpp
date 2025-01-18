#pragma once

#include "FieldStencil.hpp"
#include "TurbulentFlowField.hpp"
#include "Parameters.hpp"

namespace Stencils {

  class TurbulentFGHStencil: public FieldStencil<TurbulentFlowField> {

  public:
    TurbulentFGHStencil() = default;
    ~TurbulentFGHStencil() override = default;

#pragma omp declare target
    void apply(const Parameters& parameters, TurbulentFlowField& flowField, int i, int j) override;
    void apply(const Parameters& parameters, TurbulentFlowField& flowField, int i, int j, int k) override;
#pragma omp end declare target
  };

} // namespace Stencils
