#pragma once

#include "FieldStencil.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"

namespace Stencils {

  class FGHStencil: public FieldStencil<FlowField> {

  public:
    FGHStencil() = default;
    ~FGHStencil() override = default;

#pragma omp declare target
    void apply(const Parameters& parameters, FlowField& flowField, int i, int j) override;
    void apply(const Parameters& parameters, FlowField& flowField, int i, int j, int k) override;
#pragma omp end declare target
  };

} // namespace Stencils
