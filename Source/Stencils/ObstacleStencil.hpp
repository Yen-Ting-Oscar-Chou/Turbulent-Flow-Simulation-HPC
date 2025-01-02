#pragma once

#include "FieldStencil.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"

namespace Stencils {

  /** Compute all velocities on obstacle cells
   */
  class ObstacleStencil: public FieldStencil<FlowField> {
  public:
    ObstacleStencil() = default;
    ~ObstacleStencil() override = default;

#pragma omp declare target
    void apply(const Parameters& parameters, FlowField& flowField, int i, int j) override;
    void apply(const Parameters& parameters, FlowField& flowField, int i, int j, int k) override;
#pragma omp end declare target
  };

} // namespace Stencils
