#pragma once

#include "FieldStencil.hpp"
#include "TurbulentFlowField.hpp"
#include "Parameters.hpp"
#include "StencilFunctions.hpp"

namespace Stencils {

  /** Computes distance to nearest obstacle for each cell
   */
  class DistanceStencil: public FieldStencil<TurbulentFlowField> {
  private:
    std::vector<std::tuple<RealType, RealType>>           coordinatesList2D;
    std::vector<std::tuple<RealType, RealType, RealType>> coordinatesList3D;

  public:
    DistanceStencil(const Parameters& parameters, std::vector<std::tuple<RealType, RealType>>& coordinatesList);
    DistanceStencil(const Parameters& parameters, std::vector<std::tuple<RealType, RealType, RealType>>& coordinatesList);
    ~DistanceStencil() override = default;
    void apply(TurbulentFlowField& turbulentField, int i, int j) override;
    void apply(TurbulentFlowField& turbulentField, int i, int j, int k) override;
  };

} // namespace Stencils