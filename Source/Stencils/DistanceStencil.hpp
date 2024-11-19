#pragma once

#include <list>

#include "FieldStencil.hpp"
#include "TurbulentFlowField.hpp"
#include "Parameters.hpp"
#include "StencilFunctions.hpp"

namespace Stencils {

  /** Computes distance to nearest obstacle for each cell
   */
  class DistanceStencil: public FieldStencil<FlowField> {
  private:
    std::list<std::tuple<RealType, RealType>>           coordinatesList2D;
    std::list<std::tuple<RealType, RealType, RealType>> coordinatesList3D;

  public:
    DistanceStencil(const Parameters& parameters, std::list<std::tuple<RealType, RealType>>& coordinatesList);
    DistanceStencil(const Parameters& parameters, std::list<std::tuple<RealType, RealType, RealType>>& coordinatesList);
    ~DistanceStencil() override = default;
    using FieldStencil::apply;
    void apply(TurbulentFlowField& turbulentField, int i, int j);
    void apply(TurbulentFlowField& turbulentField, int i, int j, int k);
  };

} // namespace Stencils