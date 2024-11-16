#pragma once

#include <list>

#include "FieldStencil.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"
#include "StencilFunctions.hpp"

namespace Stencils {

  /** Compute list with coordinates of obstacle cells
   */
  class ObstacleCoordinatesStencil: public FieldStencil<FlowField> {
  private:
    std::list<std::tuple<RealType, RealType>>        coordinatesList2D;
    std::list<std::tuple<RealType, RealType, RealType>> coordinatesList3D;

  public:
    ObstacleCoordinatesStencil(const Parameters& parameters, std::list<std::tuple<RealType, RealType>>& coordinatesList);
    ObstacleCoordinatesStencil(const Parameters& parameters, std::list<std::tuple<RealType, RealType, RealType>>& coordinatesList);
    ~ObstacleCoordinatesStencil() override = default;

    void apply(FlowField& flowField, int i, int j) override;
    void apply(FlowField& flowField, int i, int j, int k) override;
  };

} // namespace Stencils
