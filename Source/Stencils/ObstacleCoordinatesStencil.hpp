#pragma once

#include "BoundaryStencil.hpp"
#include "FieldStencil.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"
#include "StencilFunctions.hpp"

namespace Stencils {

  /** Compute list with coordinates of obstacle cells
   */
  class ObstacleCoordinatesStencil: public FieldStencil<FlowField>, public BoundaryStencil<FlowField> {
  private:
    std::vector<std::tuple<RealType, RealType>>&           coordinatesList2D;
    std::vector<std::tuple<RealType, RealType, RealType>>& coordinatesList3D;

    void addToList(int i, int j);
    void addToList(int i, int j, int k);

  public:
    ObstacleCoordinatesStencil(
      const Parameters& parameters, std::vector<std::tuple<RealType, RealType>>& coordinatesList2D, std::vector<std::tuple<RealType, RealType, RealType>>& coordinatesList3D
    );
    ~ObstacleCoordinatesStencil() override = default;

    void apply(FlowField& flowField, int i, int j) override;
    void apply(FlowField& flowField, int i, int j, int k) override;
    void applyLeftWall([[maybe_unused]] FlowField& flowField, int i, int j) override;
    void applyRightWall([[maybe_unused]] FlowField& flowField, int i, int j) override;
    void applyBottomWall([[maybe_unused]] FlowField& flowField, int i, int j) override;
    void applyTopWall([[maybe_unused]] FlowField& flowField, int i, int j) override;
    void applyLeftWall([[maybe_unused]] FlowField& flowField, int i, int j, int k) override;
    void applyRightWall([[maybe_unused]] FlowField& flowField, int i, int j, int k) override;
    void applyBottomWall([[maybe_unused]] FlowField& flowField, int i, int j, int k) override;
    void applyTopWall([[maybe_unused]] FlowField& flowField, int i, int j, int k) override;
    void applyFrontWall([[maybe_unused]] FlowField& flowField, int i, int j, int k) override;
    void applyBackWall([[maybe_unused]] FlowField& flowField, int i, int j, int k) override;
  };

} // namespace Stencils
