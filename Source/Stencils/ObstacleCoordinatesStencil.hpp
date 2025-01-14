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

    void addToList(const Parameters& parameters, int i, int j);
    void addToList(const Parameters& parameters, int i, int j, int k);

  public:
  // TODO refactor stencil!
    ObstacleCoordinatesStencil(std::vector<std::tuple<RealType, RealType>>& coordinatesList2D, std::vector<std::tuple<RealType, RealType, RealType>>& coordinatesList3D);
    ~ObstacleCoordinatesStencil() override = default;

    void apply(const Parameters& parameters, FlowField& flowField, int i, int j) override;
    void apply(const Parameters& parameters, FlowField& flowField, int i, int j, int k) override;
    void applyLeftWall(const Parameters& parameters, FlowField& flowField, int i, int j) override;
    void applyRightWall(const Parameters& parameters, FlowField& flowField, int i, int j) override;
    void applyBottomWall(const Parameters& parameters, FlowField& flowField, int i, int j) override;
    void applyTopWall(const Parameters& parameters, FlowField& flowField, int i, int j) override;
    void applyLeftWall(const Parameters& parameters, FlowField& flowField, int i, int j, int k) override;
    void applyRightWall(const Parameters& parameters, FlowField& flowField, int i, int j, int k) override;
    void applyBottomWall(const Parameters& parameters, FlowField& flowField, int i, int j, int k) override;
    void applyTopWall(const Parameters& parameters, FlowField& flowField, int i, int j, int k) override;
    void applyFrontWall(const Parameters& parameters, FlowField& flowField, int i, int j, int k) override;
    void applyBackWall(const Parameters& parameters, FlowField& flowField, int i, int j, int k) override;
  };

} // namespace Stencils
