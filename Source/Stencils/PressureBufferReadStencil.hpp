#pragma once

#include "BoundaryStencil.hpp"
#include "Definitions.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"

namespace Stencils {

  class PressureBufferReadStencil: public BoundaryStencil<FlowField> {
  public:
    PressureBufferReadStencil(const Parameters& parameters, const std::vector<RealType>& pressure);
    ~PressureBufferReadStencil() override = default;

    void applyLeftWall(FlowField& flowField, int i, int j) override;
    void applyRightWall(FlowField& flowField, int i, int j) override;
    void applyBottomWall(FlowField& flowField, int i, int j) override;
    void applyTopWall(FlowField& flowField, int i, int j) override;

    void applyLeftWall(FlowField& flowField, int i, int j, int k) override;
    void applyRightWall(FlowField& flowField, int i, int j, int k) override;
    void applyBottomWall(FlowField& flowField, int i, int j, int k) override;
    void applyTopWall(FlowField& flowField, int i, int j, int k) override;
    void applyFrontWall(FlowField& flowField, int i, int j, int k) override;
    void applyBackWall(FlowField& flowField, int i, int j, int k) override;

  private:
    const std::vector<RealType>& pressure_;
  };
} // namespace Stencils