#pragma once

#include <vector>

#include "BufferBase.hpp"
#include "FlowField.hpp"
#include "Iterators.hpp"
#include "Parameters.hpp"
#include "TurbulentFlowField.hpp"

namespace Stencils {
  class ViscosityBufferFillStencil : public BufferStencilBase<TurbulentFlowField> {
  public:
    explicit ViscosityBufferFillStencil(const Parameters& parameters);

    void applyLeftWall(TurbulentFlowField& flowField, int i, int j) override;
    void applyRightWall(TurbulentFlowField& flowField, int i, int j) override;
    void applyBottomWall(TurbulentFlowField& flowField, int i, int j) override;
    void applyTopWall(TurbulentFlowField& flowField, int i, int j) override;

    void applyLeftWall(TurbulentFlowField& flowField, int i, int j, int k) override;
    void applyRightWall(TurbulentFlowField& flowField, int i, int j, int k) override;
    void applyBottomWall(TurbulentFlowField& flowField, int i, int j, int k) override;
    void applyTopWall(TurbulentFlowField& flowField, int i, int j, int k) override;
    void applyFrontWall(TurbulentFlowField& flowField, int i, int j, int k) override;
    void applyBackWall(TurbulentFlowField& flowField, int i, int j, int k) override;
  };
} // namespace Stencils