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

    void applyLeftWall(const Parameters& parameters, TurbulentFlowField& flowField, int i, int j) override;
    void applyRightWall(const Parameters& parameters, TurbulentFlowField& flowField, int i, int j) override;
    void applyBottomWall(const Parameters& parameters, TurbulentFlowField& flowField, int i, int j) override;
    void applyTopWall(const Parameters& parameters, TurbulentFlowField& flowField, int i, int j) override;

    void applyLeftWall(const Parameters& parameters, TurbulentFlowField& flowField, int i, int j, int k) override;
    void applyRightWall(const Parameters& parameters, TurbulentFlowField& flowField, int i, int j, int k) override;
    void applyBottomWall(const Parameters& parameters, TurbulentFlowField& flowField, int i, int j, int k) override;
    void applyTopWall(const Parameters& parameters, TurbulentFlowField& flowField, int i, int j, int k) override;
    void applyFrontWall(const Parameters& parameters, TurbulentFlowField& flowField, int i, int j, int k) override;
    void applyBackWall(const Parameters& parameters, TurbulentFlowField& flowField, int i, int j, int k) override;
  };
} // namespace Stencils