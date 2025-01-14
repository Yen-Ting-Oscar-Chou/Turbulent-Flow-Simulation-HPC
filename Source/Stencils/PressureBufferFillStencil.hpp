#pragma once

#include <vector>

#include "BufferBase.hpp"
#include "FlowField.hpp"
#include "Iterators.hpp"
#include "Parameters.hpp"

namespace Stencils {
  class PressureBufferFillStencil: public BufferStencilBase<FlowField> {
  public:
    explicit PressureBufferFillStencil(const Parameters& parameters);

    ~PressureBufferFillStencil() override = default;

    /* Methods for 2D case */
    void applyLeftWall(const Parameters& parameters, FlowField& flowField, int i, int j) override;
    void applyRightWall(const Parameters& parameters, FlowField& flowField, int i, int j) override;
    void applyBottomWall(const Parameters& parameters, FlowField& flowField, int i, int j) override;
    void applyTopWall(const Parameters& parameters, FlowField& flowField, int i, int j) override;

    /* Methods for 3D case */
    void applyLeftWall(const Parameters& parameters, FlowField& flowField, int i, int j, int k) override;
    void applyRightWall(const Parameters& parameters, FlowField& flowField, int i, int j, int k) override;
    void applyBottomWall(const Parameters& parameters, FlowField& flowField, int i, int j, int k) override;
    void applyTopWall(const Parameters& parameters, FlowField& flowField, int i, int j, int k) override;
    void applyFrontWall(const Parameters& parameters, FlowField& flowField, int i, int j, int k) override;
    void applyBackWall(const Parameters& parameters, FlowField& flowField, int i, int j, int k) override;
  };
} // namespace Stencils