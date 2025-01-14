#pragma once

#include <vector>

#include "BoundaryStencil.hpp"
#include "BufferBase.hpp"
#include "Definitions.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"
namespace Stencils {

  class PressureBufferReadStencil: public BufferStencilBase<FlowField> {
  public:
    explicit PressureBufferReadStencil(const Parameters& parameters);
    ~PressureBufferReadStencil() override = default;

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