#pragma once

#include <vector>

#include "BoundaryStencil.hpp"
#include "BufferBase.hpp"
#include "Definitions.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"
#include "TurbulentFlowField.hpp"
namespace Stencils {

  class ViscosityBufferReadStencil: public BufferStencilBase<TurbulentFlowField> {
  public:
    explicit ViscosityBufferReadStencil(const Parameters& parameters);
    ~ViscosityBufferReadStencil() override = default;

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