#pragma once

#include "FieldStencil.hpp"
#include "Parameters.hpp"
#include "TurbulentFlowField.hpp"
#include "BoundaryStencil.hpp"

namespace Stencils {

  /** Compute all velocities on obstacle cells
   */
  class MaxViscStencil: public FieldStencil<TurbulentFlowField>, public BoundaryStencil<TurbulentFlowField> {
  private:
    RealType maxValue_; //! Stores the maximum module of every component

  public:
    MaxViscStencil(const Parameters& parameters);
    ~MaxViscStencil() override = default;

    void apply(TurbulentFlowField& flowField, int i, int j) override;
    void apply(TurbulentFlowField& flowField, int i, int j, int k) override;

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

    RealType getMaxValue() const;
    void           reset();
  };

} // namespace Stencils
