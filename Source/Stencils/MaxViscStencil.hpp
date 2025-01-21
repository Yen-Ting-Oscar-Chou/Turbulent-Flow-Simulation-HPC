#pragma once

#include "BoundaryStencil.hpp"
#include "FieldStencil.hpp"
#include "Parameters.hpp"
#include "TurbulentFlowField.hpp"

namespace Stencils {

  /** Compute all velocities on obstacle cells
   */
  class MaxViscStencil: public FieldStencil<TurbulentFlowField>, public BoundaryStencil<TurbulentFlowField> {

  public:
    RealType maxValue_; //! Stores the maximum module of every component

    MaxViscStencil() = default;
    ~MaxViscStencil() override = default;

#pragma omp declare target
    void apply(const Parameters& parameters, TurbulentFlowField& flowField, int i, int j) override { cellMaxValue(parameters, flowField, i, j); }
    void apply(const Parameters& parameters, TurbulentFlowField& flowField, int i, int j, int k) override { cellMaxValue(parameters, flowField, i, j, k); }

    void applyLeftWall(const Parameters& parameters, TurbulentFlowField& flowField, int i, int j) override { cellMaxValue(parameters, flowField, i, j); }
    void applyRightWall(const Parameters& parameters, TurbulentFlowField& flowField, int i, int j) override { cellMaxValue(parameters, flowField, i, j); }
    void applyBottomWall(const Parameters& parameters, TurbulentFlowField& flowField, int i, int j) override { cellMaxValue(parameters, flowField, i, j); }
    void applyTopWall(const Parameters& parameters, TurbulentFlowField& flowField, int i, int j) override { cellMaxValue(parameters, flowField, i, j); }

    void applyLeftWall(const Parameters& parameters, TurbulentFlowField& flowField, int i, int j, int k) override { cellMaxValue(parameters, flowField, i, j, k); }
    void applyRightWall(const Parameters& parameters, TurbulentFlowField& flowField, int i, int j, int k) override { cellMaxValue(parameters, flowField, i, j, k); }
    void applyBottomWall(const Parameters& parameters, TurbulentFlowField& flowField, int i, int j, int k) override { cellMaxValue(parameters, flowField, i, j, k); }
    void applyTopWall(const Parameters& parameters, TurbulentFlowField& flowField, int i, int j, int k) override { cellMaxValue(parameters, flowField, i, j, k); }
    void applyFrontWall(const Parameters& parameters, TurbulentFlowField& flowField, int i, int j, int k) override { cellMaxValue(parameters, flowField, i, j, k); }
    void applyBackWall(const Parameters& parameters, TurbulentFlowField& flowField, int i, int j, int k) override { cellMaxValue(parameters, flowField, i, j, k); }
    void cellMaxValue([[maybe_unused]] const Parameters& parameters, TurbulentFlowField& flowField, int i, int j, int k = 0) {
      const RealType visc = fabs(flowField.getViscosity().getScalar(i, j, k));
      if (visc > maxValue_) {
#pragma omp critical(maxVisc)
        if (visc > maxValue_) {
          maxValue_ = visc;
        }
      }
    }
    RealType getMaxValue() const { return maxValue_; }
    void     reset() { maxValue_ = 0.0; }
#pragma omp end declare target
  };

} // namespace Stencils
