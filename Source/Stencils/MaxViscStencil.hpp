#pragma once

#include "FieldStencil.hpp"
#include "Parameters.hpp"
#include "TurbulentFlowField.hpp"
#include "BoundaryStencil.hpp"

namespace Stencils {

  /** Compute all velocities on obstacle cells
   */
  class MaxViscStencil: public FieldStencil<TurbulentFlowField>, public BoundaryStencil<TurbulentFlowField> {

  public:
    RealType maxValue_; //! Stores the maximum module of every component
    
    MaxViscStencil();
    ~MaxViscStencil() override = default;

#pragma omp declare target
    void apply(const Parameters& parameters, TurbulentFlowField& flowField, int i, int j) override;
    void apply(const Parameters& parameters, TurbulentFlowField& flowField, int i, int j, int k) override;

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

    RealType getMaxValue() const;
    void           reset();
#pragma omp end declare target
  };

} // namespace Stencils
