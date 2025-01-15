#pragma once

#include "BoundaryStencil.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"

namespace Stencils {

  /** A stencil to set the input velocity in channel flows. Only implements the applyLeftWall(...) methods.
   */
  class BFInputVelocityStencil: public BoundaryStencil<FlowField> {
  private:
    void initStepSize(const Parameters& parameters);

  public:
    RealType stepSize_; //! Fixes the size of the step. If zero, is just channel flow.
    BFInputVelocityStencil();
    ~BFInputVelocityStencil() override = default;

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

  /** FGH stencil which corresponds to one-sided Dirichlet conditions (i.e. the channel flow profile).
   *  Only implements the applyLeftWall(...) methods.
   */
  class BFInputFGHStencil: public BoundaryStencil<FlowField> {
  private:
    void initStepSize(const Parameters& parameters);

  public:
    RealType stepSize_; //! Fixes the size of the step. If zero, is just channel flow.
    BFInputFGHStencil();
    ~BFInputFGHStencil() override = default;

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
