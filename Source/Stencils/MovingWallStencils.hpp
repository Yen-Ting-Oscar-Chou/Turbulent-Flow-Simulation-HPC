#pragma once

#include "BoundaryStencil.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"

namespace Stencils {

  /** Boundary stencil to set constant velocities at the faces.
   *
   * The values for the velocities are stored in the parameters.
   */
  class MovingWallVelocityStencil: public BoundaryStencil<FlowField> {
  public:
    ~MovingWallVelocityStencil() override = default;

#pragma omp declare target
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
#pragma omp end declare target
  };

  /** Boundary stencil to set correct values of FGH at the faces.
   *
   * The values for the velocities are stored in the parameters.
   */
  class MovingWallFGHStencil: public BoundaryStencil<FlowField> {
  public:
    ~MovingWallFGHStencil() override = default;

#pragma omp declare target
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
#pragma omp end declare target
  };

} // namespace Stencils
