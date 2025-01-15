#pragma once

#include "BoundaryStencil.hpp"
#include "DataStructures.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"

namespace Stencils {

  class NeumannVelocityBoundaryStencil: public BoundaryStencil<FlowField> {
  public:
    ~NeumannVelocityBoundaryStencil() override = default;

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

  class NeumannFGHBoundaryStencil: public BoundaryStencil<FlowField> {
  public:
    ~NeumannFGHBoundaryStencil() override = default;

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
