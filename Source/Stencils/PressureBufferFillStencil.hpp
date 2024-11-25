#pragma once

#include "Iterators.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"
#include <vector>

namespace Stencils {
  class PressureBufferFillStencil: public BoundaryStencil<FlowField> {
  public:
    PressureBufferFillStencil(const Parameters& parameters);

    ~PressureBufferFillStencil() override = default;

    /* Methods for 2D case */
    void applyLeftWall(FlowField& flowField, int i, int j) override;
    void applyRightWall(FlowField& flowField, int i, int j) override;
    void applyBottomWall(FlowField& flowField, int i, int j) override;
    void applyTopWall(FlowField& flowField, int i, int j) override;

    /* Methods for 3D case */
    void applyLeftWall(FlowField& flowField, int i, int j, int k) override;
    void applyRightWall(FlowField& flowField, int i, int j, int k) override;
    void applyBottomWall(FlowField& flowField, int i, int j, int k) override;
    void applyTopWall(FlowField& flowField, int i, int j, int k) override;
    void applyFrontWall(FlowField& flowField, int i, int j, int k) override;
    void applyBackWall(FlowField& flowField, int i, int j, int k) override;

    std::vector<RealType>& getpressureLeft();
    std::vector<RealType>& getpressureRight();
    std::vector<RealType>& getpressureBottom();
    std::vector<RealType>& getpressureTop();
    std::vector<RealType>& getpressureFront();
    std::vector<RealType>& getpressureBack();

  private:
    std::vector<RealType> pressureLeft_;
    std::vector<RealType> pressureRight_;
    std::vector<RealType> pressureBottom_;
    std::vector<RealType> pressureTop_;
    std::vector<RealType> pressureFront_;
    std::vector<RealType> pressureBack_;

  };
} // namespace Stencils