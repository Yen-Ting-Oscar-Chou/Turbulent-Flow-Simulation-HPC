#pragma once

#include "BoundaryStencil.hpp"
#include "Definitions.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"
#include <vector>
namespace Stencils {

  class PressureBufferReadStencil: public BoundaryStencil<FlowField> {
  public:
    PressureBufferReadStencil(const Parameters& parameters);
    ~PressureBufferReadStencil() override = default;

    void applyLeftWall(FlowField& flowField, int i, int j) override;
    void applyRightWall(FlowField& flowField, int i, int j) override;
    void applyBottomWall(FlowField& flowField, int i, int j) override;
    void applyTopWall(FlowField& flowField, int i, int j) override;

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
    Parameters parameters_;
  };
} // namespace Stencils