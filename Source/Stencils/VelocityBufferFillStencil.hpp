#pragma once

#include "Iterators.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"
#include <vector>

namespace Stencils {

  class VelocityBufferFillStencil : public BoundaryStencil<FlowField> {
  public:
    VelocityBufferFillStencil(const Parameters& parameters);
    ~VelocityBufferFillStencil() override = default;

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

    static std::vector<RealType>& selectVelocityLeft(const Parameters& parameters) {
      if (parameters.geometry.dim == 2) {
        static std::vector<RealType> velocityLeft2D;
        return velocityLeft2D;
      } else {
        static std::vector<RealType> velocityLeft3D;
        return velocityLeft3D;
      }
    }

    std::vector<RealType>& getvelocityLeft();
    std::vector<RealType>& getvelocityRight();
    std::vector<RealType>& getvelocityTop();
    std::vector<RealType>& getvelocityBottom();
    std::vector<RealType>& getvelocityFront();
    std::vector<RealType>& getvelocityBack();

  private:
    std::vector<RealType> velocityLeft_;
    std::vector<RealType> velocityRight_;
    std::vector<RealType> velocityBottom_;
    std::vector<RealType> velocityTop_;
    std::vector<RealType> velocityFront_;
    std::vector<RealType> velocityBack_;
  };

} // namespace Stencils
