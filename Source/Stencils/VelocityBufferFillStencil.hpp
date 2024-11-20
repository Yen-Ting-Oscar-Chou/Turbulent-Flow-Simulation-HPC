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
    // What is the purpose for selectVelocity?
    static std::vector<RealType>& selectVelocityLeft(const Parameters& parameters) {
      if (parameters.geometry.dim == 2) {
        static std::vector<RealType> velocityLeft2D;
        return velocityLeft2D;
      } else {
        static std::vector<RealType> velocityLeft3D;
        return velocityLeft3D;
      }
    }
    // bug: return pressure for a velocity stencil
    static std::vector<RealType>& selectVelocityRight(const Parameters& parameters) {
      if (parameters.geometry.dim == 2) {
        static std::vector<RealType> pressureRight2D;
        return pressureRight2D;
      } else {
        static std::vector<RealType> pressureRight3D;
        return pressureRight3D;
      }
    }
    // Why are the below sections different from the two above?
    static std::vector<RealType>& selectVelocityBottom(const Parameters& parameters) {
      static std::vector<RealType> velocityBottom;
      return velocityBottom;
    }

    static std::vector<RealType>& selectVelocityTop(const Parameters& parameters) {
      static std::vector<RealType> velocityTop;
      return velocityTop;
    }

    static std::vector<RealType>& selectVelocityFront(const Parameters& parameters) {
      static std::vector<RealType> velocityFront;
      return velocityFront;
    }

    static std::vector<RealType>& selectVelocityBack(const Parameters& parameters) {
      static std::vector<RealType> velocityBack;
      return velocityBack;
    }
  // What is the purpose for getVector?
  protected:
    static RealType &getVector(FlowField &flowField, int i, int j);
    static RealType &getVector(FlowField &flowField, int i, int j, int k);

  private:
    std::vector<RealType>& velocityLeft_;
    std::vector<RealType>& velocityRight_;
    std::vector<RealType>& velocityBottom_;
    std::vector<RealType>& velocityTop_;
    std::vector<RealType>& velocityFront_;
    std::vector<RealType>& velocityBack_;
  };

} // namespace Stencils
