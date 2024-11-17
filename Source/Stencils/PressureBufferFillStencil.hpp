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

    std::vector<RealType>& selectPressureLeft(const Parameters& parameters) {
      if (parameters.geometry.dim == 2) {
        static std::vector<RealType> pressureLeft2D;
        return pressureLeft2D;
      } else {
        static std::vector<RealType> pressureLeft3D;
        return pressureLeft3D;
      }
    }

    std::vector<RealType>& selectPressureRight(const Parameters& parameters) {
      if (parameters.geometry.dim == 2) {
        static std::vector<RealType> pressureRight2D;
        return pressureRight2D;
      } else {
        static std::vector<RealType> pressureRight3D;
        return pressureRight3D;
      }
    }

    std::vector<RealType>& selectPressureBottom(const Parameters& parameters) {
      static std::vector<RealType> pressureBottom;
      return pressureBottom;
    }

    std::vector<RealType>& selectPressureTop(const Parameters& parameters) {
      static std::vector<RealType> pressureTop;
      return pressureTop;
    }

    std::vector<RealType>& selectPressureFront(const Parameters& parameters) {
      static std::vector<RealType> pressureFront;
      return pressureFront;
    }

    std::vector<RealType>& selectPressureBack(const Parameters& parameters) {
      static std::vector<RealType> pressureBack;
      return pressureBack;
    }

  protected:
    static RealType &getScalar(FlowField &flowField, int i, int j);
    static RealType &getScalar(FlowField &flowField, int i, int j, int k);

  private:
    std::vector<RealType>& pressureLeft_;
    std::vector<RealType>& pressureRight_;
    std::vector<RealType>& pressureBottom_;
    std::vector<RealType>& pressureTop_;
    std::vector<RealType>& pressureFront_;
    std::vector<RealType>& pressureBack_;

  };
} // namespace Stencils