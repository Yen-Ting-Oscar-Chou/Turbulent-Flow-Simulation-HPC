#pragma once

#include "Iterators.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"
#include "BoundaryStencil.hpp"
#include "Definitions.hpp"
#include <vector>

namespace Stencils {
  class VelocityBufferReadStencil: public BoundaryStencil<FlowField> {
  public:
      VelocityBufferReadStencil(const Parameters& parameters);
      ~VelocityBufferReadStencil() override = default;

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
