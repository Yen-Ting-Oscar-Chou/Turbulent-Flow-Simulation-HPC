#pragma once

#include "BoundaryStencil.hpp"
#include "Definitions.hpp"
#include "TurbulentFlowField.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"
#include <vector>
namespace Stencils {

  class ViscosityBufferReadStencil: public BoundaryStencil<TurbulentFlowField> {
  public:
    ViscosityBufferReadStencil(const Parameters& parameters);
    ~ViscosityBufferReadStencil() override = default;

    void applyLeftWall(TurbulentFlowField& flowField, int i, int j) override;
    void applyRightWall(TurbulentFlowField& flowField, int i, int j) override;
    void applyBottomWall(TurbulentFlowField& flowField, int i, int j) override;
    void applyTopWall(TurbulentFlowField& flowField, int i, int j) override;

    void applyLeftWall(TurbulentFlowField& flowField, int i, int j, int k) override;
    void applyRightWall(TurbulentFlowField& flowField, int i, int j, int k) override;
    void applyBottomWall(TurbulentFlowField& flowField, int i, int j, int k) override;
    void applyTopWall(TurbulentFlowField& flowField, int i, int j, int k) override;
    void applyFrontWall(TurbulentFlowField& flowField, int i, int j, int k) override;
    void applyBackWall(TurbulentFlowField& flowField, int i, int j, int k) override;

    std::vector<RealType>& getviscosityLeft();
    std::vector<RealType>& getviscosityRight();
    std::vector<RealType>& getviscosityBottom();
    std::vector<RealType>& getviscosityTop();
    std::vector<RealType>& getviscosityFront();
    std::vector<RealType>& getviscosityBack();
  private:
    std::vector<RealType> viscosityLeft_;
    std::vector<RealType> viscosityRight_;
    std::vector<RealType> viscosityBottom_;
    std::vector<RealType> viscosityTop_;
    std::vector<RealType> viscosityFront_;
    std::vector<RealType> viscosityBack_;
  };
} // namespace Stencils