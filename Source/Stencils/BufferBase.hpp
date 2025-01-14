#pragma once

#include <vector>

#include "../Definitions.hpp"
#include "BoundaryStencil.hpp"
#include "../Parameters.hpp"

namespace Stencils {
  class BufferBase {
  public:
    virtual ~BufferBase() = default;

    std::vector<RealType>& getBufferLeft() { return bufferLeft_; }
    std::vector<RealType>& getBufferRight() { return bufferRight_; }
    std::vector<RealType>& getBufferBottom() { return bufferBottom_; }
    std::vector<RealType>& getBufferTop() { return bufferTop_; }
    std::vector<RealType>& getBufferFront() { return bufferFront_; }
    std::vector<RealType>& getBufferBack() { return bufferBack_; }

  protected:
    std::vector<RealType> bufferLeft_;
    std::vector<RealType> bufferRight_;
    std::vector<RealType> bufferBottom_;
    std::vector<RealType> bufferTop_;
    std::vector<RealType> bufferFront_;
    std::vector<RealType> bufferBack_;
  };

  template <typename FlowFieldType>
  class BufferStencilBase : public BoundaryStencil<FlowFieldType>, public BufferBase {
  public:
    explicit BufferStencilBase() {}

    ~BufferStencilBase() override = default;

    void applyLeftWall(const Parameters& parameters, FlowFieldType& flowField, int i, int j) override = 0;
    void applyRightWall(const Parameters& parameters, FlowFieldType& flowField, int i, int j) override = 0;
    void applyBottomWall(const Parameters& parameters, FlowFieldType& flowField, int i, int j) override = 0;
    void applyTopWall(const Parameters& parameters, FlowFieldType& flowField, int i, int j) override = 0;

    void applyLeftWall(const Parameters& parameters, FlowFieldType& flowField, int i, int j, int k) override = 0;
    void applyRightWall(const Parameters& parameters, FlowFieldType& flowField, int i, int j, int k) override = 0;
    void applyBottomWall(const Parameters& parameters, FlowFieldType& flowField, int i, int j, int k) override = 0;
    void applyTopWall(const Parameters& parameters, FlowFieldType& flowField, int i, int j, int k) override = 0;
    void applyFrontWall(const Parameters& parameters, FlowFieldType& flowField, int i, int j, int k) override = 0;
    void applyBackWall(const Parameters& parameters, FlowFieldType& flowField, int i, int j, int k) override = 0;
  };
}


