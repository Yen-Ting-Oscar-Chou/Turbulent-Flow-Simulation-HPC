#pragma once

#include "FlowField.hpp"
#include "Parameters.hpp"

#include "../Stencils/PressureBufferFillStencil.hpp"
#include "../Stencils/PressureBufferReadStencil.hpp"

#include "../Stencils/VelocityBufferFillStencil.hpp"
#include "../Stencils/VelocityBufferReadStencil.hpp"

class PetscParallelManager {

  public:
    FlowField flowfield_;
    const Parameters& _parameters;

    struct BufferSize {
      int LR; // Left-Right
      int TB; // Top-Bottom
      int FB; // Front-Back
    };

    BufferSize presBufSize; // Pressure buffer sizes
    BufferSize velBufSize;  // Velocity buffer sizes

    // TODO: Can we use vectors?
    struct Buffers {
      RealType* leftSend;
      RealType* rightSend;
      RealType* bottomSend;
      RealType* topSend;
      RealType* frontSend;
      RealType* backSend;

      RealType* leftRecv;
      RealType* rightRecv;
      RealType* bottomRecv;
      RealType* topRecv;
      RealType* frontRecv;
      RealType* backRecv;
    };

    Buffers pressureBuffers;
    Buffers velocityBuffers;

    Stencils::PressureBufferFillStencil* _pressureBufferFillStencil;
    Stencils::PressureBufferReadStencil* _pressureBufferReadStencil;

    ParallelBoundaryIterator<FlowField>* _parallelBoundaryPressureFillIterator;
    ParallelBoundaryIterator<FlowField>* _parallelBoundaryPressureReadIterator;
    ParallelBoundaryIterator<FlowField>* _parallelBoundaryVelocityFillIterator;
    ParallelBoundaryIterator<FlowField>* _parallelBoundaryVelocityReadIterator;

    PetscParallelManager(const Parameters & parameters, FlowField & flowField);
    ~PetscParallelManager();

    void communicatePressure();
    void communicateVelocities();
};