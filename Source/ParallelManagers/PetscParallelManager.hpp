#pragma once

#include "Definitions.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"

#include "Stencils/PressureBufferFillStencil.hpp"
#include "Stencils/PressureBufferReadStencil.hpp"

#include "Stencils/VelocityBufferFillStencil.hpp"
#include "Stencils/VelocityBufferReadStencil.hpp"


#include "Iterators.hpp"

namespace ParallelManagers {

  class PetscParallelManager {
    private: 
      Parameters& _parameters;
      FlowField& _flowField;
    public:
      /*
      PetscParallelManager(const PetscParallelManager&) = delete;
      PetscParallelManager& operator=(const PetscParallelManager&) = delete;

      // Allow move operations, if needed
      PetscParallelManager(PetscParallelManager&&) = default;
      PetscParallelManager& operator=(PetscParallelManager&&) = default;
      */
      struct Buffers {
        std::vector<RealType> leftSend;
        std::vector<RealType> rightSend;
        std::vector<RealType> bottomSend;
        std::vector<RealType> topSend;
        std::vector<RealType> frontSend;
        std::vector<RealType> backSend;

        std::vector<RealType> leftRecv;
        std::vector<RealType> rightRecv;
        std::vector<RealType> bottomRecv;
        std::vector<RealType> topRecv;
        std::vector<RealType> frontRecv;
        std::vector<RealType> backRecv;
      };

      Buffers pressureBuffers;
      Buffers velocityBuffers;

      Stencils::PressureBufferFillStencil *_pressureBufferFillStencil;
      Stencils::PressureBufferReadStencil *_pressureBufferReadStencil;

      Stencils::VelocityBufferFillStencil *_velocityBufferFillStencil;
      Stencils::VelocityBufferReadStencil *_velocityBufferReadStencil;

      ParallelBoundaryIterator<FlowField> *_parallelBoundaryPressureFillIterator;
      ParallelBoundaryIterator<FlowField> *_parallelBoundaryPressureReadIterator;
      ParallelBoundaryIterator<FlowField> *_parallelBoundaryVelocityFillIterator;
      ParallelBoundaryIterator<FlowField> *_parallelBoundaryVelocityReadIterator;

      PetscParallelManager(Parameters& parameters, FlowField& flowField);
      ~PetscParallelManager();

      void communicatePressure();
      void communicateVelocities();
      int computeRankFromIndices(int i, int j, int k) const;
  };
} // namespace ParallelManagers