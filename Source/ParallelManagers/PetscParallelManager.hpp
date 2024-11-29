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