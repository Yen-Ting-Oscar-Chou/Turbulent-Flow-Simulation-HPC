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
      int _left;
      int _right; 
      int _top;
      int _bottom;
      int _front;
      int _back;
      
      Stencils::PressureBufferFillStencil *_pressureBufferFillStencil;
      Stencils::PressureBufferReadStencil *_pressureBufferReadStencil;

      Stencils::VelocityBufferFillStencil *_velocityBufferFillStencil;
      Stencils::VelocityBufferReadStencil *_velocityBufferReadStencil;

      ParallelBoundaryIterator<FlowField> *_parallelBoundaryPressureFillIterator;
      ParallelBoundaryIterator<FlowField> *_parallelBoundaryPressureReadIterator;
      ParallelBoundaryIterator<FlowField> *_parallelBoundaryVelocityFillIterator;
      ParallelBoundaryIterator<FlowField> *_parallelBoundaryVelocityReadIterator;
    
    public:
      PetscParallelManager(Parameters& parameters, FlowField& flowField);
      ~PetscParallelManager();

      void communicatePressure();
      void communicateVelocities();
  };
} // namespace ParallelManagers