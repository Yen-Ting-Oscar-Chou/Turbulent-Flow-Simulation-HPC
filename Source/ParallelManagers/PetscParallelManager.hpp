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
      FlowField& _flowField;

      Stencils::PressureBufferFillStencil _pressureBufferFillStencil;
      Stencils::PressureBufferReadStencil _pressureBufferReadStencil;

      Stencils::VelocityBufferFillStencil _velocityBufferFillStencil;
      Stencils::VelocityBufferReadStencil _velocityBufferReadStencil;

      ParallelBoundaryIterator<FlowField> _parallelBoundaryPressureFillIterator;
      ParallelBoundaryIterator<FlowField> _parallelBoundaryPressureReadIterator;
      ParallelBoundaryIterator<FlowField> _parallelBoundaryVelocityFillIterator;
      ParallelBoundaryIterator<FlowField> _parallelBoundaryVelocityReadIterator;

    protected:
      int _left;

      int _front;

      int _back;

      int _top;

      int _bottom;

      int _right;

      Parameters& _parameters;

    public:
      PetscParallelManager(Parameters& parameters, FlowField& flowField);
      ~PetscParallelManager() = default;

      void communicatePressure();
      void communicateVelocities();
  };
} // namespace ParallelManagers