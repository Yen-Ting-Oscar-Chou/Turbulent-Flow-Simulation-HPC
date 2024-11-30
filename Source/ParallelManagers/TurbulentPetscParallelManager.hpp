#pragma once

#include "Definitions.hpp"
#include "TurbulentFlowField.hpp"
#include "Parameters.hpp"

#include "Stencils/ViscosityBufferFillStencil.hpp"
#include "Stencils/ViscosityBufferReadStencil.hpp"


#include "Iterators.hpp"

namespace ParallelManagers {

  class TurbulentPetscParallelManager {
    private: 
      Parameters& _parameters;
      TurbulentFlowField& _flowField;
      int _left;
      int _right; 
      int _top;
      int _bottom;
      int _front;
      int _back;
      
      Stencils::ViscosityBufferFillStencil *_viscosityBufferFillStencil;
      Stencils::ViscosityBufferReadStencil *_viscosityBufferReadStencil;

      ParallelBoundaryIterator<TurbulentFlowField> *_parallelBoundaryViscosityFillIterator;
      ParallelBoundaryIterator<TurbulentFlowField> *_parallelBoundaryViscosityReadIterator;
    
    public:
      TurbulentPetscParallelManager(Parameters& parameters, TurbulentFlowField& flowField);
      ~TurbulentPetscParallelManager();

      void communicateViscosity();
  };
} // namespace ParallelManagers