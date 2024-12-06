#pragma once

#include "Definitions.hpp"
#include "Iterators.hpp"
#include "Parameters.hpp"
#include "PetscParallelManager.hpp"
#include "Stencils/ViscosityBufferFillStencil.hpp"
#include "Stencils/ViscosityBufferReadStencil.hpp"
#include "TurbulentFlowField.hpp"

namespace ParallelManagers {

  class TurbulentPetscParallelManager: public PetscParallelManager {
  private:
    TurbulentFlowField& _turbulentFlowField; // Reference specific to TurbulentFlowField

    Stencils::ViscosityBufferFillStencil _viscosityBufferFillStencil;
    Stencils::ViscosityBufferReadStencil _viscosityBufferReadStencil;

    ParallelBoundaryIterator<TurbulentFlowField> _parallelBoundaryViscosityFillIterator;
    ParallelBoundaryIterator<TurbulentFlowField> _parallelBoundaryViscosityReadIterator;

  public:
    TurbulentPetscParallelManager(Parameters& parameters, TurbulentFlowField& flowField);
    ~TurbulentPetscParallelManager() = default;

    void communicateViscosity();

    void communicateObstacleCoordinates(std::vector<std::tuple<RealType, RealType>>& coordinatesList2D, std::vector<std::tuple<RealType, RealType, RealType>>& coordinatesList3D);
  };
} // namespace ParallelManagers