#pragma once

#include "VTKStencil.hpp"
#include "TurbulentFlowField.hpp"

namespace Stencils {
  class TurbulentVTKStencil : public VTKStencil<TurbulentFlowField>{
  private:
    std::stringstream viscosityStream_; //! Stream for the viscosity data
    std::stringstream distanceStream_;  //! Stream for the distance data

  public:
    TurbulentVTKStencil(const Parameters& parameters);
    ~TurbulentVTKStencil() override = default;
    
    void apply(TurbulentFlowField& turbulentField, int i, int j) override;
    void apply(TurbulentFlowField& turbulentField, int i, int j, int k) override;

    void write(TurbulentFlowField& turbulentField, int timeStep, RealType simulationTime);
  };
} // namespace Stencils