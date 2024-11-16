#include "StdAfx.hpp"

#include "TurbulentVTKStencil.hpp"

Stencils::TurbulentVTKStencil::TurbulentVTKStencil(const Parameters& parameters):
  VTKStencil(parameters) {};

void Stencils::TurbulentVTKStencil::apply(TurbulentFlowField& turbulentField, int i, int j) {
  ASSERTION(FieldStencil<FlowField>::parameters_.geometry.dim == 2);

  RealType pressure    = 0.0;
  RealType viscosity   = 0.0;
  RealType distance    = 0.0;
  RealType velocity[2] = {0.0, 0.0};

  if ((turbulentField.getFlags().getValue(i, j) & OBSTACLE_SELF) == 0) {
    turbulentField.getPressureAndVelocity(pressure, velocity, i, j);
    turbulentField.getDistanceAndViscosity(distance, viscosity, i, j);

    pressureStream_ << pressure << std::endl;
    velocityStream_ << velocity[0] << " " << velocity[1] << " 0" << std::endl;
    viscosityStream_ << viscosity << std::endl;
    distanceStream_ << distance << std::endl;
  } else {
    pressureStream_ << "0.0" << std::endl;
    velocityStream_ << "0.0 0.0 0.0" << std::endl;
    viscosityStream_ << "0.0" << std::endl;
    distanceStream_ << "0.0" << std::endl;
  }
}

void Stencils::TurbulentVTKStencil::apply(TurbulentFlowField& turbulentField, int i, int j, int k) {
  ASSERTION(FieldStencil<FlowField>::parameters_.geometry.dim == 3);

  RealType pressure    = 0.0;
  RealType viscosity = 0.0;
  RealType distance = 0.0;
  RealType velocity[3] = {0.0, 0.0, 0.0};

  if ((turbulentField.getFlags().getValue(i, j, k) & OBSTACLE_SELF) == 0) {
    turbulentField.getPressureAndVelocity(pressure, velocity, i, j, k);
    turbulentField.getDistanceAndViscosity(distance, viscosity, i, j, k);

    pressureStream_ << pressure << std::endl;
    velocityStream_ << velocity[0] << " " << velocity[1] << " " << velocity[2] << std::endl;
    viscosityStream_ << viscosity << std::endl;
    distanceStream_ << distance << std::endl;
  } else {
    pressureStream_ << "0.0" << std::endl;
    velocityStream_ << "0.0 0.0 0.0" << std::endl;
    viscosityStream_ << "0.0" << std::endl;
    distanceStream_ << "0.0" << std::endl;
  }
}
