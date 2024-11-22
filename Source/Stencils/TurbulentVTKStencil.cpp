#include "StdAfx.hpp"

#include "TurbulentVTKStencil.hpp"

Stencils::TurbulentVTKStencil::TurbulentVTKStencil(const Parameters& parameters):
  VTKStencil(parameters) {};

void Stencils::TurbulentVTKStencil::apply(TurbulentFlowField& turbulentField, int i, int j) {
  ASSERTION(FieldStencil<TurbulentFlowField>::parameters_.geometry.dim == 2);

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
  ASSERTION(FieldStencil<TurbulentFlowField>::parameters_.geometry.dim == 3);

  RealType pressure    = 0.0;
  RealType viscosity   = 0.0;
  RealType distance    = 0.0;
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

void Stencils::TurbulentVTKStencil::write(TurbulentFlowField& turbulentField, int timeStep, RealType simulationTime) {
  openFile(timeStep, simulationTime);

  if (FieldStencil<TurbulentFlowField>::parameters_.geometry.dim == 2) {
    // Write pressure
    ofile_ << "CELL_DATA " << turbulentField.getNx() * turbulentField.getNy() << std::endl << "SCALARS pressure float 1" << std::endl << "LOOKUP_TABLE default" << std::endl;
    ofile_ << pressureStream_.str() << std::endl;
    pressureStream_.str("");

    ofile_ << "SCALARS distance float 1" << std::endl << "LOOKUP_TABLE default" << std::endl;
    ofile_ << distanceStream_.str() << std::endl;
    distanceStream_.str("");

    ofile_ << "SCALARS turbulent_viscosity float 1" << std::endl << "LOOKUP_TABLE default" << std::endl;
    ofile_ << viscosityStream_.str() << std::endl;
    viscosityStream_.str("");

    // Write velocity
    ofile_ << "VECTORS velocity float" << std::endl;
    ofile_ << velocityStream_.str() << std::endl;
    velocityStream_.str("");
  }

  if (FieldStencil<TurbulentFlowField>::parameters_.geometry.dim == 3) {
    // Write pressure
    ofile_
      << "CELL_DATA " << turbulentField.getNx() * turbulentField.getNy() * turbulentField.getNz() << std::endl
      << "SCALARS pressure float 1" << std::endl
      << "LOOKUP_TABLE default" << std::endl;
    ofile_ << pressureStream_.str() << std::endl;
    pressureStream_.str("");

    ofile_ << "SCALARS distance float 1" << std::endl << "LOOKUP_TABLE default" << std::endl;
    ofile_ << distanceStream_.str() << std::endl;
    distanceStream_.str("");

    ofile_ << "SCALARS turbulent_viscosity float 1" << std::endl << "LOOKUP_TABLE default" << std::endl;
    ofile_ << viscosityStream_.str() << std::endl;
    viscosityStream_.str("");

    // Write velocity
    ofile_ << "VECTORS velocity float" << std::endl;
    ofile_ << velocityStream_.str() << std::endl;
    velocityStream_.str("");
  }

  written_ = true;
  closeFile();
}
