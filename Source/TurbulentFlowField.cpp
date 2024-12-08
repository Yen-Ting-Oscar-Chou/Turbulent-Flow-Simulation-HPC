#include "StdAfx.hpp"
#include "TurbulentFlowField.hpp"

TurbulentFlowField::TurbulentFlowField(const Parameters& parameters):
  FlowField(parameters),
  viscosity_(
    parameters.geometry.dim == 2
      ? ScalarField(parameters.parallel.localSize[0] + 3, parameters.parallel.localSize[1] + 3)
      : ScalarField(parameters.parallel.localSize[0] + 3, parameters.parallel.localSize[1] + 3, parameters.parallel.localSize[2] + 3)
  ),
  distance_(
    parameters.geometry.dim == 2
      ? ScalarField(parameters.parallel.localSize[0] + 3, parameters.parallel.localSize[1] + 3)
      : ScalarField(parameters.parallel.localSize[0] + 3, parameters.parallel.localSize[1] + 3, parameters.parallel.localSize[2] + 3)
  ) {}

ScalarField& TurbulentFlowField::getViscosity() { return viscosity_; }
ScalarField& TurbulentFlowField::getDistance() { return distance_; }

void TurbulentFlowField::getDistanceAndViscosity(RealType& distance, RealType& viscosity, int i, int j) {
  distance  = getDistance().getScalar(i, j);
  viscosity = getViscosity().getScalar(i, j);
}
void TurbulentFlowField::getDistanceAndViscosity(RealType& distance, RealType& viscosity, int i, int j, int k) {
  distance  = getDistance().getScalar(i, j, k);
  viscosity = getViscosity().getScalar(i, j, k);
}