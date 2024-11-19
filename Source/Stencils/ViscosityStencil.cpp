#include "StdAfx.hpp"

#include "ViscosityStencil.hpp"
#include "StencilFunctions.hpp"


Stencils::ViscosityStencil::ViscosityStencil(const Parameters& parameters):
  FieldStencil<TurbulentFlowField>(parameters) {}

void Stencils::ViscosityStencil::apply(TurbulentFlowField& turbulentField, int i, int j) {
  const int obstacle = turbulentField.getFlags().getValue(i, j);
  if ((obstacle & OBSTACLE_SELF) == 1) {
    return;
  }
  loadLocalVelocity2D(turbulentField, localVelocity_, i, j);
  loadLocalMeshsize2D(parameters_, localMeshsize_, i, j);
  
  RealType distance = turbulentField.getDistance().getScalar(i, j);

  RealType coords[2] = {0.0, 0.0};
  computeGlobalCoordinates(coords, parameters_, i, j);
  const RealType x = coords[0];

  const RealType re           = parameters_.flow.Re;
  const RealType mixingLength = std::min(KAPPA * distance, 0.09 * parameters_.turbulence.deltaMixLen(x, re));

  // 1/2 * (dudy + dvdx)
  // 0.5 * ((dudy + dudz) + (dvdx + dwdx))
  const RealType S_ij = 0.5 * (dudy(localVelocity_, localMeshsize_) + dvdx(localVelocity_, localMeshsize_));
  RealType& v_t = turbulentField.getViscosity().getScalar(i, j);
  v_t = mixingLength * sqrt(2 * S_ij * S_ij);
}

void Stencils::ViscosityStencil::apply(TurbulentFlowField& turbulentField, int i, int j, int k) {}