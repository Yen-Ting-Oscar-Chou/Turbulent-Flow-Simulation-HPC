#include "StdAfx.hpp"

#include "ViscosityStencil.hpp"

#include "StencilFunctions.hpp"

void Stencils::ViscosityStencil::apply(const Parameters& parameters, TurbulentFlowField& turbulentField, int i, int j) {
  const int obstacle = turbulentField.getFlags().getValue(i, j);
  if ((obstacle & OBSTACLE_SELF) == 1) {
    return;
  }
  loadLocalVelocity2D(turbulentField, localVelocity_, i, j);
  loadLocalMeshsize2D(parameters, localMeshsize_, i, j);

  const RealType distance = turbulentField.getDistance().getScalar(i, j);

  RealType coords[2] = {0.0, 0.0};
  computeGlobalCoordinates(coords, parameters, i, j);
  const RealType x = std::fabs(coords[0]);

  const RealType       re             = parameters.flow.Re;
  const TurbulenceType turbulenceType = parameters.turbulence.deltaMixLen;
  const RealType       deltaMixLen    = (turbulenceType == ZERO) * deltaZero(x, re) + (turbulenceType == LAMINAR) * deltaLaminar(x, re)
                               + (turbulenceType == TURBULENT) * deltaTurbulent(x, re);
  const RealType mixingLength = std::min(KAPPA * distance, 0.09 * deltaMixLen);

  // S_ij S_ij = S11^2 + S_22^2 + 2S_12^2
  // (dudx)^2 + (dvdy)^2 + (dudy + dvdx)^2
  // v_t = l_m^2 * sqrt(2*S_ij*S_ij)
  const RealType S_11    = dudx(localVelocity_, localMeshsize_);
  const RealType S_22    = dvdy(localVelocity_, localMeshsize_);
  const RealType S_12    = 0.5 * (dudy(localVelocity_, localMeshsize_) + dvdx(localVelocity_, localMeshsize_));
  const RealType S_ij_sq = S_11 * S_11 + S_22 * S_22 + 2 * S_12 * S_12;
  RealType&      v_t     = turbulentField.getViscosity().getScalar(i, j);
  v_t                    = mixingLength * mixingLength * sqrt(2 * S_ij_sq);

  // For BFS
  if ((obstacle & OBSTACLE_BOTTOM) == OBSTACLE_BOTTOM) {
    turbulentField.getViscosity().getScalar(i, j - 1) = -v_t;
  }

  if ((obstacle & OBSTACLE_LEFT) == OBSTACLE_LEFT) {
    turbulentField.getViscosity().getScalar(i - 1, j) = -v_t;
  }
}

void Stencils::ViscosityStencil::apply(const Parameters& parameters, TurbulentFlowField& turbulentField, int i, int j, int k) {
  const int obstacle = turbulentField.getFlags().getValue(i, j, k);
  if ((obstacle & OBSTACLE_SELF) == 1) {
    return;
  }
  loadLocalVelocity3D(turbulentField, localVelocity_, i, j, k);
  loadLocalMeshsize3D(parameters, localMeshsize_, i, j, k);

  const RealType distance = turbulentField.getDistance().getScalar(i, j, k);

  RealType coords[3] = {0.0, 0.0, 0.0};
  computeGlobalCoordinates(coords, parameters, i, j, k);
  const RealType x = std::fabs(coords[0]);

  const RealType re           = parameters.flow.Re;
  const TurbulenceType turbulenceType = parameters.turbulence.deltaMixLen;
  const RealType       deltaMixLen    = (turbulenceType == ZERO) * deltaZero(x, re) + (turbulenceType == LAMINAR) * deltaLaminar(x, re)
                               + (turbulenceType == TURBULENT) * deltaTurbulent(x, re);
  const RealType mixingLength = std::min(KAPPA * distance, 0.09 * deltaMixLen);

  const RealType S_11    = dudx(localVelocity_, localMeshsize_);
  const RealType S_22    = dvdy(localVelocity_, localMeshsize_);
  const RealType S_33    = dwdz(localVelocity_, localMeshsize_);
  const RealType S_12    = 0.5 * (dudy(localVelocity_, localMeshsize_) + dvdx(localVelocity_, localMeshsize_));
  const RealType S_13    = 0.5 * (dudz(localVelocity_, localMeshsize_) + dwdx(localVelocity_, localMeshsize_));
  const RealType S_23    = 0.5 * (dvdz(localVelocity_, localMeshsize_) + dwdy(localVelocity_, localMeshsize_));
  const RealType S_ij_sq = S_11 * S_11 + S_22 * S_22 + S_33 * S_33 + 2 * (S_12 * S_12 + S_13 * S_13 + S_23 * S_23);
  RealType&      v_t     = turbulentField.getViscosity().getScalar(i, j, k);
  v_t                    = mixingLength * mixingLength * sqrt(2 * S_ij_sq);

  // For BFS
  if ((obstacle & OBSTACLE_BOTTOM) == OBSTACLE_BOTTOM) {
    turbulentField.getViscosity().getScalar(i, j - 1, k) = -v_t;
  }

  if ((obstacle & OBSTACLE_LEFT) == OBSTACLE_LEFT) {
    turbulentField.getViscosity().getScalar(i - 1, j, k) = -v_t;
  }
}