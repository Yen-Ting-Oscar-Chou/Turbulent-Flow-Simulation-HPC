#pragma once

#include "FieldStencil.hpp"
#include "Parameters.hpp"
#include "StencilFunctions.hpp"
#include "TurbulentFlowField.hpp"

namespace Stencils {

  class ViscosityStencil: FieldStencil<TurbulentFlowField> {
  public:
    ViscosityStencil()           = default;
    ~ViscosityStencil() override = default;

#pragma omp declare target
    void apply(const Parameters& parameters, TurbulentFlowField& turbulentField, int i, int j) override {
      const int obstacle = turbulentField.getFlags().getValue(i, j);
      if ((obstacle & OBSTACLE_SELF) == 1) {
        return;
      }

      VectorField& velocity  = turbulentField.getVelocity();
      ScalarField& viscosity = turbulentField.getViscosity();

      const RealType distance = turbulentField.getDistance().getScalar(i, j);

      RealType x_;
      computeGlobalXCoordinate(x_, parameters, i, j);
      const RealType x = std::fabs(x_);

      const RealType       re             = parameters.flow.Re;
      const TurbulenceType turbulenceType = parameters.turbulence.deltaMixLen;
      const RealType       deltaMixLen    = (turbulenceType == ZERO) * deltaZero(x, re) + (turbulenceType == LAMINAR) * deltaLaminar(x, re)
                                   + (turbulenceType == TURBULENT) * deltaTurbulent(x, re);
      const RealType mixingLength = std::min(KAPPA * distance, 0.09 * deltaMixLen);

      // S_ij S_ij = S11^2 + S_22^2 + 2S_12^2
      // (dudx)^2 + (dvdy)^2 + (dudy + dvdx)^2
      // v_t = l_m^2 * sqrt(2*S_ij*S_ij)
      const RealType S_11    = dudx(velocity, parameters, i, j);
      const RealType S_22    = dvdy(velocity, parameters, i, j);
      const RealType S_12    = 0.5 * (dudy(velocity, parameters, i, j) + dvdx(velocity, parameters, i, j));
      const RealType S_ij_sq = S_11 * S_11 + S_22 * S_22 + 2 * S_12 * S_12;
      const RealType v_t     = mixingLength * mixingLength * sqrt(2 * S_ij_sq);

      viscosity.getScalar(i, j) = v_t;

      // For BFS
      if ((obstacle & OBSTACLE_BOTTOM) == OBSTACLE_BOTTOM) {
        viscosity.getScalar(i, j - 1) = -v_t;
      }

      if ((obstacle & OBSTACLE_LEFT) == OBSTACLE_LEFT) {
        viscosity.getScalar(i - 1, j) = -v_t;
      }
    }
    void apply(const Parameters& parameters, TurbulentFlowField& turbulentField, int i, int j, int k) override {
      const int obstacle = turbulentField.getFlags().getValue(i, j, k);
      if ((obstacle & OBSTACLE_SELF) == 1) {
        return;
      }
      VectorField& velocity  = turbulentField.getVelocity();
      ScalarField& viscosity = turbulentField.getViscosity();

      const RealType distance = turbulentField.getDistance().getScalar(i, j, k);

      RealType x_;
      computeGlobalXCoordinate(x_, parameters, i, j, k);
      const RealType x = std::fabs(x_);

      const RealType       re             = parameters.flow.Re;
      const TurbulenceType turbulenceType = parameters.turbulence.deltaMixLen;
      const RealType       deltaMixLen    = (turbulenceType == ZERO) * deltaZero(x, re) + (turbulenceType == LAMINAR) * deltaLaminar(x, re)
                                   + (turbulenceType == TURBULENT) * deltaTurbulent(x, re);
      const RealType mixingLength = std::min(KAPPA * distance, 0.09 * deltaMixLen);

      const RealType S_11    = dudx(velocity, parameters, i, j, k);
      const RealType S_22    = dvdy(velocity, parameters, i, j, k);
      const RealType S_33    = dwdz(velocity, parameters, i, j, k);
      const RealType S_12    = 0.5 * (dudy(velocity, parameters, i, j, k) + dvdx(velocity, parameters, i, j, k));
      const RealType S_13    = 0.5 * (dudz(velocity, parameters, i, j, k) + dwdx(velocity, parameters, i, j, k));
      const RealType S_23    = 0.5 * (dvdz(velocity, parameters, i, j, k) + dwdy(velocity, parameters, i, j, k));
      const RealType S_ij_sq = S_11 * S_11 + S_22 * S_22 + S_33 * S_33 + 2 * (S_12 * S_12 + S_13 * S_13 + S_23 * S_23);
      RealType&      v_t     = viscosity.getScalar(i, j, k);
      v_t                    = mixingLength * mixingLength * sqrt(2 * S_ij_sq);

      // For BFS
      if ((obstacle & OBSTACLE_BOTTOM) == OBSTACLE_BOTTOM) {
        viscosity.getScalar(i, j - 1, k) = -v_t;
      }

      if ((obstacle & OBSTACLE_LEFT) == OBSTACLE_LEFT) {
        viscosity.getScalar(i - 1, j, k) = -v_t;
      }
    }
#pragma omp end declare target
  };

} // namespace Stencils
