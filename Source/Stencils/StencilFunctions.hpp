#pragma once

#include "StdAfx.hpp"

#include "Definitions.hpp"
#include "Parameters.hpp"
#include "TurbulentFlowField.hpp"

namespace Stencils {

  inline RealType computeDistance(RealType* coords1, RealType* coords2, const Parameters& parameters_) {
    RealType result = 0.0;
    for (int i = 0; i < parameters_.geometry.dim; i++) {
      result += pow(coords1[i] - coords2[i], 2);
    }
    return sqrt(result);
  }

  inline void computeGlobalCoordinates(RealType* coords, const Parameters& parameters_, int i, int j) {
    coords[0] = parameters_.meshsize->getPosX(i, j, -1) + 0.5 * parameters_.meshsize->getDx(i, j, -1);
    coords[1] = parameters_.meshsize->getPosY(i, j, -1) + 0.5 * parameters_.meshsize->getDy(i, j, -1);
  }

  inline void computeGlobalCoordinates(RealType* coords, const Parameters& parameters_, int i, int j, int k) {
    coords[0] = parameters_.meshsize->getPosX(i, j, k) + 0.5 * parameters_.meshsize->getDx(i, j, k);
    coords[1] = parameters_.meshsize->getPosY(i, j, k) + 0.5 * parameters_.meshsize->getDy(i, j, k);
    coords[2] = parameters_.meshsize->getPosZ(i, j, k) + 0.5 * parameters_.meshsize->getDz(i, j, k);
  }

  // Load the local velocity cube with relevant velocities of the 2D plane
  inline void loadLocalVelocity2D(FlowField& flowField, RealType* const localVelocity, int i, int j) {
    for (int row = -1; row <= 1; row++) {
      for (int column = -1; column <= 1; column++) {
        const RealType* const point                  = flowField.getVelocity().getVector(i + column, j + row);
        localVelocity[39 + 9 * row + 3 * column]     = point[0]; // x-component
        localVelocity[39 + 9 * row + 3 * column + 1] = point[1]; // y-component
      }
    }
  }

  // Load the local velocity cube with surrounding velocities
  inline void loadLocalVelocity3D(FlowField& flowField, RealType* const localVelocity, int i, int j, int k) {
    for (int layer = -1; layer <= 1; layer++) {
      for (int row = -1; row <= 1; row++) {
        for (int column = -1; column <= 1; column++) {
          const RealType* const point                               = flowField.getVelocity().getVector(i + column, j + row, k + layer);
          localVelocity[39 + 27 * layer + 9 * row + 3 * column]     = point[0]; // x-component
          localVelocity[39 + 27 * layer + 9 * row + 3 * column + 1] = point[1]; // y-component
          localVelocity[39 + 27 * layer + 9 * row + 3 * column + 2] = point[2]; // z-component
        }
      }
    }
  }

  // Load local meshsize for 2D -> same as loadLocalVelocity2D, but invoking call to meshsize-ptr
  inline void loadLocalMeshsize2D(const Parameters& parameters, RealType* const localMeshsize, int i, int j) {
    for (int row = -1; row <= 1; row++) {
      for (int column = -1; column <= 1; column++) {
        localMeshsize[39 + 9 * row + 3 * column]     = parameters.meshsize->getDx(i + column, j + row);
        localMeshsize[39 + 9 * row + 3 * column + 1] = parameters.meshsize->getDy(i + column, j + row);
      }
    }
  }

  // Load local meshsize for 3D
  inline void loadLocalMeshsize3D(const Parameters& parameters, RealType* const localMeshsize, int i, int j, int k) {
    for (int layer = -1; layer <= 1; layer++) {
      for (int row = -1; row <= 1; row++) {
        for (int column = -1; column <= 1; column++) {
          localMeshsize[39 + 27 * layer + 9 * row + 3 * column]     = parameters.meshsize->getDx(i + column, j + row, k + layer);
          localMeshsize[39 + 27 * layer + 9 * row + 3 * column + 1] = parameters.meshsize->getDy(i + column, j + row, k + layer);
          localMeshsize[39 + 27 * layer + 9 * row + 3 * column + 2] = parameters.meshsize->getDz(i + column, j + row, k + layer);
        }
      }
    }
  }

  // Load local viscosity for 2D -> same as loadLocalVelocity2D, but invoking call to meshsize-ptr
  inline void loadLocalViscosity2D(const Parameters& parameters, TurbulentFlowField& turbulentField, RealType* const localViscosity, int i, int j) {
    for (int row = -1; row <= 1; row++) {
      for (int column = -1; column <= 1; column++) {
        localViscosity[39 + 9 * row + 3 * column]     = 1 / parameters.flow.Re + turbulentField.getViscosity().getScalar(i + column, j + row);
        localViscosity[39 + 9 * row + 3 * column + 1] = 1 / parameters.flow.Re + turbulentField.getViscosity().getScalar(i + column, j + row);
      }
    }
  }

  // Load local meshsize for 3D
  inline void loadLocalViscosity3D(const Parameters& parameters, TurbulentFlowField& turbulentField, RealType* const localViscosity, int i, int j, int k) {
    for (int layer = -1; layer <= 1; layer++) {
      for (int row = -1; row <= 1; row++) {
        for (int column = -1; column <= 1; column++) {
          localViscosity[39 + 27 * layer + 9 * row + 3 * column]     = 1 / parameters.flow.Re + turbulentField.getViscosity().getScalar(i + column, j + row, k + layer);
          localViscosity[39 + 27 * layer + 9 * row + 3 * column + 1] = 1 / parameters.flow.Re + turbulentField.getViscosity().getScalar(i + column, j + row, k + layer);
          localViscosity[39 + 27 * layer + 9 * row + 3 * column + 2] = 1 / parameters.flow.Re + turbulentField.getViscosity().getScalar(i + column, j + row, k + layer);
        }
      }
    }
  }

  // Maps an index and a component to the corresponding value in the cube.
  inline int mapd(int i, int j, int k, int component) { return 39 + 27 * k + 9 * j + 3 * i + component; }

  inline RealType firstDerivative(const RealType* const lv, const RealType* const lm, const COMP comp, const DERIV deriv) {
    const bool isX = DERIVX == deriv;
    const bool isY = DERIVY == deriv;
    const bool isZ = DERIVZ == deriv;

    const int indexThis  = mapd(0, 0, 0, comp);
    const int indexMinus = mapd(isX * -1, isY * -1, isZ * -1, comp);
    const int indexPlus  = mapd(isX * 1, isY * 1, isZ * 1, comp);
    const int indexDm    = mapd(isX * -1, isY * -1, isZ * -1, deriv);
    const int indexDp    = mapd(0, 0, 0, deriv);

    const RealType dp = lm[indexDp];
    const RealType dm = lm[indexDm];

    return 1 / (2 * (dp * dm)) * (dm * lv[indexPlus] + (dp - dm) * lv[indexThis] - dp * lv[indexMinus]);
  }

  inline RealType secondDerivative(const RealType* const lv, const RealType* const lm, const COMP comp, const DERIV deriv) {
    const bool isX = DERIVX == deriv;
    const bool isY = DERIVY == deriv;
    const bool isZ = DERIVZ == deriv;

    const int indexThis  = mapd(0, 0, 0, comp);
    const int indexMinus = mapd(isX * -1, isY * -1, isZ * -1, comp);
    const int indexPlus  = mapd(isX * 1, isY * 1, isZ * 1, comp);
    const int indexDm    = mapd(isX * -1, isY * -1, isZ * -1, deriv);
    const int indexDp    = mapd(0, 0, 0, deriv);

    const RealType dp = lm[indexDp];
    const RealType dm = lm[indexDm];

    const RealType fac = 2 / (dm * dp * dp + dp * dm * dm);
    return fac * (dm * (lv[indexPlus] - lv[indexThis]) + dp * (lv[indexMinus] - lv[indexThis]));
  }

  // Derivative functions. They are applied to a cube of 3x3x3 cells. lv stands for the local velocity, lm represents
  // the local mesh sizes dudx <-> first derivative of u-component of velocity field w.r.t. x-direction.
  inline RealType dudx(const RealType* const lv, const RealType* const lm) { return firstDerivative(lv, lm, COMPX, DERIVX); }

  inline RealType dudy(const RealType* const lv, const RealType* const lm) { return firstDerivative(lv, lm, COMPX, DERIVY); }

  inline RealType dudz(const RealType* const lv, const RealType* const lm) { return firstDerivative(lv, lm, COMPX, DERIVZ); }

  inline RealType dvdx(const RealType* const lv, const RealType* const lm) { return firstDerivative(lv, lm, COMPY, DERIVX); }

  inline RealType dvdy(const RealType* const lv, const RealType* const lm) { return firstDerivative(lv, lm, COMPY, DERIVY); }

  inline RealType dvdz(const RealType* const lv, const RealType* const lm) { return firstDerivative(lv, lm, COMPY, DERIVZ); }

  inline RealType dwdx(const RealType* const lv, const RealType* const lm) { return firstDerivative(lv, lm, COMPZ, DERIVX); }

  inline RealType dwdy(const RealType* const lv, const RealType* const lm) { return firstDerivative(lv, lm, COMPZ, DERIVY); }

  inline RealType dwdz(const RealType* const lv, const RealType* const lm) { return firstDerivative(lv, lm, COMPZ, DERIVZ); }

  inline RealType interpolateViscosity(const RealType visc_bl, const RealType visc_br, const RealType visc_tl, const RealType visc_tr) {
    return 0.25 * (visc_bl + visc_br + visc_tl + visc_tr);
  }

  inline RealType viscositySingleDerivative(const RealType* const lVel, const RealType* const lVis, const RealType* const lMesh, const RealType re, const COMP comp) {
    const bool  isX   = COMPX == comp;
    const bool  isY   = COMPY == comp;
    const bool  isZ   = COMPZ == comp;
    const DERIV deriv = static_cast<DERIV>(comp);

    const int indexThis  = mapd(0, 0, 0, comp);
    const int indexMinus = mapd(isX * -1, isY * -1, isZ * -1, comp);
    const int indexPlus  = mapd(isX * 1, isY * 1, isZ * 1, comp);
    const int indexDm    = mapd(isX * -1, isY * -1, isZ * -1, deriv);
    const int indexDp    = mapd(0, 0, 0, deriv);

    const RealType dm       = lMesh[indexDm];
    const RealType dp       = lMesh[indexDp];
    const RealType d        = (dm + dp) / 2;
    const RealType forward  = (lVel[indexPlus] - lVel[indexThis]) / dp;
    const RealType backward = (lVel[indexThis] - lVel[indexMinus]) / dm;

    const RealType vis_plus = lVis[indexPlus] + (1 / re);
    const RealType vis_this = lVis[indexThis] + (1 / re);

    return 1 / d * (vis_plus * forward - vis_this * backward);
  }

  inline RealType viscosityDoubleDerivative(
    const RealType* const lVel, const RealType* const lVis, const RealType* const lMesh, const RealType re, const COMP comp1, const COMP comp2
  ) {
    const bool  isX1   = COMPX == comp1;
    const bool  isY1   = COMPY == comp1;
    const bool  isZ1   = COMPZ == comp1;
    const bool  isX2   = COMPX == comp2;
    const bool  isY2   = COMPY == comp2;
    const bool  isZ2   = COMPZ == comp2;
    const DERIV deriv1 = static_cast<DERIV>(comp1);
    const DERIV deriv2 = static_cast<DERIV>(comp2);

    const int indexThis1       = mapd(0, 0, 0, comp2);
    const int indexThis2       = mapd(0, 0, 0, comp1);
    const int indexMinus1      = mapd(isX1 * -1, isY1 * -1, isZ1 * -1, comp2);
    const int indexPlus2Minus1 = mapd(isX2 + isX1 * -1, isY2 + isY1 * -1, isZ2 + isZ1 * -1, comp1);
    const int indexPlus1       = mapd(isX1 * 1, isY1 * 1, isZ1 * 1, comp2);
    const int indexPlus2       = mapd(isX2 * 1, isY2 * 1, isZ2 * 1, comp1);
    const int indexDm1         = mapd(isX1 * -1, isY1 * -1, isZ1 * -1, deriv1);
    const int indexDm2         = mapd(isX1 * -1, isY1 * -1, isZ1 * -1, deriv2);
    const int indexDp1         = mapd(0, 0, 0, deriv1);
    const int indexDp2         = mapd(0, 0, 0, deriv2);

    const RealType dm1 = lMesh[indexDm1];
    const RealType dm2 = lMesh[indexDm2];
    const RealType dp1 = lMesh[indexDp1];
    const RealType dp2 = lMesh[indexDp2];
    const RealType d   = (dm1 + dp1) / 2; // (j-1+j+1)/2

    const int      indexBottomLeft1  = indexThis1;
    const int      indexBottomRight1 = indexPlus2;
    const int      indexTopLeft1     = mapd(isX1, isY1, isZ1, 0);
    const int      indexTopRight1    = mapd(isX1 + isX2, isY1 + isY2, isZ1 + isZ2, 0);
    const RealType visc1             = interpolateViscosity(lVis[indexBottomLeft1], lVis[indexBottomRight1], lVis[indexTopLeft1], lVis[indexTopRight1]) + (1 / re);

    const int      indexBottomLeft2  = indexMinus1;
    const int      indexBottomRight2 = mapd(isX1 * -1 + isX2, isY1 * -1 + isY2, isZ1 * -1 + isZ2, 0);
    const int      indexTopLeft2     = indexThis2;
    const int      indexTopRight2    = mapd(isX2, isY2, isZ2, 0);
    const RealType visc2             = interpolateViscosity(lVis[indexBottomLeft2], lVis[indexBottomRight2], lVis[indexTopLeft2], lVis[indexTopRight2]) + (1 / re);

    const RealType forward1  = (lVel[indexPlus1] - lVel[indexThis1]) / dp1;
    const RealType forward2  = (lVel[indexPlus2] - lVel[indexThis2]) / dp2;
    const RealType backward1 = (lVel[indexThis1] - lVel[indexMinus1]) / dm1;
    const RealType backward2 = (lVel[indexPlus2Minus1] - lVel[indexMinus1]) / dm2;

    return 1 / d * (visc1 * (forward1 + forward2) - visc2 * (backward1 + backward2));
  }

  // d/dx * (v * du/dx)
  inline RealType dVdxdudx(const RealType* const lVel, const RealType* const lVis, const RealType* const lMesh, const RealType re) {
    return viscositySingleDerivative(lVel, lVis, lMesh, re, COMPX);
  }

  // d/dy * (v * dv/dy)
  inline RealType dVdydvdy(const RealType* const lVel, const RealType* const lVis, const RealType* const lMesh, const RealType re) {
    return viscositySingleDerivative(lVel, lVis, lMesh, re, COMPY);
  }

  // d/dz * (v * dw/dz)
  inline RealType dVdzdwdz(const RealType* const lVel, const RealType* const lVis, const RealType* const lMesh, const RealType re) {
    return viscositySingleDerivative(lVel, lVis, lMesh, re, COMPZ);
  }

  // d/dy * (v * (du/dy + dv/dx))
  inline RealType dVdydudydvdx(const RealType* const lVel, const RealType* const lVis, const RealType* const lMesh, const RealType re) {
    return viscosityDoubleDerivative(lVel, lVis, lMesh, re, COMPY, COMPX);
  }

  // d/dz * (v * (du/dz + dw/dx))
  inline RealType dVdzdudzdwdx(const RealType* const lVel, const RealType* const lVis, const RealType* const lMesh, const RealType re) {
    return viscosityDoubleDerivative(lVel, lVis, lMesh, re, COMPZ, COMPX);
  }

  // d/dx * (v * (dv/dx + du/dy))
  inline RealType dVdxdvdxdudy(const RealType* const lVel, const RealType* const lVis, const RealType* const lMesh, const RealType re) {
    return viscosityDoubleDerivative(lVel, lVis, lMesh, re, COMPX, COMPY);
  }

  // d/dz * (v * (dv/dz + dw/dy))
  inline RealType dVdzdvdzdwdy(const RealType* const lVel, const RealType* const lVis, const RealType* const lMesh, const RealType re) {
    return viscosityDoubleDerivative(lVel, lVis, lMesh, re, COMPZ, COMPY);
  }

  // d/dx * (v * (dw/dx + du/dz))
  inline RealType dVdxdwdxdudz(const RealType* const lVel, const RealType* const lVis, const RealType* const lMesh, const RealType re) {
    return viscosityDoubleDerivative(lVel, lVis, lMesh, re, COMPX, COMPZ);
  }

  // d/dy * (v * (dw/dy + dv/dz))
  inline RealType dVdydwdydvdz(const RealType* const lVel, const RealType* const lVis, const RealType* const lMesh, const RealType re) {
    return viscosityDoubleDerivative(lVel, lVis, lMesh, re, COMPY, COMPZ);
  }

  // TODO WS1: Second derivatives
  // Second derivative of u w.r.t. x-direction
  inline RealType d2udx2(const RealType* const lv, const RealType* const lm) { return secondDerivative(lv, lm, COMPX, DERIVX); }

  // Second derivative of u w.r.t. y-direction
  inline RealType d2udy2(const RealType* const lv, const RealType* const lm) { return secondDerivative(lv, lm, COMPX, DERIVY); }

  // Second derivative of u w.r.t. z-direction
  inline RealType d2udz2(const RealType* const lv, const RealType* const lm) { return secondDerivative(lv, lm, COMPX, DERIVZ); }

  // Second derivative of v w.r.t. x-direction
  inline RealType d2vdx2(const RealType* const lv, const RealType* const lm) { return secondDerivative(lv, lm, COMPY, DERIVX); }

  // Second derivative of v w.r.t. y-direction
  inline RealType d2vdy2(const RealType* const lv, const RealType* const lm) { return secondDerivative(lv, lm, COMPY, DERIVY); }

  // Second derivative of v w.r.t. z-direction
  inline RealType d2vdz2(const RealType* const lv, const RealType* const lm) { return secondDerivative(lv, lm, COMPY, DERIVZ); }

  // Second derivative of w w.r.t. x-direction
  inline RealType d2wdx2(const RealType* const lv, const RealType* const lm) { return secondDerivative(lv, lm, COMPZ, DERIVX); }

  // Second derivative of w w.r.t. y-direction
  inline RealType d2wdy2(const RealType* const lv, const RealType* const lm) { return secondDerivative(lv, lm, COMPZ, DERIVY); }

  // Second derivative of w w.r.t. z-direction
  inline RealType d2wdz2(const RealType* const lv, const RealType* const lm) { return secondDerivative(lv, lm, COMPZ, DERIVZ); }

  // First derivative of product (u*v), evaluated at the location of the v-component.
  inline RealType duvdx(const RealType* const lv, const Parameters& parameters, const RealType* const lm) {
#ifndef NDEBUG
    const RealType tmp1 = 1.0 / 4.0 * ((((lv[mapd(0, 0, 0, 0)] + lv[mapd(0, 1, 0, 0)]) *
        (lv[mapd(0, 0, 0, 1)] + lv[mapd(1, 0, 0, 1)])) -
        ((lv[mapd(-1, 0, 0, 0)] + lv[mapd(-1, 1, 0, 0)]) *
            (lv[mapd(-1, 0, 0, 1)] + lv[mapd(0, 0, 0, 1)])))
        + parameters.solver.gamma * ((fabs(lv[mapd(0, 0, 0, 0)] + lv[mapd(0, 1, 0, 0)]) *
            (lv[mapd(0, 0, 0, 1)] - lv[mapd(1, 0, 0, 1)])) -
            (fabs(lv[mapd(-1, 0, 0, 0)] + lv[mapd(-1, 1, 0, 0)]) *
                (lv[mapd(-1, 0, 0, 1)] - lv[mapd(0, 0, 0, 1)])))
        ) / lm[mapd(0, 0, 0, 0)];
#endif

    const RealType hxShort = 0.5 * lm[mapd(0, 0, 0, 0)];                           // Distance of corner points in x-direction from center v-value
    const RealType hxLong0 = 0.5 * (lm[mapd(0, 0, 0, 0)] + lm[mapd(-1, 0, 0, 0)]); // Distance between center and west
                                                                                   // v-value
    const RealType hxLong1 = 0.5 * (lm[mapd(0, 0, 0, 0)] + lm[mapd(1, 0, 0, 0)]);  // Distance between center and east
                                                                                   // v-value
    const RealType hyShort = 0.5 * lm[mapd(0, 0, 0, 1)];                           // Distance of center u-value from upper edge of cell
    const RealType hyLong  = 0.5 * (lm[mapd(0, 0, 0, 1)] + lm[mapd(0, 1, 0, 1)]);  // Distance of north and center u-value

    const RealType u00 = lv[mapd(0, 0, 0, 0)];
    const RealType u01 = lv[mapd(0, 1, 0, 0)];
    const RealType v00 = lv[mapd(0, 0, 0, 1)];
    const RealType v10 = lv[mapd(1, 0, 0, 1)];

    const RealType uM10 = lv[mapd(-1, 0, 0, 0)];
    const RealType uM11 = lv[mapd(-1, 1, 0, 0)];
    const RealType vM10 = lv[mapd(-1, 0, 0, 1)];

    // This a central difference expression for the first-derivative. We therefore linearly interpolate u*v onto the
    // surface of the current cell (in 2D: upper left and upper right corner) and then take the central difference.
    const RealType secondOrder = (((hyLong - hyShort) / hyLong * u00 + hyShort / hyLong * u01) * ((hxLong1 - hxShort) / hxLong1 * v00 + hxShort / hxLong1 * v10)
                                  - ((hyLong - hyShort) / hyLong * uM10 + hyShort / hyLong * uM11) * ((hxLong0 - hxShort) / hxLong0 * v00 + hxShort / hxLong0 * vM10))
                                 / (2.0 * hxShort);

    // This is a forward-difference in donor-cell style. We apply donor cell and again interpolate the velocity values
    // (u-comp.) onto the surface of the cell. We then apply the standard donor cell scheme. This will, however, result
    // in non-equal mesh spacing evaluations (in case of stretched meshes).
    const RealType kr = (hyLong - hyShort) / hyLong * u00 + hyShort / hyLong * u01;
    const RealType kl = (hyLong - hyShort) / hyLong * uM10 + hyShort / hyLong * uM11;

    const RealType firstOrder = 1.0 / (4.0 * hxShort) * (kr * (v00 + v10) - kl * (vM10 + v00) + fabs(kr) * (v00 - v10) - fabs(kl) * (vM10 - v00));

    // Return linear combination of central and donor-cell difference
    const RealType tmp2 = (1.0 - parameters.solver.gamma) * secondOrder + parameters.solver.gamma * firstOrder;

#ifndef NDEBUG
    if (fabs(tmp1 - tmp2) > 1.0e-12) {
      throw std::runtime_error("Error in duvdx");
    }
#endif

    return tmp2;
  }

  // Evaluates first derivative w.r.t. y for u*v at location of u-component. For details on implementation, see duvdx.
  inline RealType duvdy(const RealType* const lv, const Parameters& parameters, const RealType* const lm) {
#ifndef NDEBUG
    const RealType tmp1 = 1.0 / 4.0 * ((((lv[mapd(0, 0, 0, 1)] + lv[mapd(1, 0, 0, 1)]) *
        (lv[mapd(0, 0, 0, 0)] + lv[mapd(0, 1, 0, 0)])) -
        ((lv[mapd(0, -1, 0, 1)] + lv[mapd(1, -1, 0, 1)]) *
            (lv[mapd(0, -1, 0, 0)] + lv[mapd(0, 0, 0, 0)]))) +
        parameters.solver.gamma * ((fabs(lv[mapd(0, 0, 0, 1)] + lv[mapd(1, 0, 0, 1)]) *
            (lv[mapd(0, 0, 0, 0)] - lv[mapd(0, 1, 0, 0)])) -
            (fabs(lv[mapd(0, -1, 0, 1)] + lv[mapd(1, -1, 0, 1)]) *
                (lv[mapd(0, -1, 0, 0)] - lv[mapd(0, 0, 0, 0)])))) /
        lm[mapd(0, 0, 0, 1)];
#endif

    const RealType hyShort = 0.5 * lm[mapd(0, 0, 0, 1)];                           // Distance of corner points in x-direction from center v-value
    const RealType hyLong0 = 0.5 * (lm[mapd(0, 0, 0, 1)] + lm[mapd(0, -1, 0, 1)]); // Distance between center and west
                                                                                   // v-value
    const RealType hyLong1 = 0.5 * (lm[mapd(0, 0, 0, 1)] + lm[mapd(0, 1, 0, 1)]);  // Distance between center and east
                                                                                   // v-value
    const RealType hxShort = 0.5 * lm[mapd(0, 0, 0, 0)];                           // Distance of center u-value from upper edge of cell
    const RealType hxLong  = 0.5 * (lm[mapd(0, 0, 0, 0)] + lm[mapd(1, 0, 0, 0)]);  // Distance of north and center u-value

    const RealType v00 = lv[mapd(0, 0, 0, 1)];
    const RealType v10 = lv[mapd(1, 0, 0, 1)];
    const RealType u00 = lv[mapd(0, 0, 0, 0)];
    const RealType u01 = lv[mapd(0, 1, 0, 0)];

    const RealType v0M1 = lv[mapd(0, -1, 0, 1)];
    const RealType v1M1 = lv[mapd(1, -1, 0, 1)];
    const RealType u0M1 = lv[mapd(0, -1, 0, 0)];

    const RealType secondOrder = (((hxLong - hxShort) / hxLong * v00 + hxShort / hxLong * v10) * ((hyLong1 - hyShort) / hyLong1 * u00 + hyShort / hyLong1 * u01)
                                  - ((hxLong - hxShort) / hxLong * v0M1 + hxShort / hxLong * v1M1) * ((hyLong0 - hyShort) / hyLong0 * u00 + hyShort / hyLong0 * u0M1))
                                 / (2.0 * hyShort);

    const RealType kr = (hxLong - hxShort) / hxLong * v00 + hxShort / hxLong * v10;
    const RealType kl = (hxLong - hxShort) / hxLong * v0M1 + hxShort / hxLong * v1M1;

    const RealType firstOrder = 1.0 / (4.0 * hyShort) * (kr * (u00 + u01) - kl * (u0M1 + u00) + fabs(kr) * (u00 - u01) - fabs(kl) * (u0M1 - u00));

    const RealType tmp2 = (1.0 - parameters.solver.gamma) * secondOrder + parameters.solver.gamma * firstOrder;

#ifndef NDEBUG
    if (fabs(tmp1 - tmp2) > 1.0e-12) {
      throw std::runtime_error("Error in duvdy");
    }
#endif

    return tmp2;
  }

  // Evaluates first derivative w.r.t. x for u*w at location of w-component. For details on implementation, see duvdx.
  inline RealType duwdx(const RealType* const lv, const Parameters& parameters, const RealType* const lm) {
#ifndef NDEBUG
    const RealType tmp1 = 1.0 / 4.0 * ((((lv[mapd(0, 0, 0, 0)] + lv[mapd(0, 0, 1, 0)]) *
        (lv[mapd(0, 0, 0, 2)] + lv[mapd(1, 0, 0, 2)])) -
        ((lv[mapd(-1, 0, 0, 0)] + lv[mapd(-1, 0, 1, 0)]) *
            (lv[mapd(-1, 0, 0, 2)] + lv[mapd(0, 0, 0, 2)]))) +
        parameters.solver.gamma * ((fabs(lv[mapd(0, 0, 0, 0)] + lv[mapd(0, 0, 1, 0)]) *
            (lv[mapd(0, 0, 0, 2)] - lv[mapd(1, 0, 0, 2)])) -
            (fabs(lv[mapd(-1, 0, 0, 0)] + lv[mapd(-1, 0, 1, 0)]) *
                (lv[mapd(-1, 0, 0, 2)] - lv[mapd(0, 0, 0, 2)])))) /
        lm[mapd(0, 0, 0, 0)];
#endif

    const RealType hxShort = 0.5 * lm[mapd(0, 0, 0, 0)];                           // Distance of corner points in x-direction from center v-value
    const RealType hxLong0 = 0.5 * (lm[mapd(0, 0, 0, 0)] + lm[mapd(-1, 0, 0, 0)]); // Distance between center and west
                                                                                   // v-value
    const RealType hxLong1 = 0.5 * (lm[mapd(0, 0, 0, 0)] + lm[mapd(1, 0, 0, 0)]);  // Distance between center and east
                                                                                   // v-value
    const RealType hzShort = 0.5 * lm[mapd(0, 0, 0, 2)];                           // Distance of center u-value from upper edge of cell
    const RealType hzLong  = 0.5 * (lm[mapd(0, 0, 0, 2)] + lm[mapd(0, 0, 1, 2)]);  // Distance of north and center u-value

    const RealType u00 = lv[mapd(0, 0, 0, 0)];
    const RealType u01 = lv[mapd(0, 0, 1, 0)];
    const RealType w00 = lv[mapd(0, 0, 0, 2)];
    const RealType w10 = lv[mapd(1, 0, 0, 2)];

    const RealType uM10 = lv[mapd(-1, 0, 0, 0)];
    const RealType uM11 = lv[mapd(-1, 0, 1, 0)];
    const RealType wM10 = lv[mapd(-1, 0, 0, 2)];

    const RealType secondOrder = (((hzLong - hzShort) / hzLong * u00 + hzShort / hzLong * u01) * ((hxLong1 - hxShort) / hxLong1 * w00 + hxShort / hxLong1 * w10)
                                  - ((hzLong - hzShort) / hzLong * uM10 + hzShort / hzLong * uM11) * ((hxLong0 - hxShort) / hxLong0 * w00 + hxShort / hxLong0 * wM10))
                                 / (2.0 * hxShort);

    const RealType kr = (hzLong - hzShort) / hzLong * u00 + hzShort / hzLong * u01;
    const RealType kl = (hzLong - hzShort) / hzLong * uM10 + hzShort / hzLong * uM11;

    const RealType firstOrder = 1.0 / (4.0 * hxShort) * (kr * (w00 + w10) - kl * (wM10 + w00) + fabs(kr) * (w00 - w10) - fabs(kl) * (wM10 - w00));

    const RealType tmp2 = (1.0 - parameters.solver.gamma) * secondOrder + parameters.solver.gamma * firstOrder;

#ifndef NDEBUG
    if (fabs(tmp1 - tmp2) > 1.0e-12) {
      throw std::runtime_error("Error in duwdx");
    }
#endif

    return tmp2;
  }

  // Evaluates first derivative w.r.t. z for u*w at location of u-component. For details on implementation, see duvdx.
  inline RealType duwdz(const RealType* const lv, const Parameters& parameters, const RealType* const lm) {
#ifndef NDEBUG
    const RealType tmp1 = 1.0 / 4.0 * ((((lv[mapd(0, 0, 0, 2)] + lv[mapd(1, 0, 0, 2)]) *
        (lv[mapd(0, 0, 0, 0)] + lv[mapd(0, 0, 1, 0)])) -
        ((lv[mapd(0, 0, -1, 2)] + lv[mapd(1, 0, -1, 2)]) *
            (lv[mapd(0, 0, -1, 0)] + lv[mapd(0, 0, 0, 0)]))) +
        parameters.solver.gamma * ((fabs(lv[mapd(0, 0, 0, 2)] + lv[mapd(1, 0, 0, 2)]) *
            (lv[mapd(0, 0, 0, 0)] - lv[mapd(0, 0, 1, 0)])) -
            (fabs(lv[mapd(0, 0, -1, 2)] + lv[mapd(1, 0, -1, 2)]) *
                (lv[mapd(0, 0, -1, 0)] - lv[mapd(0, 0, 0, 0)])))) /
        lm[mapd(0, 0, 0, 2)];
#endif

    const RealType hzShort = 0.5 * lm[mapd(0, 0, 0, 2)];                           // Distance of corner points in x-direction from center v-value
    const RealType hzLong0 = 0.5 * (lm[mapd(0, 0, 0, 2)] + lm[mapd(0, 0, -1, 2)]); // Distance between center and west
                                                                                   // v-value
    const RealType hzLong1 = 0.5 * (lm[mapd(0, 0, 0, 2)] + lm[mapd(0, 0, 1, 2)]);  // Distance between center and east
                                                                                   // v-value
    const RealType hxShort = 0.5 * lm[mapd(0, 0, 0, 0)];                           // Distance of center u-value from upper edge of cell
    const RealType hxLong  = 0.5 * (lm[mapd(0, 0, 0, 0)] + lm[mapd(1, 0, 0, 0)]);  // Distance of north and center u-value

    const RealType w00 = lv[mapd(0, 0, 0, 2)];
    const RealType w10 = lv[mapd(1, 0, 0, 2)];
    const RealType u00 = lv[mapd(0, 0, 0, 0)];
    const RealType u01 = lv[mapd(0, 0, 1, 0)];

    const RealType w0M1 = lv[mapd(0, 0, -1, 2)];
    const RealType w1M1 = lv[mapd(1, 0, -1, 2)];
    const RealType u0M1 = lv[mapd(0, 0, -1, 0)];

    const RealType secondOrder = (((hxLong - hxShort) / hxLong * w00 + hxShort / hxLong * w10) * ((hzLong1 - hzShort) / hzLong1 * u00 + hzShort / hzLong1 * u01)
                                  - ((hxLong - hxShort) / hxLong * w0M1 + hxShort / hxLong * w1M1) * ((hzLong0 - hzShort) / hzLong0 * u00 + hzShort / hzLong0 * u0M1))
                                 / (2.0 * hzShort);

    const RealType kr = (hxLong - hxShort) / hxLong * w00 + hxShort / hxLong * w10;
    const RealType kl = (hxLong - hxShort) / hxLong * w0M1 + hxShort / hxLong * w1M1;

    const RealType firstOrder = 1.0 / (4.0 * hzShort) * (kr * (u00 + u01) - kl * (u0M1 + u00) + fabs(kr) * (u00 - u01) - fabs(kl) * (u0M1 - u00));

    const RealType tmp2 = (1.0 - parameters.solver.gamma) * secondOrder + parameters.solver.gamma * firstOrder;

#ifndef NDEBUG
    if (fabs(tmp1 - tmp2) > 1.0e-12) {
      throw std::runtime_error("Error in duwdz");
    }
#endif

    return tmp2;
  }

  // Evaluates first derivative w.r.t. y for v*w at location of w-component. For details on implementation, see duvdx.
  inline RealType dvwdy(const RealType* const lv, const Parameters& parameters, const RealType* const lm) {
#ifndef NDEBUG
    const RealType tmp1 = 1.0 / 4.0 * ((((lv[mapd(0, 0, 0, 1)] + lv[mapd(0, 0, 1, 1)]) *
        (lv[mapd(0, 0, 0, 2)] + lv[mapd(0, 1, 0, 2)])) -
        ((lv[mapd(0, -1, 0, 1)] + lv[mapd(0, -1, 1, 1)]) *
            (lv[mapd(0, -1, 0, 2)] + lv[mapd(0, 0, 0, 2)]))) +
        parameters.solver.gamma * ((fabs(lv[mapd(0, 0, 0, 1)] + lv[mapd(0, 0, 1, 1)]) *
            (lv[mapd(0, 0, 0, 2)] - lv[mapd(0, 1, 0, 2)])) -
            (fabs(lv[mapd(0, -1, 0, 1)] + lv[mapd(0, -1, 1, 1)]) *
                (lv[mapd(0, -1, 0, 2)] - lv[mapd(0, 0, 0, 2)])))) /
        lm[mapd(0, 0, 0, 1)];
#endif

    const RealType hyShort = 0.5 * lm[mapd(0, 0, 0, 1)];                           // Distance of corner points in x-direction from center v-value
    const RealType hyLong0 = 0.5 * (lm[mapd(0, 0, 0, 1)] + lm[mapd(0, -1, 0, 1)]); // Distance between center and west
                                                                                   // v-value
    const RealType hyLong1 = 0.5 * (lm[mapd(0, 0, 0, 1)] + lm[mapd(0, 1, 0, 1)]);  // Distance between center and east
                                                                                   // v-value
    const RealType hzShort = 0.5 * lm[mapd(0, 0, 0, 2)];                           // Distance of center u-value from upper edge of cell
    const RealType hzLong  = 0.5 * (lm[mapd(0, 0, 0, 2)] + lm[mapd(0, 0, 1, 2)]);  // Distance of north and center u-value

    const RealType v00 = lv[mapd(0, 0, 0, 1)];
    const RealType v01 = lv[mapd(0, 0, 1, 1)];
    const RealType w00 = lv[mapd(0, 0, 0, 2)];
    const RealType w10 = lv[mapd(0, 1, 0, 2)];

    const RealType vM10 = lv[mapd(0, -1, 0, 1)];
    const RealType vM11 = lv[mapd(0, -1, 1, 1)];
    const RealType wM10 = lv[mapd(0, -1, 0, 2)];

    const RealType secondOrder = (((hzLong - hzShort) / hzLong * v00 + hzShort / hzLong * v01) * ((hyLong1 - hyShort) / hyLong1 * w00 + hyShort / hyLong1 * w10)
                                  - ((hzLong - hzShort) / hzLong * vM10 + hzShort / hzLong * vM11) * ((hyLong0 - hyShort) / hyLong0 * w00 + hyShort / hyLong0 * wM10))
                                 / (2.0 * hyShort);

    const RealType kr = (hzLong - hzShort) / hzLong * v00 + hzShort / hzLong * v01;
    const RealType kl = (hzLong - hzShort) / hzLong * vM10 + hzShort / hzLong * vM11;

    const RealType firstOrder = 1.0 / (4.0 * hyShort) * (kr * (w00 + w10) - kl * (wM10 + w00) + fabs(kr) * (w00 - w10) - fabs(kl) * (wM10 - w00));

    const RealType tmp2 = (1.0 - parameters.solver.gamma) * secondOrder + parameters.solver.gamma * firstOrder;

#ifndef NDEBUG
    if (fabs(tmp1 - tmp2) > 1.0e-12) {
      throw std::runtime_error("Error in dvwdy");
    }
#endif

    return tmp2;
  }

  // Evaluates first derivative w.r.t. z for v*w at location of v-component. For details on implementation, see duvdx.
  inline RealType dvwdz(const RealType* const lv, const Parameters& parameters, const RealType* const lm) {
#ifndef NDEBUG
    const RealType tmp1 = 1.0 / 4.0 * ((((lv[mapd(0, 0, 0, 2)] + lv[mapd(0, 1, 0, 2)]) *
        (lv[mapd(0, 0, 0, 1)] + lv[mapd(0, 0, 1, 1)])) -
        ((lv[mapd(0, 0, -1, 2)] + lv[mapd(0, 1, -1, 2)]) *
            (lv[mapd(0, 0, -1, 1)] + lv[mapd(0, 0, 0, 1)]))) +
        parameters.solver.gamma * ((fabs(lv[mapd(0, 0, 0, 2)] + lv[mapd(0, 1, 0, 2)]) *
            (lv[mapd(0, 0, 0, 1)] - lv[mapd(0, 0, 1, 1)])) -
            (fabs(lv[mapd(0, 0, -1, 2)] + lv[mapd(0, 1, -1, 2)]) *
                (lv[mapd(0, 0, -1, 1)] - lv[mapd(0, 0, 0, 1)])))) /
        lm[mapd(0, 0, 0, 2)];
#endif

    const RealType hzShort = 0.5 * lm[mapd(0, 0, 0, 2)];                           // Distance of corner points in x-direction from center v-value
    const RealType hzLong0 = 0.5 * (lm[mapd(0, 0, 0, 2)] + lm[mapd(0, 0, -1, 2)]); // Distance between center and west
                                                                                   // v-value
    const RealType hzLong1 = 0.5 * (lm[mapd(0, 0, 0, 2)] + lm[mapd(0, 0, 1, 2)]);  // Distance between center and east
                                                                                   // v-value
    const RealType hyShort = 0.5 * lm[mapd(0, 0, 0, 1)];                           // Distance of center u-value from upper edge of cell
    const RealType hyLong  = 0.5 * (lm[mapd(0, 0, 0, 1)] + lm[mapd(0, 1, 0, 1)]);  // Distance of north and center u-value

    const RealType w00 = lv[mapd(0, 0, 0, 2)];
    const RealType w10 = lv[mapd(0, 1, 0, 2)];
    const RealType v00 = lv[mapd(0, 0, 0, 1)];
    const RealType v01 = lv[mapd(0, 0, 1, 1)];

    const RealType w0M1 = lv[mapd(0, 0, -1, 2)];
    const RealType w1M1 = lv[mapd(0, 1, -1, 2)];
    const RealType v0M1 = lv[mapd(0, 0, -1, 1)];

    const RealType secondOrder = (((hyLong - hyShort) / hyLong * w00 + hyShort / hyLong * w10) * ((hzLong1 - hzShort) / hzLong1 * v00 + hzShort / hzLong1 * v01)
                                  - ((hyLong - hyShort) / hyLong * w0M1 + hyShort / hyLong * w1M1) * ((hzLong0 - hzShort) / hzLong0 * v00 + hzShort / hzLong0 * v0M1))
                                 / (2.0 * hzShort);

    const RealType kr = (hyLong - hyShort) / hyLong * w00 + hyShort / hyLong * w10;
    const RealType kl = (hyLong - hyShort) / hyLong * w0M1 + hyShort / hyLong * w1M1;

    const RealType firstOrder = 1.0 / (4.0 * hzShort) * (kr * (v00 + v01) - kl * (v0M1 + v00) + fabs(kr) * (v00 - v01) - fabs(kl) * (v0M1 - v00));

    const RealType tmp2 = (1.0 - parameters.solver.gamma) * secondOrder + parameters.solver.gamma * firstOrder;

#ifndef NDEBUG
    if (fabs(tmp1 - tmp2) > 1.0e-12) {
      spdlog::error("{}, {}", tmp1, tmp2);
      throw std::runtime_error("Error in dvwdz");
    }
#endif

    return tmp2;
  }

  // First derivative of u*u w.r.t. x, evaluated at location of u-component.
  inline RealType du2dx(const RealType* const lv, const Parameters& parameters, const RealType* const lm) {
#ifndef NDEBUG
    const RealType tmp1 = 1.0 / 4.0 * ((((lv[mapd(0, 0, 0, 0)] + lv[mapd(1, 0, 0, 0)]) *
        (lv[mapd(0, 0, 0, 0)] + lv[mapd(1, 0, 0, 0)])) -
        ((lv[mapd(-1, 0, 0, 0)] + lv[mapd(0, 0, 0, 0)]) *
            (lv[mapd(-1, 0, 0, 0)] + lv[mapd(0, 0, 0, 0)]))) +
        parameters.solver.gamma * ((fabs(lv[mapd(0, 0, 0, 0)] + lv[mapd(1, 0, 0, 0)]) *
            (lv[mapd(0, 0, 0, 0)] - lv[mapd(1, 0, 0, 0)])) -
            (fabs(lv[mapd(-1, 0, 0, 0)] + lv[mapd(0, 0, 0, 0)]) *
                (lv[mapd(-1, 0, 0, 0)] - lv[mapd(0, 0, 0, 0)])))) /
        lm[mapd(0, 0, 0, 0)];
#endif

    const RealType dxShort = 0.5 * lm[mapd(0, 0, 0, 0)];
    // const RealType dxLong0 = 0.5*(lm[mapd(-1,0,0,0)] + lm[mapd(0,0,0,0)]);
    const RealType dxLong1 = 0.5 * (lm[mapd(0, 0, 0, 0)] + lm[mapd(1, 0, 0, 0)]);

    const RealType u0  = lv[mapd(0, 0, 0, 0)];
    const RealType uM1 = lv[mapd(-1, 0, 0, 0)];
    const RealType u1  = lv[mapd(1, 0, 0, 0)];

    // const RealType kr = (dxLong1 - dxShort) / dxLong1 * u0 + dxShort / dxLong1 * u1;
    // const RealType kl = (dxLong0 - dxShort) / dxLong0 * u0 + dxShort / dxLong0 * uM1;
    const RealType kr = (u0 + u1) / 2;
    const RealType kl = (u0 + uM1) / 2;

    // Central difference expression which is second-order accurate for uniform meshes. We interpolate u half-way
    // between neighboured u-component values and afterwards build the central difference for u*u.

    /*const RealType secondOrder = (((dxLong1 - dxShort) / dxLong1 * u0 + dxShort / dxLong1 * u1) * ((dxLong1 - dxShort)
       / dxLong1 * u0 + dxShort / dxLong1 * u1)
        - ((dxLong0 - dxShort) / dxLong0 * u0 + dxShort / dxLong0 * uM1) * ((dxLong0 - dxShort) / dxLong0 * u0 + dxShort
       / dxLong0 * uM1) ) / (2.0 * dxShort);*/

    const RealType secondOrder = ((u0 + u1) * (u0 + u1) - (u0 + uM1) * (u0 + uM1)) / (4 * dxLong1);

    // Donor-cell like derivative expression. We evaluate u half-way between neighboured u-components and use this as a
    // prediction of the transport direction.
    const RealType firstOrder = 1.0 / (4.0 * dxShort) * (kr * (u0 + u1) - kl * (uM1 + u0) + fabs(kr) * (u0 - u1) - fabs(kl) * (uM1 - u0));

    // Return linear combination of central- and upwind difference
    const RealType tmp2 = (1.0 - parameters.solver.gamma) * secondOrder + parameters.solver.gamma * firstOrder;

#ifndef NDEBUG
    if (fabs(tmp1 - tmp2) > 1.0e-12) {
      spdlog::error("{}, {}", tmp1, tmp2);
      throw std::runtime_error("Error in du2dx");
    }
#endif

    return tmp2;
  }

  // First derivative of v*v w.r.t. y, evaluated at location of v-component. For details, see du2dx.
  inline RealType dv2dy(const RealType* const lv, const Parameters& parameters, const RealType* const lm) {
#ifndef NDEBUG
    const RealType tmp1 = 1.0 / 4.0 * ((((lv[mapd(0, 0, 0, 1)] + lv[mapd(0, 1, 0, 1)]) *
        (lv[mapd(0, 0, 0, 1)] + lv[mapd(0, 1, 0, 1)])) -
        ((lv[mapd(0, -1, 0, 1)] + lv[mapd(0, 0, 0, 1)]) *
            (lv[mapd(0, -1, 0, 1)] + lv[mapd(0, 0, 0, 1)]))) +
        parameters.solver.gamma * ((fabs(lv[mapd(0, 0, 0, 1)] + lv[mapd(0, 1, 0, 1)]) *
            (lv[mapd(0, 0, 0, 1)] - lv[mapd(0, 1, 0, 1)])) -
            (fabs(lv[mapd(0, -1, 0, 1)] + lv[mapd(0, 0, 0, 1)]) *
                (lv[mapd(0, -1, 0, 1)] - lv[mapd(0, 0, 0, 1)])))) /
        lm[mapd(0, 0, 0, 1)];
#endif

    const RealType dyShort = 0.5 * lm[mapd(0, 0, 0, 1)];
    // const RealType dyLong0 = 0.5*(lm[mapd(0,-1,0,1)] + lm[mapd(0,0,0,1)]);
    const RealType dyLong1 = 0.5 * (lm[mapd(0, 0, 0, 1)] + lm[mapd(0, 1, 0, 1)]);

    const RealType v0  = lv[mapd(0, 0, 0, 1)];
    const RealType vM1 = lv[mapd(0, -1, 0, 1)];
    const RealType v1  = lv[mapd(0, 1, 0, 1)];

    // const RealType kr = (dyLong1 - dyShort) / dyLong1 * v0 + dyShort / dyLong1 * v1;
    // const RealType kl = (dyLong0 - dyShort) / dyLong0 * v0 + dyShort / dyLong0 * vM1;
    const RealType kr = (v0 + v1) / 2;
    const RealType kl = (v0 + vM1) / 2;

    /*const RealType secondOrder = (((dyLong1 - dyShort) / dyLong1 * v0 + dyShort / dyLong1 * v1) * ((dyLong1 - dyShort)
       / dyLong1 * v0 + dyShort / dyLong1 * v1)
        - ((dyLong0 - dyShort) / dyLong0 * v0 + dyShort / dyLong0 * vM1) * ((dyLong0 - dyShort) / dyLong0 * v0 + dyShort
       / dyLong0 * vM1) ) / (2.0 * dyShort);*/

    const RealType secondOrder = ((v0 + v1) * (v0 + v1) - (v0 + vM1) * (v0 + vM1)) / (4 * dyLong1);

    const RealType firstOrder = 1.0 / (4.0 * dyShort) * (kr * (v0 + v1) - kl * (vM1 + v0) + fabs(kr) * (v0 - v1) - fabs(kl) * (vM1 - v0));

    const RealType tmp2 = (1.0 - parameters.solver.gamma) * secondOrder + parameters.solver.gamma * firstOrder;

#ifndef NDEBUG
    if (fabs(tmp1 - tmp2) > 1.0e-12) {
      spdlog::error("{}, {}", tmp1, tmp2);
      throw std::runtime_error("Error in dv2dy");
    }
#endif

    return tmp2;
  }

  // First derivative of w*w w.r.t. z, evaluated at location of w-component. For details, see du2dx.
  inline RealType dw2dz(const RealType* const lv, const Parameters& parameters, const RealType* const lm) {
#ifndef NDEBUG
    const RealType tmp1 = 1.0 / 4.0 * ((((lv[mapd(0, 0, 0, 2)] + lv[mapd(0, 0, 1, 2)]) *
        (lv[mapd(0, 0, 0, 2)] + lv[mapd(0, 0, 1, 2)])) -
        ((lv[mapd(0, 0, -1, 2)] + lv[mapd(0, 0, 0, 2)]) *
            (lv[mapd(0, 0, -1, 2)] + lv[mapd(0, 0, 0, 2)]))) +
        parameters.solver.gamma * ((fabs(lv[mapd(0, 0, 0, 2)] + lv[mapd(0, 0, 1, 2)]) *
            (lv[mapd(0, 0, 0, 2)] - lv[mapd(0, 0, 1, 2)])) -
            (fabs(lv[mapd(0, 0, -1, 2)] + lv[mapd(0, 0, 0, 2)]) *
                (lv[mapd(0, 0, -1, 2)] - lv[mapd(0, 0, 0, 2)])))) /
        lm[mapd(0, 0, 0, 2)];
#endif

    const RealType dzShort = 0.5 * lm[mapd(0, 0, 0, 2)];
    // const RealType dzLong0 = 0.5 * (lm[mapd(0, 0, -1, 2)] + lm[mapd(0, 0, 0, 2)]);
    const RealType dzLong1 = 0.5 * (lm[mapd(0, 0, 0, 2)] + lm[mapd(0, 0, 1, 2)]);

    const RealType w0  = lv[mapd(0, 0, 0, 2)];
    const RealType wM1 = lv[mapd(0, 0, -1, 2)];
    const RealType w1  = lv[mapd(0, 0, 1, 2)];

    // const RealType kr = (dzLong1 - dzShort) / dzLong1 * w0 + dzShort / dzLong1 * w1;
    // const RealType kl = (dzLong0 - dzShort) / dzLong0 * w0 + dzShort / dzLong0 * wM1;
    const RealType kr = (w0 + w1) / 2;
    const RealType kl = (w0 + wM1) / 2;

    /*const RealType secondOrder = (((dzLong1 - dzShort) / dzLong1 * w0 + dzShort / dzLong1 * w1) * ((dzLong1 - dzShort)
       / dzLong1 * w0 + dzShort / dzLong1 * w1)
        - ((dzLong0 - dzShort) / dzLong0 * w0 + dzShort / dzLong0 * wM1) * ((dzLong0 - dzShort) / dzLong0 * w0 + dzShort
       / dzLong0 * wM1) ) / (2.0 * dzShort);*/

    const RealType secondOrder = ((w0 + w1) * (w0 + w1) - (w0 + wM1) * (w0 + wM1)) / (4 * dzLong1);

    const RealType firstOrder = 1.0 / (4.0 * dzShort) * (kr * (w0 + w1) - kl * (wM1 + w0) + fabs(kr) * (w0 - w1) - fabs(kl) * (wM1 - w0));

    const RealType tmp2 = (1.0 - parameters.solver.gamma) * secondOrder + parameters.solver.gamma * firstOrder;

#ifndef NDEBUG
    if (fabs(tmp1 - tmp2) > 1.0e-12) {
      spdlog::error("{}, {}", tmp1, tmp2);
      throw std::runtime_error("Error in dw2dz");
    }
#endif

    return tmp2;
  }

  inline RealType computeF2D(const RealType* const localVelocity, const RealType* const localMeshsize, const Parameters& parameters, RealType dt) {
    return localVelocity[mapd(0, 0, 0, 0)]
           + dt
               * (-du2dx(localVelocity, parameters, localMeshsize) - duvdy(localVelocity, parameters, localMeshsize) + 1 / parameters.flow.Re * (d2udx2(localVelocity, localMeshsize) + d2udy2(localVelocity, localMeshsize)));
  }

  inline RealType computeTurbulentF2D(
    const RealType* const localVelocity, const RealType* const localViscosity, const RealType* const localMeshsize, const Parameters& parameters, RealType dt
  ) {
    const RealType re = parameters.flow.Re;
    return localVelocity[mapd(0, 0, 0, 0)]
           + dt
               * (-du2dx(localVelocity, parameters, localMeshsize) - duvdy(localVelocity, parameters, localMeshsize) + 2 * dVdxdudx(localVelocity, localViscosity, localMeshsize, re) + dVdydudydvdx(localVelocity, localViscosity, localMeshsize, re));
  }

  inline RealType computeF3D(const RealType* const localVelocity, const RealType* const localMeshsize, const Parameters& parameters, RealType dt) {
    return localVelocity[mapd(0, 0, 0, 0)]
        + dt * (
            -du2dx(localVelocity, parameters, localMeshsize) 
            - duvdy(localVelocity, parameters, localMeshsize) 
            - duwdz(localVelocity, parameters, localMeshsize) 
            + 1 / parameters.flow.Re * (
              d2udx2(localVelocity, localMeshsize) 
              + d2udy2(localVelocity, localMeshsize) 
              + d2udz2(localVelocity, localMeshsize)
            )
            );
  }

  inline RealType computeTurbulentF3D(
    const RealType* const localVelocity, const RealType* const localViscosity, const RealType* const localMeshsize, const Parameters& parameters, RealType dt
  ) {
    const RealType re = parameters.flow.Re;
    return localVelocity[mapd(0, 0, 0, 0)]
           + dt
               * (-du2dx(localVelocity, parameters, localMeshsize) - duvdy(localVelocity, parameters, localMeshsize) - duwdz(localVelocity, parameters, localMeshsize) + 2 * dVdxdudx(localVelocity, localViscosity, localMeshsize, re) + dVdydudydvdx(localVelocity, localViscosity, localMeshsize, re) + dVdzdudzdwdx(localVelocity, localViscosity, localMeshsize, re));
  }

  inline RealType computeG2D(const RealType* const localVelocity, const RealType* const localMeshsize, const Parameters& parameters, RealType dt) {
    return localVelocity[mapd(0, 0, 0, 1)]
           + dt
               * (-duvdx(localVelocity, parameters, localMeshsize) - dv2dy(localVelocity, parameters, localMeshsize) + 1 / parameters.flow.Re * (d2vdx2(localVelocity, localMeshsize) + d2vdy2(localVelocity, localMeshsize)));
  }

  inline RealType computeTurbulentG2D(
    const RealType* const localVelocity, const RealType* const localViscosity, const RealType* const localMeshsize, const Parameters& parameters, RealType dt
  ) {
    const RealType re = parameters.flow.Re;
    return localVelocity[mapd(0, 0, 0, 1)]
           + dt
               * (-duvdx(localVelocity, parameters, localMeshsize) - dv2dy(localVelocity, parameters, localMeshsize) + dVdxdvdxdudy(localVelocity, localViscosity, localMeshsize, re) + 2 * dVdydvdy(localVelocity, localViscosity, localMeshsize, re));
  }

  inline RealType computeG3D(const RealType* const localVelocity, const RealType* const localMeshsize, const Parameters& parameters, RealType dt) {
    return localVelocity[mapd(0, 0, 0, 1)]
        + dt * (
          -dv2dy(localVelocity, parameters, localMeshsize) 
          - duvdx(localVelocity, parameters, localMeshsize)
          - dvwdz(localVelocity, parameters, localMeshsize) 
          + 1 / parameters.flow.Re * (
            d2vdx2(localVelocity, localMeshsize)
            + d2vdy2(localVelocity, localMeshsize)
            + d2vdz2(localVelocity, localMeshsize)
          )
          );
  }

  inline RealType computeTurbulentG3D(
    const RealType* const localVelocity, const RealType* const localViscosity, const RealType* const localMeshsize, const Parameters& parameters, RealType dt
  ) {
    const RealType re = parameters.flow.Re;
    return localVelocity[mapd(0, 0, 0, 1)]
           + dt
               * (-dv2dy(localVelocity, parameters, localMeshsize) - duvdx(localVelocity, parameters, localMeshsize) - dvwdz(localVelocity, parameters, localMeshsize) + dVdxdvdxdudy(localVelocity, localViscosity, localMeshsize, re) + 2 * dVdydvdy(localVelocity, localViscosity, localMeshsize, re) + dVdzdvdzdwdy(localVelocity, localViscosity, localMeshsize, re));
  }

  inline RealType computeH3D(const RealType* const localVelocity, const RealType* const localMeshsize, const Parameters& parameters, RealType dt) {
    return localVelocity[mapd(0, 0, 0, 2)]
        + dt * (
          -dw2dz(localVelocity, parameters, localMeshsize) 
          - duwdx(localVelocity, parameters, localMeshsize)
          - dvwdy(localVelocity, parameters, localMeshsize) 
          + 1 / parameters.flow.Re * (
            d2wdx2(localVelocity, localMeshsize)
            + d2wdy2(localVelocity, localMeshsize)
            + d2wdz2(localVelocity, localMeshsize)
          )
          );
  }

  inline RealType computeTurbulentH3D(
    const RealType* const localVelocity, const RealType* const localViscosity, const RealType* const localMeshsize, const Parameters& parameters, RealType dt
  ) {
    const RealType re = parameters.flow.Re;
    return localVelocity[mapd(0, 0, 0, 2)]
           + dt
               * (-dw2dz(localVelocity, parameters, localMeshsize) - duwdx(localVelocity, parameters, localMeshsize) - dvwdy(localVelocity, parameters, localMeshsize) + (dVdxdwdxdudz(localVelocity, localViscosity, localMeshsize, re) + dVdydwdydvdz(localVelocity, localViscosity, localMeshsize, re) + 2 * dVdzdwdz(localVelocity, localViscosity, localMeshsize, re)));
  }

} // namespace Stencils
