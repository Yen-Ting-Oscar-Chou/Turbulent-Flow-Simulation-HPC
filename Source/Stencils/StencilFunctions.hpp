#pragma once

#include "StdAfx.hpp"

#include "Definitions.hpp"
#include "Parameters.hpp"
#include "TurbulentFlowField.hpp"

namespace Stencils {

  inline RealType computeDistance(RealType* coords1, RealType* coords2, const Parameters& parameters_) {
    RealType result = 0.0;
    for (int i = 0; i < parameters_.geometry.dim; i++) {
      result += (coords1[i] - coords2[i]) * (coords1[i] - coords2[i]);
    }
    return sqrt(result);
  }

  inline void computeGlobalCoordinates(RealType* coords, const Parameters& parameters_, int i, int j) {
    coords[0] = parameters_.meshsize->getPosX(i, j) + 0.5 * parameters_.meshsize->getDx(i, j);
    coords[1] = parameters_.meshsize->getPosY(i, j) + 0.5 * parameters_.meshsize->getDy(i, j);
  }

  inline void computeGlobalXCoordinate(RealType& x, const Parameters& parameters, int i, int j, int k = 0) {
    x = parameters.meshsize->getPosX(i, j) + 0.5 * parameters.meshsize->getDx(i, j, k);
  }

  inline void computeGlobalCoordinates(RealType* coords, const Parameters& parameters_, int i, int j, int k) {
    coords[0] = parameters_.meshsize->getPosX(i, j, k) + 0.5 * parameters_.meshsize->getDx(i, j, k);
    coords[1] = parameters_.meshsize->getPosY(i, j, k) + 0.5 * parameters_.meshsize->getDy(i, j, k);
    coords[2] = parameters_.meshsize->getPosZ(i, j, k) + 0.5 * parameters_.meshsize->getDz(i, j, k);
  }

  inline RealType firstDerivative(VectorField& velocity, const Parameters& parameters, const COMP comp, const DERIV deriv, int i, int j, int k = 0) {
    const bool isX = DERIVX == deriv;
    const bool isY = DERIVY == deriv;
    const bool isZ = DERIVZ == deriv;

    RealType dp;
    RealType dm;
    switch (deriv) {
    case DERIVX:
      dm = parameters.meshsize->getDx(i - 1, j, k);
      dp = parameters.meshsize->getDx(i, j, k);
      break;
    case DERIVY:
      dm = parameters.meshsize->getDy(i, j - 1, k);
      dp = parameters.meshsize->getDy(i, j, k);
      break;
    case DERIVZ:
      dm = parameters.meshsize->getDz(i, j, k - 1);
      dp = parameters.meshsize->getDz(i, j, k);
      break;
    }

    const RealType lvPlus  = velocity.getVectorElement(i + isX, j + isY, k + isZ, comp);
    const RealType lvThis  = velocity.getVectorElement(i, j, k, comp);
    const RealType lvMinus = velocity.getVectorElement(i - isX, j - isY, k - isZ, comp);

    return 1 / (2 * (dp * dm)) * (dm * lvPlus + (dp - dm) * lvThis - dp * lvMinus);
  }

  inline RealType secondDerivative(VectorField& velocity, const Parameters& parameters, const COMP comp, const DERIV deriv, int i, int j, int k = 0) {
    const bool isX = DERIVX == deriv;
    const bool isY = DERIVY == deriv;
    const bool isZ = DERIVZ == deriv;

    RealType dp;
    RealType dm;
    switch (deriv) {
    case DERIVX:
      dm = parameters.meshsize->getDx(i - 1, j, k);
      dp = parameters.meshsize->getDx(i, j, k);
      break;
    case DERIVY:
      dm = parameters.meshsize->getDy(i, j - 1, k);
      dp = parameters.meshsize->getDy(i, j, k);
      break;
    case DERIVZ:
      dm = parameters.meshsize->getDz(i, j, k - 1);
      dp = parameters.meshsize->getDz(i, j, k);
      break;
    }

    const RealType fac     = 2 / (dm * dp * dp + dp * dm * dm);
    const RealType lvPlus  = velocity.getVectorElement(i + isX, j + isY, k + isZ, comp);
    const RealType lvThis  = velocity.getVectorElement(i, j, k, comp);
    const RealType lvMinus = velocity.getVectorElement(i - isX, j - isY, k - isZ, comp);
    return fac * (dm * (lvPlus - lvThis) + dp * (lvMinus - lvThis));
  }

  // Derivative functions. They are applied to a cube of 3x3x3 cells. lv stands for the local velocity, lm represents
  // the local mesh sizes dudx <-> first derivative of u-component of velocity field w.r.t. x-direction.
  inline RealType dudx(VectorField& velocity, const Parameters& parameters, int i, int j, int k = 0) { return firstDerivative(velocity, parameters, COMPX, DERIVX, i, j, k); }

  inline RealType dudy(VectorField& velocity, const Parameters& parameters, int i, int j, int k = 0) { return firstDerivative(velocity, parameters, COMPX, DERIVY, i, j, k); }

  inline RealType dudz(VectorField& velocity, const Parameters& parameters, int i, int j, int k = 0) { return firstDerivative(velocity, parameters, COMPX, DERIVZ, i, j, k); }

  inline RealType dvdx(VectorField& velocity, const Parameters& parameters, int i, int j, int k = 0) { return firstDerivative(velocity, parameters, COMPY, DERIVX, i, j, k); }

  inline RealType dvdy(VectorField& velocity, const Parameters& parameters, int i, int j, int k = 0) { return firstDerivative(velocity, parameters, COMPY, DERIVY, i, j, k); }

  inline RealType dvdz(VectorField& velocity, const Parameters& parameters, int i, int j, int k = 0) { return firstDerivative(velocity, parameters, COMPY, DERIVZ, i, j, k); }

  inline RealType dwdx(VectorField& velocity, const Parameters& parameters, int i, int j, int k = 0) { return firstDerivative(velocity, parameters, COMPZ, DERIVX, i, j, k); }

  inline RealType dwdy(VectorField& velocity, const Parameters& parameters, int i, int j, int k = 0) { return firstDerivative(velocity, parameters, COMPZ, DERIVY, i, j, k); }

  inline RealType dwdz(VectorField& velocity, const Parameters& parameters, int i, int j, int k = 0) { return firstDerivative(velocity, parameters, COMPZ, DERIVZ, i, j, k); }

  inline RealType interpolateViscosity(const RealType visc_bl, const RealType visc_br, const RealType visc_tl, const RealType visc_tr) {
    return 0.25 * (visc_bl + visc_br + visc_tl + visc_tr);
  }

  inline RealType viscositySingleDerivative(VectorField& velocity, ScalarField& viscosity, const Parameters& parameters, const RealType re, const COMP comp, int i, int j, int k) {
    const bool  isX   = COMPX == comp;
    const bool  isY   = COMPY == comp;
    const bool  isZ   = COMPZ == comp;
    const DERIV deriv = static_cast<DERIV>(comp);

    RealType dm;
    RealType dp;
    switch (deriv) {
    case DERIVX:
      dm = parameters.meshsize->getDx(i - 1, j, k);
      dp = parameters.meshsize->getDx(i, j, k);
      break;
    case DERIVY:
      dm = parameters.meshsize->getDy(i, j - 1, k);
      dp = parameters.meshsize->getDy(i, j, k);
      break;
    case DERIVZ:
      dm = parameters.meshsize->getDz(i, j, k - 1);
      dp = parameters.meshsize->getDz(i, j, k);
      break;
    }

    const RealType d = (dm + dp) / 2;

    const RealType lVelPlus  = velocity.getVectorElement(i + isX, j + isY, k + isZ, comp);
    const RealType lVelThis  = velocity.getVectorElement(i, j, k, comp);
    const RealType lVelMinus = velocity.getVectorElement(i - isX, j - isY, k - isZ, comp);
    const RealType forward   = (lVelPlus - lVelThis) / dp;
    const RealType backward  = (lVelThis - lVelMinus) / dm;

    const RealType vis_plus = viscosity.getScalar(i + isX * 1, j + isY * 1, k + isZ * 1) + (1 / re);
    const RealType vis_this = viscosity.getScalar(i, j, k) + (1 / re);

    return 1 / d * (vis_plus * forward - vis_this * backward);
  }

  inline RealType viscosityDoubleDerivative(
    VectorField& velocity, ScalarField& viscosity, const Parameters& parameters, const RealType re, const COMP comp1, const COMP comp2, int i, int j, int k
  ) {
    // first component is the outer derivative and first inner derivative
    const bool isX1 = COMPX == comp1;
    const bool isY1 = COMPY == comp1;
    const bool isZ1 = COMPZ == comp1;
    // second component is the second inner derivative
    const bool  isX2   = COMPX == comp2;
    const bool  isY2   = COMPY == comp2;
    const bool  isZ2   = COMPZ == comp2;
    const DERIV deriv1 = static_cast<DERIV>(comp1);
    const DERIV deriv2 = static_cast<DERIV>(comp2);

    RealType dm1;
    RealType dm2;
    RealType dp1;
    RealType dp2;

    switch (deriv1) {
    case DERIVX:
      dm1 = parameters.meshsize->getDx(i - isX1, j - isY1, k - isZ1);
      dp1 = parameters.meshsize->getDx(i, j, k);
      break;
    case DERIVY:
      dm1 = parameters.meshsize->getDy(i - isX1, j - isY1, k - isZ1);
      dp1 = parameters.meshsize->getDy(i, j, k);
      break;
    case DERIVZ:
      dm1 = parameters.meshsize->getDz(i - isX1, j - isY1, k - isZ1);
      dp1 = parameters.meshsize->getDz(i, j, k);
      break;
    }

    switch (deriv2) {
    case DERIVX:
      dm2 = parameters.meshsize->getDx(i - isX2, j - isY2, k - isZ2);
      dp2 = parameters.meshsize->getDx(i, j, k);
      break;
    case DERIVY:
      dm2 = parameters.meshsize->getDy(i - isX2, j - isY2, k - isZ2);
      dp2 = parameters.meshsize->getDy(i, j, k);
      break;
    case DERIVZ:
      dm2 = parameters.meshsize->getDz(i - isX2, j - isY2, k - isZ2);
      dp2 = parameters.meshsize->getDz(i, j, k);
      break;
    }


    const RealType d = (dm1 + dp1) / 2; // stretched discretization length in first component

    // comments are for y, x
    const RealType viscBottomLeft1  = viscosity.getScalar(i, j, k);                                           // (i,j)
    const RealType viscBottomRight1 = viscosity.getScalar(i + isX2, j + isY2, k + isZ2);                      // (i+1, j)
    const RealType viscTopLeft1     = viscosity.getScalar(i + isX1, j + isY1, k + isZ1);                      // (i, j+1)
    const RealType viscTopRight1    = viscosity.getScalar(i + isX1 + isX2, j + isY1 + isY2, k + isZ1 + isZ2); // (i+1, j+1)
    const RealType visc1            = interpolateViscosity(viscBottomLeft1, viscBottomRight1, viscTopLeft1, viscTopRight1) + (1 / re);

    const RealType viscBottomLeft2  = viscosity.getScalar(i - isX1, j - isY1, k - isZ1);                      // (i, j-1)
    const RealType viscBottomRight2 = viscosity.getScalar(i - isX1 + isX2, j - isY1 + isY2, k - isZ1 + isZ2); // (i+1, j-1)
    const RealType viscTopLeft2     = viscosity.getScalar(i, j, k);                                           // (i, j)
    const RealType viscTopRight2    = viscosity.getScalar(i + isX2, j + isY2, k + isZ2);                      // (i+1, j)
    const RealType visc2            = interpolateViscosity(viscBottomLeft2, viscBottomRight2, viscTopLeft2, viscTopRight2) + (1 / re);

    const RealType velPlus1       = velocity.getVectorElement(i + isX1, j + isY1, k + isZ1, comp2);
    const RealType velPlus2       = velocity.getVectorElement(i + isX2, j + isY2, k + isZ2, comp1);
    const RealType velThis1       = velocity.getVectorElement(i, j, k, comp2);
    const RealType velThis2       = velocity.getVectorElement(i, j, k, comp1);
    const RealType velMinus1      = velocity.getVectorElement(i - isX1, j - isY1, k - isZ1, comp2);
    const RealType velPlus2Minus1 = velocity.getVectorElement(i + isX2 - isX1, j + isY2 - isY1, k + isZ2 - isZ1, comp1);
    const RealType velMinus1Comp1 = velocity.getVectorElement(i - isX1, j - isY1, k - isZ1, comp1);

    const RealType forward1  = (velPlus1 - velThis1) / dp1;
    const RealType forward2  = (velPlus2 - velThis2) / dp2;
    const RealType backward1 = (velThis1 - velMinus1) / dm1;
    const RealType backward2 = (velPlus2Minus1 - velMinus1Comp1) / dm2;

    return 1 / d * (visc1 * (forward1 + forward2) - visc2 * (backward1 + backward2));
  }

  // d/dx * (v * du/dx)
  inline RealType dVdxdudx(VectorField& velocity, ScalarField& viscosity, const Parameters& parameters, const RealType re, int i, int j, int k = 0) {
    return viscositySingleDerivative(velocity, viscosity, parameters, re, COMPX, i, j, k);
  }

  // d/dy * (v * dv/dy)
  inline RealType dVdydvdy(VectorField& velocity, ScalarField& viscosity, const Parameters& parameters, const RealType re, int i, int j, int k = 0) {
    return viscositySingleDerivative(velocity, viscosity, parameters, re, COMPY, i, j, k);
  }

  // d/dz * (v * dw/dz)
  inline RealType dVdzdwdz(VectorField& velocity, ScalarField& viscosity, const Parameters& parameters, const RealType re, int i, int j, int k) {
    return viscositySingleDerivative(velocity, viscosity, parameters, re, COMPZ, i, j, k);
  }

  // d/dy * (v * (du/dy + dv/dx))
  inline RealType dVdydudydvdx(VectorField& velocity, ScalarField& viscosity, const Parameters& parameters, const RealType re, int i, int j, int k = 0) {
    return viscosityDoubleDerivative(velocity, viscosity, parameters, re, COMPY, COMPX, i, j, k);
  }

  // d/dz * (v * (du/dz + dw/dx))
  inline RealType dVdzdudzdwdx(VectorField& velocity, ScalarField& viscosity, const Parameters& parameters, const RealType re, int i, int j, int k) {
    return viscosityDoubleDerivative(velocity, viscosity, parameters, re, COMPZ, COMPX, i, j, k);
  }

  // d/dx * (v * (dv/dx + du/dy))
  inline RealType dVdxdvdxdudy(VectorField& velocity, ScalarField& viscosity, const Parameters& parameters, const RealType re, int i, int j, int k = 0) {
    return viscosityDoubleDerivative(velocity, viscosity, parameters, re, COMPX, COMPY, i, j, k);
  }

  // d/dz * (v * (dv/dz + dw/dy))
  inline RealType dVdzdvdzdwdy(VectorField& velocity, ScalarField& viscosity, const Parameters& parameters, const RealType re, int i, int j, int k) {
    return viscosityDoubleDerivative(velocity, viscosity, parameters, re, COMPZ, COMPY, i, j, k);
  }

  // d/dx * (v * (dw/dx + du/dz))
  inline RealType dVdxdwdxdudz(VectorField& velocity, ScalarField& viscosity, const Parameters& parameters, const RealType re, int i, int j, int k) {
    return viscosityDoubleDerivative(velocity, viscosity, parameters, re, COMPX, COMPZ, i, j, k);
  }

  // d/dy * (v * (dw/dy + dv/dz))
  inline RealType dVdydwdydvdz(VectorField& velocity, ScalarField& viscosity, const Parameters& parameters, const RealType re, int i, int j, int k) {
    return viscosityDoubleDerivative(velocity, viscosity, parameters, re, COMPY, COMPZ, i, j, k);
  }

  // TODO WS1: Second derivatives
  // Second derivative of u w.r.t. x-direction
  inline RealType d2udx2(VectorField& velocity, const Parameters& parameters, int i, int j, int k = 0) { return secondDerivative(velocity, parameters, COMPX, DERIVX, i, j, k); }

  // Second derivative of u w.r.t. y-direction
  inline RealType d2udy2(VectorField& velocity, const Parameters& parameters, int i, int j, int k = 0) { return secondDerivative(velocity, parameters, COMPX, DERIVY, i, j, k); }

  // Second derivative of u w.r.t. z-direction
  inline RealType d2udz2(VectorField& velocity, const Parameters& parameters, int i, int j, int k = 0) { return secondDerivative(velocity, parameters, COMPX, DERIVZ, i, j, k); }

  // Second derivative of v w.r.t. x-direction
  inline RealType d2vdx2(VectorField& velocity, const Parameters& parameters, int i, int j, int k = 0) { return secondDerivative(velocity, parameters, COMPY, DERIVX, i, j, k); }

  // Second derivative of v w.r.t. y-direction
  inline RealType d2vdy2(VectorField& velocity, const Parameters& parameters, int i, int j, int k = 0) { return secondDerivative(velocity, parameters, COMPY, DERIVY, i, j, k); }

  // Second derivative of v w.r.t. z-direction
  inline RealType d2vdz2(VectorField& velocity, const Parameters& parameters, int i, int j, int k = 0) { return secondDerivative(velocity, parameters, COMPY, DERIVZ, i, j, k); }

  // Second derivative of w w.r.t. x-direction
  inline RealType d2wdx2(VectorField& velocity, const Parameters& parameters, int i, int j, int k = 0) { return secondDerivative(velocity, parameters, COMPZ, DERIVX, i, j, k); }

  // Second derivative of w w.r.t. y-direction
  inline RealType d2wdy2(VectorField& velocity, const Parameters& parameters, int i, int j, int k = 0) { return secondDerivative(velocity, parameters, COMPZ, DERIVY, i, j, k); }

  // Second derivative of w w.r.t. z-direction
  inline RealType d2wdz2(VectorField& velocity, const Parameters& parameters, int i, int j, int k = 0) { return secondDerivative(velocity, parameters, COMPZ, DERIVZ, i, j, k); }

  // First derivative of product (u*v), evaluated at the location of the v-component.
  inline RealType duvdx(VectorField& velocity, const Parameters& parameters, int i, int j, int k = 0) {
    const RealType mp0000  = parameters.meshsize->getDx(i, j, k);
    const RealType hxShort = 0.5 * mp0000;                                             // Distance of corner points in x-direction from center v-value
    const RealType hxLong0 = 0.5 * (mp0000 + parameters.meshsize->getDx(i - 1, j, k)); // Distance between center and west
                                                                                       // v-value
    const RealType hxLong1 = 0.5 * (mp0000 + parameters.meshsize->getDx(i + 1, j, k)); // Distance between center and east
    // v-value
    const RealType mp0001  = parameters.meshsize->getDy(i, j, k);
    const RealType hyShort = 0.5 * mp0001;                                             // Distance of center u-value from upper edge of cell
    const RealType hyLong  = 0.5 * (mp0001 + parameters.meshsize->getDy(i, j + 1, k)); // Distance of north and center u-value

    const RealType u00 = velocity.getVectorElement(i, j, k, 0);
    const RealType u01 = velocity.getVectorElement(i, j + 1, k, 0);
    const RealType v00 = velocity.getVectorElement(i, j, k, 1);
    const RealType v10 = velocity.getVectorElement(i + 1, j, k, 1);

    const RealType uM10 = velocity.getVectorElement(i - 1, j, k, 0);
    const RealType uM11 = velocity.getVectorElement(i - 1, j + 1, k, 0);
    const RealType vM10 = velocity.getVectorElement(i - 1, j, k, 1);

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
    return tmp2;
  }

  // Evaluates first derivative w.r.t. y for u*v at location of u-component. For details on implementation, see duvdx.
  inline RealType duvdy(VectorField& velocity, const Parameters& parameters, int i, int j, int k = 0) {
    const RealType mp0001  = parameters.meshsize->getDy(i, j, k);
    const RealType hyShort = 0.5 * mp0001;                                             // Distance of corner points in x-direction from center v-value
    const RealType hyLong0 = 0.5 * (mp0001 + parameters.meshsize->getDy(i, j - 1, k)); // Distance between center and west
                                                                                       // v-value
    const RealType hyLong1 = 0.5 * (mp0001 + parameters.meshsize->getDy(i, j + 1, k)); // Distance between center and east
                                                                                       // v-value
    const RealType mp0000  = parameters.meshsize->getDx(i, j, k);
    const RealType hxShort = 0.5 * mp0000;                                             // Distance of center u-value from upper edge of cell
    const RealType hxLong  = 0.5 * (mp0000 + parameters.meshsize->getDx(i + 1, j, k)); // Distance of north and center u-value

    const RealType v00 = velocity.getVectorElement(i, j, k, 1);
    const RealType v10 = velocity.getVectorElement(i + 1, j, k, 1);
    const RealType u00 = velocity.getVectorElement(i, j, k, 0);
    const RealType u01 = velocity.getVectorElement(i, j + 1, k, 0);

    const RealType v0M1 = velocity.getVectorElement(i, j - 1, k, 1);
    const RealType v1M1 = velocity.getVectorElement(i + 1, j - 1, k, 1);
    const RealType u0M1 = velocity.getVectorElement(i, j - 1, k, 0);

    const RealType secondOrder = (((hxLong - hxShort) / hxLong * v00 + hxShort / hxLong * v10) * ((hyLong1 - hyShort) / hyLong1 * u00 + hyShort / hyLong1 * u01)
                                  - ((hxLong - hxShort) / hxLong * v0M1 + hxShort / hxLong * v1M1) * ((hyLong0 - hyShort) / hyLong0 * u00 + hyShort / hyLong0 * u0M1))
                                 / (2.0 * hyShort);

    const RealType kr = (hxLong - hxShort) / hxLong * v00 + hxShort / hxLong * v10;
    const RealType kl = (hxLong - hxShort) / hxLong * v0M1 + hxShort / hxLong * v1M1;

    const RealType firstOrder = 1.0 / (4.0 * hyShort) * (kr * (u00 + u01) - kl * (u0M1 + u00) + fabs(kr) * (u00 - u01) - fabs(kl) * (u0M1 - u00));

    const RealType tmp2 = (1.0 - parameters.solver.gamma) * secondOrder + parameters.solver.gamma * firstOrder;

    return tmp2;
  }

  // Evaluates first derivative w.r.t. x for u*w at location of w-component. For details on implementation, see duvdx.
  inline RealType duwdx(VectorField& velocity, const Parameters& parameters, int i, int j, int k) {
    const RealType mp0000  = parameters.meshsize->getDx(i, j, k);
    const RealType hxShort = 0.5 * mp0000;                                             // Distance of corner points in x-direction from center v-value
    const RealType hxLong0 = 0.5 * (mp0000 + parameters.meshsize->getDx(i - 1, j, k)); // Distance between center and west
                                                                                       // v-value
    const RealType hxLong1 = 0.5 * (mp0000 + parameters.meshsize->getDx(i + 1, j, k)); // Distance between center and east
    // v-value
    const RealType mp0002  = parameters.meshsize->getDz(i, j, k);
    const RealType hzShort = 0.5 * mp0002;                                             // Distance of center u-value from upper edge of cell
    const RealType hzLong  = 0.5 * (mp0002 + parameters.meshsize->getDz(i, j, k + 1)); // Distance of north and center u-value

    const RealType u00 = velocity.getVectorElement(i, j, k, 0);
    const RealType u01 = velocity.getVectorElement(i, j, k + 1, 0);
    const RealType w00 = velocity.getVectorElement(i, j, k, 2);
    const RealType w10 = velocity.getVectorElement(i + 1, j, k, 2);

    const RealType uM10 = velocity.getVectorElement(i - 1, j, k, 0);
    const RealType uM11 = velocity.getVectorElement(i - 1, j, k + 1, 0);
    const RealType wM10 = velocity.getVectorElement(i - 1, j, k, 2);

    const RealType secondOrder = (((hzLong - hzShort) / hzLong * u00 + hzShort / hzLong * u01) * ((hxLong1 - hxShort) / hxLong1 * w00 + hxShort / hxLong1 * w10)
                                  - ((hzLong - hzShort) / hzLong * uM10 + hzShort / hzLong * uM11) * ((hxLong0 - hxShort) / hxLong0 * w00 + hxShort / hxLong0 * wM10))
                                 / (2.0 * hxShort);

    const RealType kr = (hzLong - hzShort) / hzLong * u00 + hzShort / hzLong * u01;
    const RealType kl = (hzLong - hzShort) / hzLong * uM10 + hzShort / hzLong * uM11;

    const RealType firstOrder = 1.0 / (4.0 * hxShort) * (kr * (w00 + w10) - kl * (wM10 + w00) + fabs(kr) * (w00 - w10) - fabs(kl) * (wM10 - w00));

    const RealType tmp2 = (1.0 - parameters.solver.gamma) * secondOrder + parameters.solver.gamma * firstOrder;

    return tmp2;
  }

  // Evaluates first derivative w.r.t. z for u*w at location of u-component. For details on implementation, see duvdx.
  inline RealType duwdz(VectorField& velocity, const Parameters& parameters, int i, int j, int k) {
    const RealType mp0002  = parameters.meshsize->getDz(i, j, k);
    const RealType hzShort = 0.5 * mp0002;                                             // Distance of corner points in x-direction from center v-value
    const RealType hzLong0 = 0.5 * (mp0002 + parameters.meshsize->getDz(i, j, k - 1)); // Distance between center and west
                                                                                       // v-value
    const RealType hzLong1 = 0.5 * (mp0002 + parameters.meshsize->getDz(i, j, k + 1)); // Distance between center and east
    // v-value
    const RealType mp0000  = parameters.meshsize->getDx(i, j, k);
    const RealType hxShort = 0.5 * mp0000;                                             // Distance of center u-value from upper edge of cell
    const RealType hxLong  = 0.5 * (mp0000 + parameters.meshsize->getDx(i + 1, j, k)); // Distance of north and center u-value

    const RealType w00 = velocity.getVectorElement(i, j, k, 2);
    const RealType w10 = velocity.getVectorElement(i + 1, j, k, 2);
    const RealType u00 = velocity.getVectorElement(i, j, k, 0);
    const RealType u01 = velocity.getVectorElement(i, j, k + 1, 0);

    const RealType w0M1 = velocity.getVectorElement(i, j, k - 1, 2);
    const RealType w1M1 = velocity.getVectorElement(i + 1, j, k - 1, 2);
    const RealType u0M1 = velocity.getVectorElement(i, j, k - 1, 0);

    const RealType secondOrder = (((hxLong - hxShort) / hxLong * w00 + hxShort / hxLong * w10) * ((hzLong1 - hzShort) / hzLong1 * u00 + hzShort / hzLong1 * u01)
                                  - ((hxLong - hxShort) / hxLong * w0M1 + hxShort / hxLong * w1M1) * ((hzLong0 - hzShort) / hzLong0 * u00 + hzShort / hzLong0 * u0M1))
                                 / (2.0 * hzShort);

    const RealType kr = (hxLong - hxShort) / hxLong * w00 + hxShort / hxLong * w10;
    const RealType kl = (hxLong - hxShort) / hxLong * w0M1 + hxShort / hxLong * w1M1;

    const RealType firstOrder = 1.0 / (4.0 * hzShort) * (kr * (u00 + u01) - kl * (u0M1 + u00) + fabs(kr) * (u00 - u01) - fabs(kl) * (u0M1 - u00));

    const RealType tmp2 = (1.0 - parameters.solver.gamma) * secondOrder + parameters.solver.gamma * firstOrder;

    return tmp2;
  }

  // Evaluates first derivative w.r.t. y for v*w at location of w-component. For details on implementation, see duvdx.
  inline RealType dvwdy(VectorField& velocity, const Parameters& parameters, int i, int j, int k) {
    const RealType mp0001  = parameters.meshsize->getDy(i, j, k);
    const RealType hyShort = 0.5 * mp0001;                                             // Distance of corner points in x-direction from center v-value
    const RealType hyLong0 = 0.5 * (mp0001 + parameters.meshsize->getDy(i, j - 1, k)); // Distance between center and west
                                                                                       // v-value
    const RealType hyLong1 = 0.5 * (mp0001 + parameters.meshsize->getDy(i, j + 1, k)); // Distance between center and east
                                                                                       // v-value
    const RealType mp0002  = parameters.meshsize->getDz(i, j, k);
    const RealType hzShort = 0.5 * mp0002;                                             // Distance of center u-value from upper edge of cell
    const RealType hzLong  = 0.5 * (mp0002 + parameters.meshsize->getDz(i, j, k + 1)); // Distance of north and center u-value

    const RealType v00 = velocity.getVectorElement(i, j, k, 1);
    const RealType v01 = velocity.getVectorElement(i, j, k + 1, 1);
    const RealType w00 = velocity.getVectorElement(i, j, k, 2);
    const RealType w10 = velocity.getVectorElement(i, j + 1, k, 2);

    const RealType vM10 = velocity.getVectorElement(i, j - 1, k, 1);
    const RealType vM11 = velocity.getVectorElement(i, j - 1, k + 1, 1);
    const RealType wM10 = velocity.getVectorElement(i, j - 1, k, 2);

    const RealType secondOrder = (((hzLong - hzShort) / hzLong * v00 + hzShort / hzLong * v01) * ((hyLong1 - hyShort) / hyLong1 * w00 + hyShort / hyLong1 * w10)
                                  - ((hzLong - hzShort) / hzLong * vM10 + hzShort / hzLong * vM11) * ((hyLong0 - hyShort) / hyLong0 * w00 + hyShort / hyLong0 * wM10))
                                 / (2.0 * hyShort);

    const RealType kr = (hzLong - hzShort) / hzLong * v00 + hzShort / hzLong * v01;
    const RealType kl = (hzLong - hzShort) / hzLong * vM10 + hzShort / hzLong * vM11;

    const RealType firstOrder = 1.0 / (4.0 * hyShort) * (kr * (w00 + w10) - kl * (wM10 + w00) + fabs(kr) * (w00 - w10) - fabs(kl) * (wM10 - w00));

    const RealType tmp2 = (1.0 - parameters.solver.gamma) * secondOrder + parameters.solver.gamma * firstOrder;

    return tmp2;
  }

  // Evaluates first derivative w.r.t. z for v*w at location of v-component. For details on implementation, see duvdx.
  inline RealType dvwdz(VectorField& velocity, const Parameters& parameters, int i, int j, int k) {
    const RealType mp0002  = parameters.meshsize->getDz(i, j, k);
    const RealType hzShort = 0.5 * mp0002;                                             // Distance of corner points in x-direction from center v-value
    const RealType hzLong0 = 0.5 * (mp0002 + parameters.meshsize->getDz(i, j, k - 1)); // Distance between center and west
                                                                                       // v-value
    const RealType hzLong1 = 0.5 * (mp0002 + parameters.meshsize->getDz(i, j, k + 1)); // Distance between center and east
    // v-value
    const RealType mp0001  = parameters.meshsize->getDy(i, j, k);
    const RealType hyShort = 0.5 * mp0001;                                             // Distance of center u-value from upper edge of cell
    const RealType hyLong  = 0.5 * (mp0001 + parameters.meshsize->getDy(i, j + 1, k)); // Distance of north and center u-value

    const RealType w00 = velocity.getVectorElement(i, j, k, 2);
    const RealType w10 = velocity.getVectorElement(i, j + 1, k, 2);
    const RealType v00 = velocity.getVectorElement(i, j, k, 1);
    const RealType v01 = velocity.getVectorElement(i, j, k + 1, 1);

    const RealType w0M1 = velocity.getVectorElement(i, j, k - 1, 2);
    const RealType w1M1 = velocity.getVectorElement(i, j + 1, k - 1, 2);
    const RealType v0M1 = velocity.getVectorElement(i, j, k - 1, 1);

    const RealType secondOrder = (((hyLong - hyShort) / hyLong * w00 + hyShort / hyLong * w10) * ((hzLong1 - hzShort) / hzLong1 * v00 + hzShort / hzLong1 * v01)
                                  - ((hyLong - hyShort) / hyLong * w0M1 + hyShort / hyLong * w1M1) * ((hzLong0 - hzShort) / hzLong0 * v00 + hzShort / hzLong0 * v0M1))
                                 / (2.0 * hzShort);

    const RealType kr = (hyLong - hyShort) / hyLong * w00 + hyShort / hyLong * w10;
    const RealType kl = (hyLong - hyShort) / hyLong * w0M1 + hyShort / hyLong * w1M1;

    const RealType firstOrder = 1.0 / (4.0 * hzShort) * (kr * (v00 + v01) - kl * (v0M1 + v00) + fabs(kr) * (v00 - v01) - fabs(kl) * (v0M1 - v00));

    const RealType tmp2 = (1.0 - parameters.solver.gamma) * secondOrder + parameters.solver.gamma * firstOrder;

    return tmp2;
  }

  // First derivative of u*u w.r.t. x, evaluated at location of u-component.
  inline RealType du2dx(VectorField& velocity, const Parameters& parameters, int i, int j, int k = 0) {
    const RealType mp0000  = parameters.meshsize->getDx(i, j, k);
    const RealType dxShort = 0.5 * mp0000;
    const RealType dxLong1 = 0.5 * (mp0000 + parameters.meshsize->getDx(i + 1, j, k));

    const RealType u0  = velocity.getVectorElement(i, j, k, 0);
    const RealType uM1 = velocity.getVectorElement(i - 1, j, k, 0);
    const RealType u1  = velocity.getVectorElement(i + 1, j, k, 0);

    const RealType kr = (u0 + u1) / 2;
    const RealType kl = (u0 + uM1) / 2;

    // Central difference expression which is second-order accurate for uniform meshes. We interpolate u half-way
    // between neighboured u-component values and afterwards build the central difference for u*u.

    const RealType secondOrder = ((u0 + u1) * (u0 + u1) - (u0 + uM1) * (u0 + uM1)) / (4 * dxLong1);

    // Donor-cell like derivative expression. We evaluate u half-way between neighboured u-components and use this as a
    // prediction of the transport direction.
    const RealType firstOrder = 1.0 / (4.0 * dxShort) * (kr * (u0 + u1) - kl * (uM1 + u0) + fabs(kr) * (u0 - u1) - fabs(kl) * (uM1 - u0));

    // Return linear combination of central- and upwind difference
    const RealType tmp2 = (1.0 - parameters.solver.gamma) * secondOrder + parameters.solver.gamma * firstOrder;

    return tmp2;
  }

  // First derivative of v*v w.r.t. y, evaluated at location of v-component. For details, see du2dx.
  inline RealType dv2dy(VectorField& velocity, const Parameters& parameters, int i, int j, int k = 0) {
    const RealType mp0001  = parameters.meshsize->getDy(i, j, k);
    const RealType dyShort = 0.5 * mp0001;
    const RealType dyLong1 = 0.5 * (mp0001 + parameters.meshsize->getDy(i, j + 1, k));

    const RealType v0  = velocity.getVectorElement(i, j, k, 1);
    const RealType vM1 = velocity.getVectorElement(i, j - 1, k, 1);
    const RealType v1  = velocity.getVectorElement(i, j + 1, k, 1);

    const RealType kr = (v0 + v1) / 2;
    const RealType kl = (v0 + vM1) / 2;

    const RealType secondOrder = ((v0 + v1) * (v0 + v1) - (v0 + vM1) * (v0 + vM1)) / (4 * dyLong1);

    const RealType firstOrder = 1.0 / (4.0 * dyShort) * (kr * (v0 + v1) - kl * (vM1 + v0) + fabs(kr) * (v0 - v1) - fabs(kl) * (vM1 - v0));

    const RealType tmp2 = (1.0 - parameters.solver.gamma) * secondOrder + parameters.solver.gamma * firstOrder;

    return tmp2;
  }

  // First derivative of w*w w.r.t. z, evaluated at location of w-component. For details, see du2dx.
  inline RealType dw2dz(VectorField& velocity, const Parameters& parameters, int i, int j, int k) {
    const RealType mp0002  = parameters.meshsize->getDz(i, j, k);
    const RealType dzShort = 0.5 * mp0002;
    const RealType dzLong1 = 0.5 * (mp0002 + parameters.meshsize->getDz(i, j, k + 1));

    const RealType w0  = velocity.getVectorElement(i, j, k, 2);
    const RealType wM1 = velocity.getVectorElement(i, j, k - 1, 2);
    const RealType w1  = velocity.getVectorElement(i, j, k + 1, 2);

    const RealType kr = (w0 + w1) / 2;
    const RealType kl = (w0 + wM1) / 2;

    const RealType secondOrder = ((w0 + w1) * (w0 + w1) - (w0 + wM1) * (w0 + wM1)) / (4 * dzLong1);

    const RealType firstOrder = 1.0 / (4.0 * dzShort) * (kr * (w0 + w1) - kl * (wM1 + w0) + fabs(kr) * (w0 - w1) - fabs(kl) * (wM1 - w0));

    const RealType tmp2 = (1.0 - parameters.solver.gamma) * secondOrder + parameters.solver.gamma * firstOrder;

    return tmp2;
  }

  inline RealType computeF2D(VectorField& velocity, const Parameters& parameters, RealType dt, int i, int j) {
    return velocity.getVectorElement(i, j, 0)
           + dt
               * (-du2dx(velocity, parameters, i, j) - duvdy(velocity, parameters, i, j) + 1 / parameters.flow.Re * (d2udx2(velocity, parameters, i, j) + d2udy2(velocity, parameters, i, j)));
  }

  inline RealType computeTurbulentF2D(VectorField& velocity, ScalarField& viscosity, const Parameters& parameters, RealType dt, int i, int j) {
    const RealType re = parameters.flow.Re;
    return velocity.getVectorElement(i, j, 0)
           + dt
               * (-du2dx(velocity, parameters, i, j) - duvdy(velocity, parameters, i, j) + 2 * dVdxdudx(velocity, viscosity, parameters, re, i, j) + dVdydudydvdx(velocity, viscosity, parameters, re, i, j));
  }

  inline RealType computeF3D(VectorField& velocity, const Parameters& parameters, RealType dt, int i, int j, int k) {
    return velocity.getVectorElement(i, j, k, 0)
           + dt
               * (-du2dx(velocity, parameters, i, j, k) - duvdy(velocity, parameters, i, j, k) - duwdz(velocity, parameters, i, j, k) + 1 / parameters.flow.Re * (d2udx2(velocity, parameters, i, j, k) + d2udy2(velocity, parameters, i, j, k) + d2udz2(velocity, parameters, i, j, k)));
  }

  inline RealType computeTurbulentF3D(VectorField& velocity, ScalarField& viscosity, const Parameters& parameters, RealType dt, int i, int j, int k) {
    const RealType re = parameters.flow.Re;
    return velocity.getVectorElement(i, j, k, 0)
           + dt
               * (-du2dx(velocity, parameters, i, j, k) - duvdy(velocity, parameters, i, j, k) - duwdz(velocity, parameters, i, j, k) + 2 * dVdxdudx(velocity, viscosity, parameters, re, i, j, k) + dVdydudydvdx(velocity, viscosity, parameters, re, i, j, k) + dVdzdudzdwdx(velocity, viscosity, parameters, re, i, j, k));
  }

  inline RealType computeG2D(VectorField& velocity, const Parameters& parameters, RealType dt, int i, int j) {
    return velocity.getVectorElement(i, j, 1)
           + dt
               * (-duvdx(velocity, parameters, i, j) - dv2dy(velocity, parameters, i, j) + 1 / parameters.flow.Re * (d2vdx2(velocity, parameters, i, j) + d2vdy2(velocity, parameters, i, j)));
  }

  inline RealType computeTurbulentG2D(VectorField& velocity, ScalarField& viscosity, const Parameters& parameters, RealType dt, int i, int j) {
    const RealType re = parameters.flow.Re;
    return velocity.getVectorElement(i, j, 1)
           + dt
               * (-duvdx(velocity, parameters, i, j) - dv2dy(velocity, parameters, i, j) + dVdxdvdxdudy(velocity, viscosity, parameters, re, i, j) + 2 * dVdydvdy(velocity, viscosity, parameters, re, i, j));
  }

  inline RealType computeG3D(VectorField& velocity, const Parameters& parameters, RealType dt, int i, int j, int k) {
    return velocity.getVectorElement(i, j, k, 1)
           + dt
               * (-dv2dy(velocity, parameters, i, j, k) - duvdx(velocity, parameters, i, j, k) - dvwdz(velocity, parameters, i, j, k) + 1 / parameters.flow.Re * (d2vdx2(velocity, parameters, i, j, k) + d2vdy2(velocity, parameters, i, j, k) + d2vdz2(velocity, parameters, i, j, k)));
  }

  inline RealType computeTurbulentG3D(VectorField& velocity, ScalarField& viscosity, const Parameters& parameters, RealType dt, int i, int j, int k) {
    const RealType re = parameters.flow.Re;
    return velocity.getVectorElement(i, j, k, 1)
           + dt
               * (-dv2dy(velocity, parameters, i, j, k) - duvdx(velocity, parameters, i, j, k) - dvwdz(velocity, parameters, i, j, k) + dVdxdvdxdudy(velocity, viscosity, parameters, re, i, j, k) + 2 * dVdydvdy(velocity, viscosity, parameters, re, i, j, k) + dVdzdvdzdwdy(velocity, viscosity, parameters, re, i, j, k));
  }

  inline RealType computeH3D(VectorField& velocity, const Parameters& parameters, RealType dt, int i, int j, int k) {
    return velocity.getVectorElement(i, j, k, 2)
           + dt
               * (-dw2dz(velocity, parameters, i, j, k) - duwdx(velocity, parameters, i, j, k) - dvwdy(velocity, parameters, i, j, k) + 1 / parameters.flow.Re * (d2wdx2(velocity, parameters, i, j, k) + d2wdy2(velocity, parameters, i, j, k) + d2wdz2(velocity, parameters, i, j, k)));
  }

  inline RealType computeTurbulentH3D(VectorField& velocity, ScalarField& viscosity, const Parameters& parameters, RealType dt, int i, int j, int k) {
    const RealType re = parameters.flow.Re;
    return velocity.getVectorElement(i, j, k, 2)
           + dt
               * (-dw2dz(velocity, parameters, i, j, k) - duwdx(velocity, parameters, i, j, k) - dvwdy(velocity, parameters, i, j, k) + (dVdxdwdxdudz(velocity, viscosity, parameters, re, i, j, k) + dVdydwdydvdz(velocity, viscosity, parameters, re, i, j, k) + 2 * dVdzdwdz(velocity, viscosity, parameters, re, i, j, k)));
  }

} // namespace Stencils
