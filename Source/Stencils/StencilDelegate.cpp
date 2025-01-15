#include "StdAfx.hpp"

#include "StencilDelegate.hpp"

StencilDelegate::StencilDelegate():
  fghStencil_(),
  turbulentFGHStencil_(),
  velocityStencil_(),
  viscosityStencil_(),
  rhsStencil_(),
  maxUStencil_(),
  obstacleStencil_(),
  neumannFGHBoundaryStencil_(),
  neumannVelocityBoundaryStencil_(),
  movingWallFGHStencil_(),
  movingWallVelocityStencil_(),
  bfInputFGHStencil_(),
  bfInputVelocityStencil_(),
  maxViscStencil_() {};