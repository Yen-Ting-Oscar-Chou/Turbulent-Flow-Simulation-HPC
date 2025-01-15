#include "StdAfx.hpp"

#include "StencilDelegate.hpp"

StencilDelegate::StencilDelegate(StencilType type):
  stencilType(type),
  fghStencil(),
  turbulentFGHStencil(),
  velocityStencil(),
  viscosityStencil(),
  rhsStencil(),
  maxUStencil(),
  obstacleStencil(),
  neumannFGHBoundaryStencil(),
  neumannVelocityBoundaryStencil(),
  movingWallFGHStencil(),
  movingWallVelocityStencil(),
  bfInputFGHStencil(),
  bfInputVelocityStencil() {};