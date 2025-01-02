#include "StdAfx.hpp"

#include "FieldStencilDelegate.hpp"

FieldStencilDelegate::FieldStencilDelegate(const StencilType type):
  stencilType(type),
  //distanceStencil(),
  fghStencil(),
  // maxUStencil(),
  // maxViscStencil(),
  // //obstacleCoordinatesStencil(),
  obstacleStencil(),
  turbulentFGHStencil(),
  velocityStencil(),
  viscosityStencil(),
  rhsStencil() {};