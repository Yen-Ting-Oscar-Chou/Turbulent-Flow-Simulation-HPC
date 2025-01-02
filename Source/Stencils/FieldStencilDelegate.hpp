#pragma once

#include "DistanceStencil.hpp"
#include "FGHStencil.hpp"
#include "MaxUStencil.hpp"
#include "MaxViscStencil.hpp"
#include "ObstacleCoordinatesStencil.hpp"
#include "ObstacleStencil.hpp"
#include "RHSStencil.hpp"
#include "TurbulentFGHStencil.hpp"
#include "VelocityStencil.hpp"
#include "ViscosityStencil.hpp"

enum StencilType { DISTANCE, FGH, MAXVELO, MAXVISC, OBSTACLE, OBSTCOORD, RHS, TURBFGH, VELOCITY, VISCOSITY };

class FieldStencilDelegate {
public:
  StencilType stencilType;
  // Stencils::DistanceStencil            distanceStencil;
  Stencils::FGHStencil     fghStencil;
  // Stencils::MaxUStencil    maxUStencil;
  // Stencils::MaxViscStencil maxViscStencil;
  // Stencils::ObstacleCoordinatesStencil obstacleCoordinatesStencil;
  Stencils::ObstacleStencil     obstacleStencil;
  Stencils::RHSStencil          rhsStencil;
  Stencils::TurbulentFGHStencil turbulentFGHStencil;
  Stencils::VelocityStencil     velocityStencil;
  Stencils::ViscosityStencil    viscosityStencil;
  
  FieldStencilDelegate(const StencilType type);
  ~FieldStencilDelegate() = default;

  inline void apply(const Parameters& parameters, FlowField& flowField, int i, int j) {
    if (stencilType == FGH) {
      fghStencil.apply(parameters, flowField, i, j);
    } else if(stencilType == MAXVELO) {
      // maxUStencil.apply(parameters, flowField, i, j);
      assert(false);
    } else if(stencilType == OBSTACLE) {
      obstacleStencil.apply(parameters, flowField, i, j);
    } else if(stencilType == RHS) {
      rhsStencil.apply(parameters, flowField, i, j);
    } else if(stencilType == VELOCITY) {
      velocityStencil.apply(parameters, flowField, i, j);
    } else {
      assert(false);
    }
  }

  inline void apply(const Parameters& parameters, FlowField& flowField, int i, int j, int k) {
    if (stencilType == FGH) {
      fghStencil.apply(parameters, flowField, i, j, k);
    } else if(stencilType == MAXVELO) {
      // maxUStencil.apply(parameters, flowField, i, j, k);
      assert(false);
    } else if(stencilType == OBSTACLE) {
      obstacleStencil.apply(parameters, flowField, i, j, k);
    } else if(stencilType == RHS) {
      rhsStencil.apply(parameters, flowField, i, j, k);
    } else if(stencilType == VELOCITY) {
      velocityStencil.apply(parameters, flowField, i, j, k);
    } else {
      assert(false);
    }
  }

  inline void apply(const Parameters& parameters, TurbulentFlowField& flowField, int i, int j) {
    if (stencilType == DISTANCE){
      //   distanceStencil.apply(parameters, flowField, i, j);
      assert(false);
    } else if(stencilType == MAXVISC) {
      // maxViscStencil.apply(parameters, flowField, i, j);
      assert(false);
    } else if (stencilType == TURBFGH) {
      turbulentFGHStencil.apply(parameters, flowField, i, j);
    } else if (stencilType == OBSTCOORD) {
      assert(false);
    }
  }

  inline void apply(const Parameters& parameters, TurbulentFlowField& flowField, int i, int j, int k) {
     if (stencilType == DISTANCE){
      //   distanceStencil.apply(parameters, flowField, i, j, k);
      assert(false);
    } else if(stencilType == MAXVISC) {
      // maxViscStencil.apply(parameters, flowField, i, j, k);
      assert(false);
    } else if (stencilType == TURBFGH) {
      turbulentFGHStencil.apply(parameters, flowField, i, j, k);
    } else if (stencilType == OBSTCOORD) {
      // TODO
      assert(false);
    }
  }

  inline StencilType getType() const { return stencilType; }
};
#pragma omp declare mapper(FieldStencilDelegate fsd) \
  map(fsd) \
  map(fsd.stencilType) \
  map(fsd.rhsStencil) \
  map(fsd.fghStencil) \
  map(fsd.obstacleStencil) \
  map(fsd.turbulentFGHStencil) \
  map(fsd.velocityStencil) \
  map(fsd.viscosityStencil)
  // map(fsd.maxUStencil) 
  // map(fsd.maxViscStencil) 