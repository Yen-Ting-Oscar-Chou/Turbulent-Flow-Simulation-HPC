#include "StdAfx.hpp"

#include "Parameters.hpp"

Parameters::Parameters(const GeometricParameters& geometricParameters, const ParallelParameters& parallelParameters):
  simulation{},
  timestep{},
  environment{},
  flow{},
  solver{},
  geometry(geometricParameters),
  walls{},
  vtk{},
  parallel(parallelParameters),
  stdOut{},
  bfStep{},
  turbulence{},
  meshsize(new MeshsizeDelegate(geometry, parallel)) {
}
