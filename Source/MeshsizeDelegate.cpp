#include "StdAfx.hpp"

#include "MeshsizeDelegate.hpp"
#include "Parameters.hpp"

MeshsizeDelegate::MeshsizeDelegate(const GeometricParameters& geometricParameters, const ParallelParameters& parallelParameters):
  uniform(geometricParameters, parallelParameters),
  tanh(geometricParameters, parallelParameters),
  type(static_cast<MeshsizeType>(geometricParameters.meshsizeType)) {};