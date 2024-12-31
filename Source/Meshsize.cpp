#include "StdAfx.hpp"

#include "Meshsize.hpp"

#include "Parameters.hpp"

UniformMeshsize::UniformMeshsize(const GeometricParameters& geometricParameters, const ParallelParameters& parallelParameters):
  dx_(geometricParameters.lengthX / geometricParameters.sizeX),
  dy_(geometricParameters.lengthY / geometricParameters.sizeY),
  dz_(geometricParameters.dim == 3 ? geometricParameters.lengthZ / geometricParameters.sizeZ : 0.0),
  firstCornerX_(parallelParameters.firstCorner[0]),
  firstCornerY_(parallelParameters.firstCorner[1]),
  firstCornerZ_(geometricParameters.dim == 3 ? parallelParameters.firstCorner[2] : 0) {

  if (dx_ <= 0.0) {
    throw std::runtime_error("dx <= 0.0!");
  }
  if (dy_ <= 0.0) {
    throw std::runtime_error("dy <= 0.0!");
  }
  if (geometricParameters.dim == 3) {
    if (dz_ <= 0.0) {
      throw std::runtime_error("dz <= 0.0!");
    }
  }
}

TanhMeshStretching::TanhMeshStretching(const GeometricParameters& geometricParameters, const ParallelParameters& parallelParameters):
  uniformMeshsize_(geometricParameters, parallelParameters),
  lengthX_(geometricParameters.lengthX),
  lengthY_(geometricParameters.lengthY),
  lengthZ_(geometricParameters.dim == 3 ? geometricParameters.lengthZ : 0.0),
  sizeX_(geometricParameters.sizeX),
  sizeY_(geometricParameters.sizeY),
  sizeZ_(geometricParameters.dim == 3 ? geometricParameters.sizeZ : 1),
  firstCornerX_(parallelParameters.firstCorner[0]),
  firstCornerY_(parallelParameters.firstCorner[1]),
  firstCornerZ_(geometricParameters.dim == 3 ? parallelParameters.firstCorner[2] : 0),
  stretchX_(geometricParameters.meshsizeType == TanhStretching ? static_cast<bool>(geometricParameters.stretchX) : false),
  stretchY_(geometricParameters.meshsizeType == TanhStretching ? static_cast<bool>(geometricParameters.stretchY) : false),
  stretchZ_(geometricParameters.meshsizeType == TanhStretching ? static_cast<bool>(geometricParameters.stretchZ) : false),
  deltaS_(2.7),
  tanhDeltaS_(tanh(2.7)) // This parameters is chosen as 2.7 as used also in the dissertation by Tobias Neckel
  ,
  dxMin_(
    (geometricParameters.meshsizeType == TanhStretching ? static_cast<bool>(geometricParameters.stretchX) : false)
      ? 0.5 * geometricParameters.lengthX * (1.0 + tanh(deltaS_ * (2.0 / sizeX_ - 1.0)) / tanhDeltaS_)
      : uniformMeshsize_.getDx(0, 0)
  ),
  dyMin_(
    (geometricParameters.meshsizeType == TanhStretching ? static_cast<bool>(geometricParameters.stretchX) : false)
      ? 0.5 * geometricParameters.lengthY * (1.0 + tanh(deltaS_ * (2.0 / sizeY_ - 1.0)) / tanhDeltaS_)
      : uniformMeshsize_.getDy(0, 0)
  ),
  dzMin_(
    (geometricParameters.meshsizeType == TanhStretching ? static_cast<bool>(geometricParameters.stretchX) : false)
      ? 0.5 * geometricParameters.lengthZ * (1.0 + tanh(deltaS_ * (2.0 / sizeZ_ - 1.0)) / tanhDeltaS_)
      : uniformMeshsize_.getDz(0, 0, 0)
  ) {}
