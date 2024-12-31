#pragma once

#include "Definitions.hpp"

// Forward declaration of Parameters
class GeometricParameters;
class ParallelParameters;

enum MeshsizeType { Uniform = 0, TanhStretching = 1 };

/** Implements a uniform, equidistant grid spacing */
class UniformMeshsize {
private:
  const RealType dx_;
  const RealType dy_;
  const RealType dz_;
  const int      firstCornerX_;
  const int      firstCornerY_;
  const int      firstCornerZ_;

public:
  UniformMeshsize(const GeometricParameters& geometricParameters, const ParallelParameters& parallelParameters);
  ~UniformMeshsize() = default;

  inline RealType getDx([[maybe_unused]] int i, [[maybe_unused]] int j) const { return dx_; }
  inline RealType getDy([[maybe_unused]] int i, [[maybe_unused]] int j) const { return dy_; }

  inline RealType getDx([[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k) const { return dx_; }
  inline RealType getDy([[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k) const { return dy_; }
  inline RealType getDz([[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k) const { return dz_; }

  inline RealType getPosX([[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k) const { return dx_ * (firstCornerX_ - 2 + i); }
  inline RealType getPosY([[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k) const { return dy_ * (firstCornerY_ - 2 + j); }
  inline RealType getPosZ([[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k) const { return dz_ * (firstCornerZ_ - 2 + k); }

  inline RealType getPosX(int i, int j) const { return getPosX(i, j, 0); }
  inline RealType getPosY(int i, int j) const { return getPosY(i, j, 0); }

  inline RealType getDxMin() const { return dx_; }
  inline RealType getDyMin() const { return dy_; }
  inline RealType getDzMin() const { return dz_; }
};

/**
 * Implements a stretched mesh for e.g. channel flow. For each dimension, a stretching of the mesh can be introduced
 * towards the outer boundaries, i.e. if stretchX is true (in constructor), then the mesh will be finer close to the
 * left and right boundary. The stretching is based on a formular involving tanh-functions, as e.g. used in the
 * dissertation by Tobias Neckel, Chair of Scientific Computing in Computer Science (TUM SCCS). For non-stretched
 * meshes, the UniformMeshsize implementation is used to create a uniform mesh.
 */
class TanhMeshStretching {
private:
  const UniformMeshsize uniformMeshsize_;
  const RealType        lengthX_;
  const RealType        lengthY_;
  const RealType        lengthZ_;
  const int             sizeX_;
  const int             sizeY_;
  const int             sizeZ_;
  const int             firstCornerX_;
  const int             firstCornerY_;
  const int             firstCornerZ_;
  const bool            stretchX_;
  const bool            stretchY_;
  const bool            stretchZ_;
  const RealType        deltaS_;
  const RealType        tanhDeltaS_;
  const RealType        dxMin_;
  const RealType        dyMin_;
  const RealType        dzMin_;

  // Computes the coordinate of the lower/left/front corner of the 1D-cell at index i w.r.t. having "size" cells along
  // an interval of length "length". We refer to local indexing, so "firstCorner" denotes the first non-ghost cell index
  // of this process. We use a stretched mesh for all nodes inside the comput. bounding box, and a regular mesh outside
  // this box, using the meshsize of the next inner cell.
  inline RealType computeCoordinate(int i, int firstCorner, int size, RealType length, RealType dxMin) const {
    const int index = i - 2 + firstCorner;
    if (index < 0) {
      // Equidistant mesh on lower/left part
      return dxMin * index;
    } else if (index > size - 1) {
      // Equidistant mesh on upper/right part
      return length + dxMin * (index - size);
    } else {
      // Stretched mesh on lower half of channel -> we check if we are in lower 50% and then use stretching for 2.0 * p
      RealType p = (static_cast<RealType>(index)) / size;
      if (p < 0.5) {
        return 0.5 * length * (1.0 + tanh(deltaS_ * (2.0 * p - 1.0)) / tanhDeltaS_);
      } else {
        // Stretched mesh on upper half of channel -> we mirror the stretching
        p = (static_cast<RealType>(size) - index) / size;
        return length - 0.5 * length * (1.0 + tanh(deltaS_ * (2.0 * p - 1.0)) / tanhDeltaS_);
      }
    }
  }

  // Returns the meshsize based on vertex coordinates that span the respective 1D-cell
  inline RealType getMeshsize(int i, int firstCorner, int size, RealType length, RealType dxMin) const {
    const RealType pos0 = computeCoordinate(i, firstCorner, size, length, dxMin);
    const RealType pos1 = computeCoordinate(i + 1, firstCorner, size, length, dxMin);
#ifndef NDEBUG
    if (pos1 - pos0 < 1.0e-12) {
      throw std::runtime_error("Error TanhMeshStretching::getMeshsize(): dx < 1.0e-12!");
    }
#endif
    return pos1 - pos0;
  }

public:
  TanhMeshStretching(const GeometricParameters& geometricParameters, const ParallelParameters& parallelParameters);
  ~TanhMeshStretching() = default;

  inline RealType getDx(int i, int j) const {
    if (stretchX_) {
      return getMeshsize(i, firstCornerX_, sizeX_, lengthX_, dxMin_);
    } else {
      return uniformMeshsize_.getDx(i, j);
    }
  }

  inline RealType getDy(int i, int j) const {
    if (stretchY_) {
      return getMeshsize(j, firstCornerY_, sizeY_, lengthY_, dyMin_);
    } else {
      return uniformMeshsize_.getDy(i, j);
    }
  }

  inline RealType getDx(int i, int j, [[maybe_unused]] int k) const { return getDx(i, j); }

  inline RealType getDy(int i, int j, [[maybe_unused]] int k) const { return getDy(i, j); }

  inline RealType getDz(int i, int j, int k) const {
    if (stretchZ_) {
      return getMeshsize(k, firstCornerZ_, sizeZ_, lengthZ_, dzMin_);
    } else {
      return uniformMeshsize_.getDz(i, j, k);
    }
  }

  inline RealType getPosX(int i, int j, int k) const {
    if (stretchX_) {
      return computeCoordinate(i, firstCornerX_, sizeX_, lengthX_, dxMin_);
    } else {
      return uniformMeshsize_.getPosX(i, j, k);
    }
  }

  inline RealType getPosY(int i, int j, int k) const {
    if (stretchY_) {
      return computeCoordinate(j, firstCornerY_, sizeY_, lengthY_, dyMin_);
    } else {
      return uniformMeshsize_.getPosY(i, j, k);
    }
  }

  inline RealType getPosZ(int i, int j, int k) const {
    if (stretchZ_) {
      return computeCoordinate(k, firstCornerZ_, sizeZ_, lengthZ_, dzMin_);
    } else {
      return uniformMeshsize_.getPosZ(i, j, k);
    }
  }

  inline RealType getPosX(int i, int j) const { return getPosX(i, j, 0); }
  inline RealType getPosY(int i, int j) const { return getPosY(i, j, 0); }

  inline RealType getDxMin() const { return dxMin_; }
  inline RealType getDyMin() const { return dyMin_; }
  inline RealType getDzMin() const { return dzMin_; }
};
