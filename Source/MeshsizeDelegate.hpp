#pragma once

#include "Meshsize.hpp"

class MeshsizeDelegate {
private:
  UniformMeshsize    uniform;
  TanhMeshStretching tanh;
  const MeshsizeType type;

public:
  MeshsizeDelegate(const GeometricParameters& geometricParameters, const ParallelParameters& parallelParameters);

  ~MeshsizeDelegate() = default;

  inline RealType getDx(int i, int j) const { return (type == Uniform) ? uniform.getDx(i, j) : tanh.getDx(i, j); }

  inline RealType getDy(int i, int j) const { return (type == Uniform) ? uniform.getDy(i, j) : tanh.getDy(i, j); }

  inline RealType getDx(int i, int j, [[maybe_unused]] int k) const { return getDx(i, j); }

  inline RealType getDy(int i, int j, [[maybe_unused]] int k) const { return getDy(i, j); }

  inline RealType getDz(int i, int j, int k) const { return (type == Uniform) ? uniform.getDz(i, j, k) : tanh.getDz(i, j, k); }

  inline RealType getPosX(int i, int j, int k) const { return (type == Uniform) ? uniform.getPosX(i, j, k) : tanh.getPosX(i, j, k); }

  inline RealType getPosY(int i, int j, int k) const { return (type == Uniform) ? uniform.getPosY(i, j, k) : tanh.getPosY(i, j, k); }

  inline RealType getPosZ(int i, int j, int k) const { return (type == Uniform) ? uniform.getPosZ(i, j, k) : tanh.getPosZ(i, j, k); }

  inline RealType getPosX(int i, int j) const { return getPosX(i, j, 0); }
  inline RealType getPosY(int i, int j) const { return getPosY(i, j, 0); }

  inline RealType getDxMin() const { return (type == Uniform) ? uniform.getDxMin() : tanh.getDxMin(); }
  inline RealType getDyMin() const { return (type == Uniform) ? uniform.getDyMin() : tanh.getDyMin(); }
  inline RealType getDzMin() const { return (type == Uniform) ? uniform.getDzMin() : tanh.getDzMin(); }

  inline MeshsizeType getType() const { return type; }
};
