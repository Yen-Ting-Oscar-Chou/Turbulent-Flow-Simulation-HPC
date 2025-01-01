#pragma once

#include "Meshsize.hpp"

class MeshsizeDelegate: public Meshsize {
private:
  UniformMeshsize    uniform;
  TanhMeshStretching tanh;

public:
  const MeshsizeType type;
  MeshsizeDelegate(const GeometricParameters& geometricParameters, const ParallelParameters& parallelParameters);

  ~MeshsizeDelegate() = default;

  inline virtual RealType getDx(int i, int j) const override { return (type == Uniform) * uniform.getDx(i, j) + (type == TanhStretching) * tanh.getDx(i, j); }

  inline virtual RealType getDy(int i, int j) const override { return (type == Uniform) * uniform.getDy(i, j) + (type == TanhStretching) * tanh.getDy(i, j); }

  inline virtual RealType getDx(int i, int j, [[maybe_unused]] int k) const override { return getDx(i, j); }

  inline virtual RealType getDy(int i, int j, [[maybe_unused]] int k) const override { return getDy(i, j); }

  inline virtual RealType getDz(int i, int j, int k) const override { return (type == Uniform) * uniform.getDz(i, j, k) + (type == TanhStretching) * tanh.getDz(i, j, k); }

  inline virtual RealType getPosX(int i, int j, int k) const override { return (type == Uniform) * uniform.getPosX(i, j, k) + (type == TanhStretching) * tanh.getPosX(i, j, k); }

  inline virtual RealType getPosY(int i, int j, int k) const override { return (type == Uniform) * uniform.getPosY(i, j, k) + (type == TanhStretching) * tanh.getPosY(i, j, k); }

  inline virtual RealType getPosZ(int i, int j, int k) const override { return (type == Uniform) * uniform.getPosZ(i, j, k) + (type == TanhStretching) * tanh.getPosZ(i, j, k); }

  inline virtual RealType getPosX(int i, int j) const override { return getPosX(i, j, 0); }
  inline virtual RealType getPosY(int i, int j) const override { return getPosY(i, j, 0); }

  inline virtual RealType getDxMin() const override { return (type == Uniform) * uniform.getDxMin() + (type == TanhStretching) * tanh.getDxMin(); }
  inline virtual RealType getDyMin() const override { return (type == Uniform) * uniform.getDyMin() + (type == TanhStretching) * tanh.getDyMin(); }
  inline virtual RealType getDzMin() const override { return (type == Uniform) * uniform.getDzMin() + (type == TanhStretching) * tanh.getDzMin(); }
};
