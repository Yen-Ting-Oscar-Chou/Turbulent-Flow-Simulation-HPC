#pragma once

#include "Meshsize.hpp"

class MeshsizeDelegate {

public:
  UniformMeshsize    uniform;
  TanhMeshStretching tanh;
  MeshsizeType       type;

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

  static MeshsizeDelegate* mapToGPU(int hostDevice, int targetDevice, MeshsizeDelegate& meshsizeDelegate) {
    size_t            meshsizeDelegateSize = sizeof(meshsizeDelegate);
    MeshsizeDelegate* meshsizeDelegateGPU  = static_cast<MeshsizeDelegate*>(omp_target_alloc(meshsizeDelegateSize, targetDevice));
    if (!meshsizeDelegateGPU) {
      std::cerr << "Error: Allocation failed for meshsizeDelegateGPU pointer." << std::endl;
    }

    bool associatedMeshsize = omp_target_associate_ptr(&meshsizeDelegate, meshsizeDelegateGPU, meshsizeDelegateSize, 0, targetDevice) == 0;
    if (!associatedMeshsize) {
      std::cout << "Error: MeshsizeDelegate could not be associated to GPU pointer." << std::endl;
    }

    omp_target_memcpy(&(meshsizeDelegateGPU->type), &(meshsizeDelegate.type), sizeof(MeshsizeType), 0, 0, targetDevice, hostDevice);

    // meshsizeDelegate.uniform
    {
      omp_target_memcpy(&(meshsizeDelegateGPU->uniform.dx_), &(meshsizeDelegate.uniform.dx_), sizeof(RealType), 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&(meshsizeDelegateGPU->uniform.dy_), &(meshsizeDelegate.uniform.dy_), sizeof(RealType), 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&(meshsizeDelegateGPU->uniform.dz_), &(meshsizeDelegate.uniform.dz_), sizeof(RealType), 0, 0, targetDevice, hostDevice);

      omp_target_memcpy(&(meshsizeDelegateGPU->uniform.firstCornerX_), &(meshsizeDelegate.uniform.firstCornerX_), sizeof(int), 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&(meshsizeDelegateGPU->uniform.firstCornerY_), &(meshsizeDelegate.uniform.firstCornerY_), sizeof(int), 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&(meshsizeDelegateGPU->uniform.firstCornerZ_), &(meshsizeDelegate.uniform.firstCornerZ_), sizeof(int), 0, 0, targetDevice, hostDevice);
    }

    // meshsizeDelegate.tanh
    {
      // meshsizeDelegate.tanh.uniformMeshsize_
      {
        omp_target_memcpy(&(meshsizeDelegateGPU->tanh.uniformMeshsize_.dx_), &(meshsizeDelegate.tanh.uniformMeshsize_.dx_), sizeof(RealType), 0, 0, targetDevice, hostDevice);
        omp_target_memcpy(&(meshsizeDelegateGPU->tanh.uniformMeshsize_.dy_), &(meshsizeDelegate.tanh.uniformMeshsize_.dy_), sizeof(RealType), 0, 0, targetDevice, hostDevice);
        omp_target_memcpy(&(meshsizeDelegateGPU->tanh.uniformMeshsize_.dz_), &(meshsizeDelegate.tanh.uniformMeshsize_.dz_), sizeof(RealType), 0, 0, targetDevice, hostDevice);

        omp_target_memcpy(
          &(meshsizeDelegateGPU->tanh.uniformMeshsize_.firstCornerX_), &(meshsizeDelegate.tanh.uniformMeshsize_.firstCornerX_), sizeof(int), 0, 0, targetDevice, hostDevice
        );
        omp_target_memcpy(
          &(meshsizeDelegateGPU->tanh.uniformMeshsize_.firstCornerY_), &(meshsizeDelegate.tanh.uniformMeshsize_.firstCornerY_), sizeof(int), 0, 0, targetDevice, hostDevice
        );
        omp_target_memcpy(
          &(meshsizeDelegateGPU->tanh.uniformMeshsize_.firstCornerZ_), &(meshsizeDelegate.tanh.uniformMeshsize_.firstCornerZ_), sizeof(int), 0, 0, targetDevice, hostDevice
        );
      }

      omp_target_memcpy(&(meshsizeDelegateGPU->tanh.lengthX_), &(meshsizeDelegate.tanh.lengthX_), sizeof(RealType), 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&(meshsizeDelegateGPU->tanh.lengthY_), &(meshsizeDelegate.tanh.lengthY_), sizeof(RealType), 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&(meshsizeDelegateGPU->tanh.lengthZ_), &(meshsizeDelegate.tanh.lengthZ_), sizeof(RealType), 0, 0, targetDevice, hostDevice);

      omp_target_memcpy(&(meshsizeDelegateGPU->tanh.sizeX_), &(meshsizeDelegate.tanh.sizeX_), sizeof(int), 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&(meshsizeDelegateGPU->tanh.sizeY_), &(meshsizeDelegate.tanh.sizeY_), sizeof(int), 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&(meshsizeDelegateGPU->tanh.sizeZ_), &(meshsizeDelegate.tanh.sizeZ_), sizeof(int), 0, 0, targetDevice, hostDevice);

      omp_target_memcpy(&(meshsizeDelegateGPU->tanh.firstCornerX_), &(meshsizeDelegate.tanh.firstCornerX_), sizeof(int), 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&(meshsizeDelegateGPU->tanh.firstCornerY_), &(meshsizeDelegate.tanh.firstCornerY_), sizeof(int), 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&(meshsizeDelegateGPU->tanh.firstCornerZ_), &(meshsizeDelegate.tanh.firstCornerZ_), sizeof(int), 0, 0, targetDevice, hostDevice);

      omp_target_memcpy(&(meshsizeDelegateGPU->tanh.stretchX_), &(meshsizeDelegate.tanh.stretchX_), sizeof(bool), 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&(meshsizeDelegateGPU->tanh.stretchY_), &(meshsizeDelegate.tanh.stretchY_), sizeof(bool), 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&(meshsizeDelegateGPU->tanh.stretchZ_), &(meshsizeDelegate.tanh.stretchZ_), sizeof(bool), 0, 0, targetDevice, hostDevice);

      omp_target_memcpy(&(meshsizeDelegateGPU->tanh.deltaS_), &(meshsizeDelegate.tanh.deltaS_), sizeof(RealType), 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&(meshsizeDelegateGPU->tanh.tanhDeltaS_), &(meshsizeDelegate.tanh.tanhDeltaS_), sizeof(RealType), 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&(meshsizeDelegateGPU->tanh.dxMin_), &(meshsizeDelegate.tanh.dxMin_), sizeof(RealType), 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&(meshsizeDelegateGPU->tanh.dyMin_), &(meshsizeDelegate.tanh.dyMin_), sizeof(RealType), 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&(meshsizeDelegateGPU->tanh.dzMin_), &(meshsizeDelegate.tanh.dzMin_), sizeof(RealType), 0, 0, targetDevice, hostDevice);
    }
    return meshsizeDelegateGPU;
  }

  static void freeGPU(int hostDevice, int targetDevice, MeshsizeDelegate* meshsizeDelegateGPU, MeshsizeDelegate& meshsizeDelegate) {
    bool disassociatedMeshsizeDelegate = omp_target_disassociate_ptr(&meshsizeDelegate, targetDevice) == 0;
    if (!disassociatedMeshsizeDelegate) {
      std::cout << "Error: Could not disassociate MeshsizeDelegate pointer." << std::endl;
    }

    omp_target_free(meshsizeDelegateGPU, targetDevice);
  }
};
