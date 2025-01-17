#pragma once

#include "FlowField.hpp"

class TurbulentFlowField;

struct TurbulentFlowFieldPtrs {
  TurbulentFlowField* flowFieldPtr_;
  RealType*           pressureDataPtr_;
  RealType*           velocityDataPtr_;
  RealType*           rhsDataPtr_;
  RealType*           fghDataPtr_;
  int*                flagsDataPtr_;
  RealType*           viscosityDataPtr_;
  RealType*           distanceDataPtr_;

  TurbulentFlowFieldPtrs(
    TurbulentFlowField* flowFieldPtr,
    RealType*           pressureDataPtr,
    RealType*           velocityDataPtr,
    RealType*           rhsDataPtr,
    RealType*           fghDataPtr,
    int*                flagsDataPtr,
    RealType*           viscosityDataPtr,
    RealType*           distanceDataPtr
  ) {
    flowFieldPtr_     = flowFieldPtr;
    pressureDataPtr_  = pressureDataPtr;
    velocityDataPtr_  = velocityDataPtr;
    rhsDataPtr_       = rhsDataPtr;
    fghDataPtr_       = fghDataPtr;
    flagsDataPtr_     = flagsDataPtr;
    viscosityDataPtr_ = viscosityDataPtr;
    distanceDataPtr_  = distanceDataPtr;
  }
};

class TurbulentFlowField: public FlowField {
private:
  static TurbulentFlowField* mapTurbulentFlowFieldToGPU(int hostDevice, int targetDevice, TurbulentFlowField& flowField) {
    size_t flowFieldSize = sizeof(flowField);

    TurbulentFlowField* flowFieldGPU = static_cast<TurbulentFlowField*>(omp_target_alloc(flowFieldSize, targetDevice));
    if (!flowFieldGPU) {
      std::cout << "Error: Allocation failed for turbulentFlowFieldGPU pointer." << std::endl;
      return nullptr;
    }

    bool associatedFlowField = omp_target_associate_ptr(&flowField, flowFieldGPU, flowFieldSize, 0, targetDevice) == 0;
    if (!associatedFlowField) {
      std::cout << "Error: TurbulentFlowField could not be associated to GPU pointer." << std::endl;
      return nullptr;
    }

    omp_target_memcpy(&(flowFieldGPU->cellsX_), &(flowField.cellsX_), sizeof(int), 0, 0, targetDevice, hostDevice);
    omp_target_memcpy(&(flowFieldGPU->cellsY_), &(flowField.cellsY_), sizeof(int), 0, 0, targetDevice, hostDevice);
    omp_target_memcpy(&(flowFieldGPU->cellsZ_), &(flowField.cellsZ_), sizeof(int), 0, 0, targetDevice, hostDevice);
    omp_target_memcpy(&(flowFieldGPU->sizeX_), &(flowField.sizeX_), sizeof(int), 0, 0, targetDevice, hostDevice);
    omp_target_memcpy(&(flowFieldGPU->sizeY_), &(flowField.sizeY_), sizeof(int), 0, 0, targetDevice, hostDevice);
    omp_target_memcpy(&(flowFieldGPU->sizeZ_), &(flowField.sizeZ_), sizeof(int), 0, 0, targetDevice, hostDevice);

    return flowFieldGPU;
  }

  static RealType* mapRHSToGPU(int hostDevice, int targetDevice, TurbulentFlowField& flowField, TurbulentFlowField* flowFieldGPU) {
    size_t scalarFieldDataSize = flowField.RHS_.size_ * sizeof(RealType);

    RealType* rhsDataGPU = static_cast<RealType*>(omp_target_alloc(scalarFieldDataSize, targetDevice));
    if (!rhsDataGPU) {
      std::cout << "Error: Allocation failed for rhsDataGPU pointer." << std::endl;
      return nullptr;
    }

    bool copiedRHSData = omp_target_memcpy(rhsDataGPU, flowField.RHS_.data_, scalarFieldDataSize, 0, 0, targetDevice, hostDevice) == 0;
    if (!copiedRHSData) {
      std::cout << "Error: Copying RHS data to GPU not successful." << std::endl;
      return nullptr;
    }

    bool copiedRHSDataPtr = omp_target_memcpy(&flowFieldGPU->RHS_.data_, &rhsDataGPU, sizeof(RealType*), 0, 0, targetDevice, hostDevice) == 0;
    if (!copiedRHSDataPtr) {
      std::cout << "Error: Copying RHS data pointer to rhs object not successful." << std::endl;
      return nullptr;
    }

    omp_target_memcpy(&(flowFieldGPU->RHS_.sizeX_), &(flowField.RHS_.sizeX_), sizeof(int), 0, 0, targetDevice, hostDevice);
    omp_target_memcpy(&(flowFieldGPU->RHS_.sizeY_), &(flowField.RHS_.sizeY_), sizeof(int), 0, 0, targetDevice, hostDevice);
    omp_target_memcpy(&(flowFieldGPU->RHS_.sizeZ_), &(flowField.RHS_.sizeZ_), sizeof(int), 0, 0, targetDevice, hostDevice);
    omp_target_memcpy(&(flowFieldGPU->RHS_.components_), &(flowField.RHS_.components_), sizeof(int), 0, 0, targetDevice, hostDevice);
    omp_target_memcpy(&(flowFieldGPU->RHS_.size_), &(flowField.RHS_.size_), sizeof(int), 0, 0, targetDevice, hostDevice);

    return rhsDataGPU;
  }

  static RealType* mapPressureToGPU(int hostDevice, int targetDevice, TurbulentFlowField& flowField, TurbulentFlowField* flowFieldGPU) {
    size_t scalarFieldDataSize = flowField.pressure_.size_ * sizeof(RealType);

    RealType* pressureDataGPU = static_cast<RealType*>(omp_target_alloc(scalarFieldDataSize, targetDevice));
    if (!pressureDataGPU) {
      std::cout << "Error: Allocation failed for pressureDataGPU pointer." << std::endl;
      return nullptr;
    }

    bool copiedPressureData = omp_target_memcpy(pressureDataGPU, flowField.pressure_.data_, scalarFieldDataSize, 0, 0, targetDevice, hostDevice) == 0;
    if (!copiedPressureData) {
      std::cout << "Error: Copying pressure data to GPU not successful." << std::endl;
      return nullptr;
    }

    bool copiedPressureDataPtr = omp_target_memcpy(&flowFieldGPU->pressure_.data_, &pressureDataGPU, sizeof(RealType*), 0, 0, targetDevice, hostDevice) == 0;
    if (!copiedPressureDataPtr) {
      std::cout << "Error: Copying pressure data pointer to pressure object not successful." << std::endl;
      return nullptr;
    }

    omp_target_memcpy(&(flowFieldGPU->pressure_.sizeX_), &(flowField.pressure_.sizeX_), sizeof(int), 0, 0, targetDevice, hostDevice);
    omp_target_memcpy(&(flowFieldGPU->pressure_.sizeY_), &(flowField.pressure_.sizeY_), sizeof(int), 0, 0, targetDevice, hostDevice);
    omp_target_memcpy(&(flowFieldGPU->pressure_.sizeZ_), &(flowField.pressure_.sizeZ_), sizeof(int), 0, 0, targetDevice, hostDevice);
    omp_target_memcpy(&(flowFieldGPU->pressure_.components_), &(flowField.pressure_.components_), sizeof(int), 0, 0, targetDevice, hostDevice);
    omp_target_memcpy(&(flowFieldGPU->pressure_.size_), &(flowField.pressure_.size_), sizeof(int), 0, 0, targetDevice, hostDevice);

    return pressureDataGPU;
  }

  static int* mapFlagsToGPU(int hostDevice, int targetDevice, TurbulentFlowField& flowField, TurbulentFlowField* flowFieldGPU) {
    size_t intScalarFieldDataSize = flowField.flags_.size_ * sizeof(int);

    int* flagsDataGPU = static_cast<int*>(omp_target_alloc(intScalarFieldDataSize, targetDevice));
    if (!flagsDataGPU) {
      std::cout << "Error: Allocation failed for flagsDataGPU pointer." << std::endl;
      return nullptr;
    }

    bool copiedFlagsData = omp_target_memcpy(flagsDataGPU, flowField.flags_.data_, intScalarFieldDataSize, 0, 0, targetDevice, hostDevice) == 0;
    if (!copiedFlagsData) {
      std::cout << "Error: Copying flags data to GPU not successful." << std::endl;
      return nullptr;
    }

    bool copiedFlagsDataPtr = omp_target_memcpy(&flowFieldGPU->flags_.data_, &flagsDataGPU, sizeof(int*), 0, 0, targetDevice, hostDevice) == 0;
    if (!copiedFlagsDataPtr) {
      std::cout << "Error: Copying flags data pointer to flags object not successful." << std::endl;
      return nullptr;
    }

    omp_target_memcpy(&(flowFieldGPU->flags_.sizeX_), &(flowField.flags_.sizeX_), sizeof(int), 0, 0, targetDevice, hostDevice);
    omp_target_memcpy(&(flowFieldGPU->flags_.sizeY_), &(flowField.flags_.sizeY_), sizeof(int), 0, 0, targetDevice, hostDevice);
    omp_target_memcpy(&(flowFieldGPU->flags_.sizeZ_), &(flowField.flags_.sizeZ_), sizeof(int), 0, 0, targetDevice, hostDevice);
    omp_target_memcpy(&(flowFieldGPU->flags_.components_), &(flowField.flags_.components_), sizeof(int), 0, 0, targetDevice, hostDevice);
    omp_target_memcpy(&(flowFieldGPU->flags_.size_), &(flowField.flags_.size_), sizeof(int), 0, 0, targetDevice, hostDevice);

    return flagsDataGPU;
  }

  static RealType* mapVelocityToGPU(int hostDevice, int targetDevice, TurbulentFlowField& flowField, TurbulentFlowField* flowFieldGPU) {
    size_t vectorFieldDataSize = flowField.velocity_.size_ * sizeof(RealType);

    RealType* velocityDataGPU = static_cast<RealType*>(omp_target_alloc(vectorFieldDataSize, targetDevice));
    if (!velocityDataGPU) {
      std::cout << "Error: Allocation failed for velocityDataGPU pointer." << std::endl;
    }

    bool copiedVelocityData = omp_target_memcpy(velocityDataGPU, flowField.velocity_.data_, vectorFieldDataSize, 0, 0, targetDevice, hostDevice) == 0;
    if (!copiedVelocityData) {
      std::cout << "Error: Copying velocity data to GPU not successful." << std::endl;
      return nullptr;
    }

    bool copiedVelocityDataPtr = omp_target_memcpy(&flowFieldGPU->velocity_.data_, &velocityDataGPU, sizeof(RealType*), 0, 0, targetDevice, hostDevice) == 0;
    if (!copiedVelocityDataPtr) {
      std::cout << "Error: Copying velocity data pointer to velocity object not successful." << std::endl;
      return nullptr;
    }

    omp_target_memcpy(&(flowFieldGPU->velocity_.sizeX_), &(flowField.velocity_.sizeX_), sizeof(int), 0, 0, targetDevice, hostDevice);
    omp_target_memcpy(&(flowFieldGPU->velocity_.sizeY_), &(flowField.velocity_.sizeY_), sizeof(int), 0, 0, targetDevice, hostDevice);
    omp_target_memcpy(&(flowFieldGPU->velocity_.sizeZ_), &(flowField.velocity_.sizeZ_), sizeof(int), 0, 0, targetDevice, hostDevice);
    omp_target_memcpy(&(flowFieldGPU->velocity_.components_), &(flowField.velocity_.components_), sizeof(int), 0, 0, targetDevice, hostDevice);
    omp_target_memcpy(&(flowFieldGPU->velocity_.size_), &(flowField.velocity_.size_), sizeof(int), 0, 0, targetDevice, hostDevice);

    return velocityDataGPU;
  }

  static RealType* mapFGHToGPU(int hostDevice, int targetDevice, TurbulentFlowField& flowField, TurbulentFlowField* flowFieldGPU) {
    size_t vectorFieldDataSize = flowField.FGH_.size_ * sizeof(RealType);

    RealType* fghDataGPU = static_cast<RealType*>(omp_target_alloc(vectorFieldDataSize, targetDevice));
    if (!fghDataGPU) {
      std::cout << "Error: Allocation failed for fghDataGPU pointer." << std::endl;
      return nullptr;
    }

    bool copiedFGHData = omp_target_memcpy(fghDataGPU, flowField.FGH_.data_, vectorFieldDataSize, 0, 0, targetDevice, hostDevice) == 0;
    if (!copiedFGHData) {
      std::cout << "Error: Copying FGH data to GPU not successful." << std::endl;
      return nullptr;
    }

    bool copiedFGHDataPtr = omp_target_memcpy(&flowFieldGPU->FGH_.data_, &fghDataGPU, sizeof(RealType*), 0, 0, targetDevice, hostDevice) == 0;
    if (!copiedFGHDataPtr) {
      std::cout << "Error: Copying FGH data pointer to fgh object not successful." << std::endl;
      return nullptr;
    }

    omp_target_memcpy(&(flowFieldGPU->FGH_.sizeX_), &(flowField.FGH_.sizeX_), sizeof(int), 0, 0, targetDevice, hostDevice);
    omp_target_memcpy(&(flowFieldGPU->FGH_.sizeY_), &(flowField.FGH_.sizeY_), sizeof(int), 0, 0, targetDevice, hostDevice);
    omp_target_memcpy(&(flowFieldGPU->FGH_.sizeZ_), &(flowField.FGH_.sizeZ_), sizeof(int), 0, 0, targetDevice, hostDevice);
    omp_target_memcpy(&(flowFieldGPU->FGH_.components_), &(flowField.FGH_.components_), sizeof(int), 0, 0, targetDevice, hostDevice);
    omp_target_memcpy(&(flowFieldGPU->FGH_.size_), &(flowField.FGH_.size_), sizeof(int), 0, 0, targetDevice, hostDevice);

    return fghDataGPU;
  }

  static RealType* mapViscosityToGPU(int hostDevice, int targetDevice, TurbulentFlowField& flowField, TurbulentFlowField* flowFieldGPU) {
    size_t scalarFieldDataSize = flowField.viscosity_.size_ * sizeof(RealType);

    RealType* viscosityDataGPU = static_cast<RealType*>(omp_target_alloc(scalarFieldDataSize, targetDevice));
    if (!viscosityDataGPU) {
      std::cout << "Error: Allocation failed for viscosityDataGPU pointer." << std::endl;
      return nullptr;
    }

    bool copiedViscosityData = omp_target_memcpy(viscosityDataGPU, flowField.viscosity_.data_, scalarFieldDataSize, 0, 0, targetDevice, hostDevice) == 0;
    if (!copiedViscosityData) {
      std::cout << "Error: Copying viscosity data to GPU not successful." << std::endl;
      return nullptr;
    }

    bool copiedViscosityDataPtr = omp_target_memcpy(&flowFieldGPU->viscosity_.data_, &viscosityDataGPU, sizeof(RealType*), 0, 0, targetDevice, hostDevice) == 0;
    if (!copiedViscosityDataPtr) {
      std::cout << "Error: Copying viscosity data pointer to viscosity object not successful." << std::endl;
      return nullptr;
    }

    omp_target_memcpy(&(flowFieldGPU->viscosity_.sizeX_), &(flowField.viscosity_.sizeX_), sizeof(int), 0, 0, targetDevice, hostDevice);
    omp_target_memcpy(&(flowFieldGPU->viscosity_.sizeY_), &(flowField.viscosity_.sizeY_), sizeof(int), 0, 0, targetDevice, hostDevice);
    omp_target_memcpy(&(flowFieldGPU->viscosity_.sizeZ_), &(flowField.viscosity_.sizeZ_), sizeof(int), 0, 0, targetDevice, hostDevice);
    omp_target_memcpy(&(flowFieldGPU->viscosity_.components_), &(flowField.viscosity_.components_), sizeof(int), 0, 0, targetDevice, hostDevice);
    omp_target_memcpy(&(flowFieldGPU->viscosity_.size_), &(flowField.viscosity_.size_), sizeof(int), 0, 0, targetDevice, hostDevice);

    return viscosityDataGPU;
  }

  static RealType* mapDistanceToGPU(int hostDevice, int targetDevice, TurbulentFlowField& flowField, TurbulentFlowField* flowFieldGPU) {
    size_t scalarFieldDataSize = flowField.distance_.size_ * sizeof(RealType);

    RealType* distanceDataGPU = static_cast<RealType*>(omp_target_alloc(scalarFieldDataSize, targetDevice));
    if (!distanceDataGPU) {
      std::cout << "Error: Allocation failed for distanceDataGPU pointer." << std::endl;
      return nullptr;
    }

    bool copiedDistanceData = omp_target_memcpy(distanceDataGPU, flowField.distance_.data_, scalarFieldDataSize, 0, 0, targetDevice, hostDevice) == 0;
    if (!copiedDistanceData) {
      std::cout << "Error: Copying Distance data to GPU not successful." << std::endl;
      return nullptr;
    }

    bool copiedDistanceDataPtr = omp_target_memcpy(&flowFieldGPU->distance_.data_, &distanceDataGPU, sizeof(RealType*), 0, 0, targetDevice, hostDevice) == 0;
    if (!copiedDistanceDataPtr) {
      std::cout << "Error: Copying distance data pointer to distance object not successful." << std::endl;
      return nullptr;
    }

    omp_target_memcpy(&(flowFieldGPU->distance_.sizeX_), &(flowField.distance_.sizeX_), sizeof(int), 0, 0, targetDevice, hostDevice);
    omp_target_memcpy(&(flowFieldGPU->distance_.sizeY_), &(flowField.distance_.sizeY_), sizeof(int), 0, 0, targetDevice, hostDevice);
    omp_target_memcpy(&(flowFieldGPU->distance_.sizeZ_), &(flowField.distance_.sizeZ_), sizeof(int), 0, 0, targetDevice, hostDevice);
    omp_target_memcpy(&(flowFieldGPU->distance_.components_), &(flowField.distance_.components_), sizeof(int), 0, 0, targetDevice, hostDevice);
    omp_target_memcpy(&(flowFieldGPU->distance_.size_), &(flowField.distance_.size_), sizeof(int), 0, 0, targetDevice, hostDevice);

    return distanceDataGPU;
  }

public:
  ScalarField viscosity_;
  ScalarField distance_;
  /** Constructs a field from parameters object
   *
   * Constructs a field from a parameters object, so that it dimensionality can be defined in
   * the configuration file.
   *
   * @param parameters Parameters object with geometric information
   */
  TurbulentFlowField(const Parameters& parameters);

  virtual ~TurbulentFlowField() = default;

#pragma omp declare target
  ScalarField& getViscosity();
  ScalarField& getDistance();

  void getDistanceAndViscosity(RealType& distance, RealType& viscosity, int i, int j);
  void getDistanceAndViscosity(RealType& distance, RealType& viscosity, int i, int j, int k);
#pragma omp end declare target

  static TurbulentFlowFieldPtrs mapToGPU(int hostDevice, int targetDevice, TurbulentFlowField& flowField) {
    TurbulentFlowField* flowFieldGPU     = mapTurbulentFlowFieldToGPU(hostDevice, targetDevice, flowField);
    RealType*           rhsDataGPU       = mapRHSToGPU(hostDevice, targetDevice, flowField, flowFieldGPU);
    RealType*           pressureDataGPU  = mapPressureToGPU(hostDevice, targetDevice, flowField, flowFieldGPU);
    int*                flagsDataGPU     = mapFlagsToGPU(hostDevice, targetDevice, flowField, flowFieldGPU);
    RealType*           velocityDataGPU  = mapVelocityToGPU(hostDevice, targetDevice, flowField, flowFieldGPU);
    RealType*           fghDataGPU       = mapFGHToGPU(hostDevice, targetDevice, flowField, flowFieldGPU);
    RealType*           viscosityDataGPU = mapViscosityToGPU(hostDevice, targetDevice, flowField, flowFieldGPU);
    RealType*           distanceDataGPU  = mapDistanceToGPU(hostDevice, targetDevice, flowField, flowFieldGPU);

    return TurbulentFlowFieldPtrs(flowFieldGPU, pressureDataGPU, velocityDataGPU, rhsDataGPU, fghDataGPU, flagsDataGPU, viscosityDataGPU, distanceDataGPU);
  }

  static void mapToCPU(int hostDevice, int targetDevice, TurbulentFlowField& flowField, TurbulentFlowFieldPtrs ptrs) {
    size_t flowFieldSize          = sizeof(flowField);
    size_t scalarFieldDataSize    = flowField.RHS_.size_ * sizeof(RealType);
    size_t vectorFieldDataSize    = flowField.velocity_.size_ * sizeof(RealType);
    size_t intScalarFieldDataSize = flowField.flags_.size_ * sizeof(int);

    if (!omp_target_is_present(&flowField, targetDevice)) {
      std::cout << "Error: FlowField is not a valid target pointer." << std::endl;
      exit(EXIT_FAILURE);
    }

    bool copiedPressure = omp_target_memcpy(flowField.pressure_.data_, ptrs.pressureDataPtr_, scalarFieldDataSize, 0, 0, hostDevice, targetDevice) == 0;
    if (!copiedPressure) {
      std::cout << "Error: Copying pressure data from GPU to CPU not successful." << std::endl;
    }

    bool copiedRHS = omp_target_memcpy(flowField.RHS_.data_, ptrs.rhsDataPtr_, scalarFieldDataSize, 0, 0, hostDevice, targetDevice) == 0;
    if (!copiedRHS) {
      std::cout << "Error: Copying RHS data from GPU to CPU not successful." << std::endl;
    }

    bool copiedFlags = omp_target_memcpy(flowField.flags_.data_, ptrs.flagsDataPtr_, intScalarFieldDataSize, 0, 0, hostDevice, targetDevice) == 0;
    if (!copiedFlags) {
      std::cout << "Error: Copying flags data from GPU to CPU not successful." << std::endl;
    }

    bool copiedVelocity = omp_target_memcpy(flowField.velocity_.data_, ptrs.velocityDataPtr_, vectorFieldDataSize, 0, 0, hostDevice, targetDevice) == 0;
    if (!copiedVelocity) {
      std::cout << "Error: Copying velocity data from GPU to CPU not successful." << std::endl;
    }

    bool copiedFGH = omp_target_memcpy(flowField.FGH_.data_, ptrs.fghDataPtr_, vectorFieldDataSize, 0, 0, hostDevice, targetDevice) == 0;
    if (!copiedFGH) {
      std::cout << "Error: Copying FGH data from GPU to CPU not successful." << std::endl;
    }

    bool copiedViscosity = omp_target_memcpy(flowField.viscosity_.data_, ptrs.viscosityDataPtr_, scalarFieldDataSize, 0, 0, hostDevice, targetDevice) == 0;
    if (!copiedViscosity) {
      std::cout << "Error: Copying Viscosity data from GPU to CPU not successful." << std::endl;
    }

    bool copiedDistance = omp_target_memcpy(flowField.distance_.data_, ptrs.distanceDataPtr_, scalarFieldDataSize, 0, 0, hostDevice, targetDevice) == 0;
    if (!copiedDistance) {
      std::cout << "Error: Copying distance data from GPU to CPU not successful." << std::endl;
    }
  }

  static void mapToCPUAndFree(int hostDevice, int targetDevice, TurbulentFlowField& flowField, TurbulentFlowFieldPtrs& ptrs) {
    TurbulentFlowField::mapToCPU(hostDevice, targetDevice, flowField, ptrs);

    bool disassociatedFlowField = omp_target_disassociate_ptr(&flowField, targetDevice) == 0;
    if (!disassociatedFlowField) {
      std::cout << "Error: Could not disassociate FlowField pointer." << std::endl;
    }

    omp_target_free(ptrs.pressureDataPtr_, targetDevice);
    omp_target_free(ptrs.velocityDataPtr_, targetDevice);
    omp_target_free(ptrs.rhsDataPtr_, targetDevice);
    omp_target_free(ptrs.fghDataPtr_, targetDevice);
    omp_target_free(ptrs.flagsDataPtr_, targetDevice);
    omp_target_free(ptrs.viscosityDataPtr_, targetDevice);
    omp_target_free(ptrs.distanceDataPtr_, targetDevice);
    omp_target_free(ptrs.flowFieldPtr_, targetDevice);
  }

  static void mapToCPUVTK(int hostDevice, int targetDevice, TurbulentFlowField& flowField, TurbulentFlowFieldPtrs& ptrs) {
    size_t scalarFieldDataSize = flowField.RHS_.size_ * sizeof(RealType);
    size_t vectorFieldDataSize = flowField.velocity_.size_ * sizeof(RealType);

    if (!omp_target_is_present(&flowField, targetDevice)) {
      std::cout << "Error: FlowField is not a valid target pointer." << std::endl;
      exit(EXIT_FAILURE);
    }

    bool copiedPressure = omp_target_memcpy(flowField.pressure_.data_, ptrs.pressureDataPtr_, scalarFieldDataSize, 0, 0, hostDevice, targetDevice) == 0;
    if (!copiedPressure) {
      std::cout << "Error: Copying pressure data from GPU to CPU not successful." << std::endl;
    }

    bool copiedVelocity = omp_target_memcpy(flowField.velocity_.data_, ptrs.velocityDataPtr_, vectorFieldDataSize, 0, 0, hostDevice, targetDevice) == 0;
    if (!copiedVelocity) {
      std::cout << "Error: Copying velocity data from GPU to CPU not successful." << std::endl;
    }

    bool copiedViscosity = omp_target_memcpy(flowField.viscosity_.data_, ptrs.viscosityDataPtr_, scalarFieldDataSize, 0, 0, hostDevice, targetDevice) == 0;
    if (!copiedViscosity) {
      std::cout << "Error: Copying viscosity data from GPU to CPU not successful." << std::endl;
    }
  }
};
