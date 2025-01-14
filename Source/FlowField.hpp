#pragma once

#include "DataStructures.hpp"
#include "Parameters.hpp"

class FlowField;

struct FlowFieldGPUPtrs {
  FlowField* flowFieldPtr_;
  RealType*  pressureDataPtr_;
  RealType*  velocityDataPtr_;
  RealType*  rhsDataPtr_;
  RealType*  fghDataPtr_;
  int*       flagsDataPtr_;

  FlowFieldGPUPtrs(FlowField* flowFieldPtr, RealType* pressureDataPtr, RealType* velocityDataPtr, RealType* rhsDataPtr, RealType* fghDataPtr, int* flagsDataPtr) {
    flowFieldPtr_    = flowFieldPtr;
    pressureDataPtr_ = pressureDataPtr;
    velocityDataPtr_ = velocityDataPtr;
    rhsDataPtr_      = rhsDataPtr;
    fghDataPtr_      = fghDataPtr;
    flagsDataPtr_    = flagsDataPtr;
  }
};

/** Flow field
 *
 * Class intended to contain the state of the domain.
 */
class FlowField {
public:
  int sizeX_; //! Size in the X direction
  int sizeY_; //! Size in the Y direction
  int sizeZ_; //! Size in the Z direction
  int cellsX_;
  int cellsY_;
  int cellsZ_;

  ScalarField pressure_; //! Scalar field representing the pressure
  VectorField velocity_; //! Multicomponent field representing velocity

  IntScalarField flags_; //! Integer field for the flags

  VectorField FGH_;
  ScalarField RHS_; //! Right hand side for the Poisson equation

  /** Constructor for the 2D flow field
   *
   * Constructor for the flow field. Allocates all the fields and sets
   * the sizes. Currently, this contructor is only used for testing purposes.
   *
   * @param Nx Size of the fuild domain (non-ghost cells), in the X direction
   * @param Ny Size of the fuild domain (non-ghost cells), in the Y direction
   */
  FlowField(int Nx, int Ny);

  /** Constructor for the 3D flow field
   *
   * Constructor for the flow field. Allocates all the fields and sets
   * the sizes. Currently, this contructor is only used for testing purposes.
   *
   * @param Nx Size of the fuild domain (non-ghost cells), in the X direction
   * @param Ny Size of the fuild domain (non-ghost cells), in the Y direction
   * @param Nz Size of the fuild domain (non-ghost cells), in the Z direction
   */
  FlowField(int Nx, int Ny, int Nz);

  /** Constructs a field from parameters object
   *
   * Constructs a field from a parameters object, so that it dimensionality can be defined in
   * the configuration file.
   *
   * @param parameters Parameters object with geometric information
   */
  FlowField(const Parameters& parameters);

  virtual ~FlowField() = default;

  FlowField(FlowField& flowField):
    sizeX_(flowField.getNx()),
    sizeY_(flowField.getNy()),
    sizeZ_(flowField.getNz()),
    cellsX_(flowField.getCellsX()),
    cellsY_(flowField.getCellsY()),
    cellsZ_(flowField.getCellsZ()),
    pressure_(flowField.getPressure()),
    velocity_(flowField.getVelocity()),
    flags_(flowField.getFlags()),
    FGH_(flowField.getFGH()),
    RHS_(flowField.getRHS()) {}

#pragma omp declare target

  /** Obtain size in the X direction
   *
   * @return Number of cells in the X direction
   */
  int getNx() const;

  /** Obtain size in the Y direction
   *
   * @return Number of cells in the Y direction
   */
  int getNy() const;

  /** Obtain size in the Z direction
   *
   * @return Number of cells in the Z direction
   */
  int getNz() const;

  int getCellsX() const;
  int getCellsY() const;
  int getCellsZ() const;

  ScalarField& getPressure();
  VectorField& getVelocity();

  IntScalarField& getFlags();

  VectorField& getFGH();

  ScalarField& getRHS();

  void getPressureAndVelocity(RealType& pressure, RealType* const velocity, int i, int j);
  void getPressureAndVelocity(RealType& pressure, RealType* const velocity, int i, int j, int k);

  static FlowFieldGPUPtrs mapToGPU(int hostDevice, int targetDevice, FlowField& flowField) {
    size_t flowFieldSize          = sizeof(flowField);
    size_t scalarFieldDataSize    = flowField.RHS_.size_ * sizeof(RealType);
    size_t vectorFieldDataSize    = flowField.velocity_.size_ * sizeof(RealType);
    size_t intScalarFieldDataSize = flowField.flags_.size_ * sizeof(int);

    RealType* rhsDataGPU = static_cast<RealType*>(omp_target_alloc(scalarFieldDataSize, targetDevice));
    if (!rhsDataGPU) {
      std::cout << "Error: Allocation failed for rhsDataGPU pointer." << std::endl;
    }

    RealType* pressureDataGPU = static_cast<RealType*>(omp_target_alloc(scalarFieldDataSize, targetDevice));
    if (!pressureDataGPU) {
      std::cout << "Error: Allocation failed for pressureDataGPU pointer." << std::endl;
    }

    int* flagsDataGPU = static_cast<int*>(omp_target_alloc(intScalarFieldDataSize, targetDevice));
    if (!flagsDataGPU) {
      std::cout << "Error: Allocation failed for flagsDataGPU pointer." << std::endl;
    }

    RealType* velocityDataGPU = static_cast<RealType*>(omp_target_alloc(vectorFieldDataSize, targetDevice));
    if (!velocityDataGPU) {
      std::cout << "Error: Allocation failed for velocityDataGPU pointer." << std::endl;
    }

    RealType* fghDataGPU = static_cast<RealType*>(omp_target_alloc(vectorFieldDataSize, targetDevice));
    if (!fghDataGPU) {
      std::cout << "Error: Allocation failed for fghDataGPU pointer." << std::endl;
    }

    FlowField* flowFieldGPU = static_cast<FlowField*>(omp_target_alloc(flowFieldSize, targetDevice));
    if (!flowFieldGPU) {
      std::cout << "Error: Allocation failed for flowFieldGPU pointer." << std::endl;
    }

    bool associatedFlowField = omp_target_associate_ptr(&flowField, flowFieldGPU, flowFieldSize, 0, targetDevice) == 0;
    if (!associatedFlowField) {
      std::cout << "Error: FlowField could not be associated to GPU pointer." << std::endl;
    }

    {
      bool copiedPressureData = omp_target_memcpy(pressureDataGPU, flowField.pressure_.data_, scalarFieldDataSize, 0, 0, targetDevice, hostDevice) == 0;
      if (!copiedPressureData) {
        std::cout << "Error: Copying pressure data to GPU not successful." << std::endl;
      }

      bool copiedRHSData = omp_target_memcpy(rhsDataGPU, flowField.RHS_.data_, scalarFieldDataSize, 0, 0, targetDevice, hostDevice) == 0;
      if (!copiedRHSData) {
        std::cout << "Error: Copying RHS data to GPU not successful." << std::endl;
      }

      bool copiedFGHData = omp_target_memcpy(fghDataGPU, flowField.FGH_.data_, vectorFieldDataSize, 0, 0, targetDevice, hostDevice) == 0;
      if (!copiedFGHData) {
        std::cout << "Error: Copying FGH data to GPU not successful." << std::endl;
      }

      bool copiedFlagsData = omp_target_memcpy(flagsDataGPU, flowField.flags_.data_, intScalarFieldDataSize, 0, 0, targetDevice, hostDevice) == 0;
      if (!copiedFlagsData) {
        std::cout << "Error: Copying flags data to GPU not successful." << std::endl;
      }

      bool copiedVelocityData = omp_target_memcpy(velocityDataGPU, flowField.velocity_.data_, vectorFieldDataSize, 0, 0, targetDevice, hostDevice) == 0;
      if (!copiedVelocityData) {
        std::cout << "Error: Copying velocity data to GPU not successful." << std::endl;
      }
    }

    {
      bool copiedPressureDataPtr = omp_target_memcpy(&flowFieldGPU->pressure_.data_, &pressureDataGPU, sizeof(RealType*), 0, 0, targetDevice, hostDevice) == 0;
      if (!copiedPressureDataPtr) {
        std::cout << "Error: Copying pressure data pointer to pressure object not successful." << std::endl;
      }

      bool copiedRHSDataPtr = omp_target_memcpy(&flowFieldGPU->RHS_.data_, &rhsDataGPU, sizeof(RealType*), 0, 0, targetDevice, hostDevice) == 0;
      if (!copiedRHSDataPtr) {
        std::cout << "Error: Copying RHS data pointer to rhs object not successful." << std::endl;
      }

      bool copiedFGHDataPtr = omp_target_memcpy(&flowFieldGPU->FGH_.data_, &fghDataGPU, sizeof(RealType*), 0, 0, targetDevice, hostDevice) == 0;
      if (!copiedFGHDataPtr) {
        std::cout << "Error: Copying FGH data pointer to fgh object not successful." << std::endl;
      }

      bool copiedFlagsDataPtr = omp_target_memcpy(&flowFieldGPU->flags_.data_, &flagsDataGPU, sizeof(int*), 0, 0, targetDevice, hostDevice) == 0;
      if (!copiedFlagsDataPtr) {
        std::cout << "Error: Copying flags data pointer to flags object not successful." << std::endl;
      }

      bool copiedVelocityDataPtr = omp_target_memcpy(&flowFieldGPU->velocity_.data_, &velocityDataGPU, sizeof(RealType*), 0, 0, targetDevice, hostDevice) == 0;
      if (!copiedVelocityDataPtr) {
        std::cout << "Error: Copying velocity data pointer to velocity object not successful." << std::endl;
      }
    }

    {
      omp_target_memcpy(&(flowFieldGPU->cellsX_), &(flowField.cellsX_), sizeof(int), 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&(flowFieldGPU->cellsY_), &(flowField.cellsY_), sizeof(int), 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&(flowFieldGPU->cellsZ_), &(flowField.cellsZ_), sizeof(int), 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&(flowFieldGPU->sizeX_), &(flowField.sizeX_), sizeof(int), 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&(flowFieldGPU->sizeY_), &(flowField.sizeY_), sizeof(int), 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&(flowFieldGPU->sizeZ_), &(flowField.sizeZ_), sizeof(int), 0, 0, targetDevice, hostDevice);

      omp_target_memcpy(&(flowFieldGPU->flags_.sizeX_), &(flowField.flags_.sizeX_), sizeof(int), 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&(flowFieldGPU->flags_.sizeY_), &(flowField.flags_.sizeY_), sizeof(int), 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&(flowFieldGPU->flags_.sizeZ_), &(flowField.flags_.sizeZ_), sizeof(int), 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&(flowFieldGPU->flags_.components_), &(flowField.flags_.components_), sizeof(int), 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&(flowFieldGPU->flags_.size_), &(flowField.flags_.size_), sizeof(int), 0, 0, targetDevice, hostDevice);

      omp_target_memcpy(&(flowFieldGPU->FGH_.sizeX_), &(flowField.FGH_.sizeX_), sizeof(int), 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&(flowFieldGPU->FGH_.sizeY_), &(flowField.FGH_.sizeY_), sizeof(int), 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&(flowFieldGPU->FGH_.sizeZ_), &(flowField.FGH_.sizeZ_), sizeof(int), 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&(flowFieldGPU->FGH_.components_), &(flowField.FGH_.components_), sizeof(int), 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&(flowFieldGPU->FGH_.size_), &(flowField.FGH_.size_), sizeof(int), 0, 0, targetDevice, hostDevice);

      omp_target_memcpy(&(flowFieldGPU->velocity_.sizeX_), &(flowField.velocity_.sizeX_), sizeof(int), 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&(flowFieldGPU->velocity_.sizeY_), &(flowField.velocity_.sizeY_), sizeof(int), 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&(flowFieldGPU->velocity_.sizeZ_), &(flowField.velocity_.sizeZ_), sizeof(int), 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&(flowFieldGPU->velocity_.components_), &(flowField.velocity_.components_), sizeof(int), 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&(flowFieldGPU->velocity_.size_), &(flowField.velocity_.size_), sizeof(int), 0, 0, targetDevice, hostDevice);

      omp_target_memcpy(&(flowFieldGPU->RHS_.sizeX_), &(flowField.RHS_.sizeX_), sizeof(int), 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&(flowFieldGPU->RHS_.sizeY_), &(flowField.RHS_.sizeY_), sizeof(int), 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&(flowFieldGPU->RHS_.sizeZ_), &(flowField.RHS_.sizeZ_), sizeof(int), 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&(flowFieldGPU->RHS_.components_), &(flowField.RHS_.components_), sizeof(int), 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&(flowFieldGPU->RHS_.size_), &(flowField.RHS_.size_), sizeof(int), 0, 0, targetDevice, hostDevice);

      omp_target_memcpy(&(flowFieldGPU->pressure_.sizeX_), &(flowField.pressure_.sizeX_), sizeof(int), 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&(flowFieldGPU->pressure_.sizeY_), &(flowField.pressure_.sizeY_), sizeof(int), 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&(flowFieldGPU->pressure_.sizeZ_), &(flowField.pressure_.sizeZ_), sizeof(int), 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&(flowFieldGPU->pressure_.components_), &(flowField.pressure_.components_), sizeof(int), 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&(flowFieldGPU->pressure_.size_), &(flowField.pressure_.size_), sizeof(int), 0, 0, targetDevice, hostDevice);
    }

    return FlowFieldGPUPtrs(flowFieldGPU, pressureDataGPU, velocityDataGPU, rhsDataGPU, fghDataGPU, flagsDataGPU);
  }

  static void mapToCPU(int hostDevice, int targetDevice, FlowField& flowField, FlowFieldGPUPtrs ptrs) {
    size_t flowFieldSize          = sizeof(flowField);
    size_t scalarFieldDataSize    = flowField.RHS_.size_ * sizeof(RealType);
    size_t vectorFieldDataSize    = flowField.velocity_.size_ * sizeof(RealType);
    size_t intScalarFieldDataSize = flowField.flags_.size_ * sizeof(int);

    std::cout << "FlowField size: " << flowFieldSize << " bytes" << std::endl;
    std::cout << "ScalarField data size: " << scalarFieldDataSize << " bytes" << std::endl;
    std::cout << "VectorField data size: " << vectorFieldDataSize << " bytes" << std::endl;
    std::cout << "IntScalarField data size: " << intScalarFieldDataSize << " bytes" << std::endl;

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
  }

  static void mapToCPUAndFree(int hostDevice, int targetDevice, FlowField& flowField, FlowFieldGPUPtrs ptrs) {
    FlowField::mapToCPU(hostDevice, targetDevice, flowField, ptrs);

    bool disassociatedFlowField = omp_target_disassociate_ptr(&flowField, targetDevice) == 0;
    if (!disassociatedFlowField) {
      std::cout << "Error: Could not disassociate FlowField pointer." << std::endl;
    }

    omp_target_free(ptrs.pressureDataPtr_, targetDevice);
    omp_target_free(ptrs.velocityDataPtr_, targetDevice);
    omp_target_free(ptrs.rhsDataPtr_, targetDevice);
    omp_target_free(ptrs.fghDataPtr_, targetDevice);
    omp_target_free(ptrs.flagsDataPtr_, targetDevice);
    omp_target_free(ptrs.flowFieldPtr_, targetDevice);
  }
#pragma omp end declare target
};
