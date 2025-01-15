#pragma once

#include "BoundaryType.hpp"
#include "Definitions.hpp"
#include "MeshsizeDelegate.hpp"
#include "MixingLengths.hpp"

enum ScenarioType { CAVITY, CHANNEL };
enum VelocityProfile { UNIFORM, PARABOLIC };

//! Classes for the parts of the parameters
//@{
class TimestepParameters {
public:
  RealType dt  = 0; //! Timestep
  RealType tau = 0; //! Security factor
};

class SimulationParameters {
public:
  RealType        finalTime = 0;             //! Final time for the simulation
  std::string     type;                      //! Type of the simulation (DNS vs. Turbulence)
  ScenarioType    scenario        = CAVITY;  //! If channel or cavity, for example
  VelocityProfile velocityProfile = UNIFORM; //! block or parabolic, block by default
};

class EnvironmentalParameters {
public:
  // Gravity components
  RealType gx = 0;
  RealType gy = 0;
  RealType gz = 0;
};

class FlowParameters {
public:
  RealType Re = 0; //! Reynolds number
};

class TurbulenceParameters {
public:
  TurbulenceType deltaMixLen = ZERO;
};

class SolverParameters {
public:
  RealType gamma         = 0;  //! Donor cell balance coefficient
  int      maxIterations = -1; //! Maximum number of iterations in the linear solver
};

class GeometricParameters {
public:
  // Dimensions
  int dim = -1;

  // Number of cells
  int sizeX = -1;
  int sizeY = -1;
  int sizeZ = -1;

  // Cell sizing
  RealType lengthX = 0;
  RealType lengthY = 0;
  RealType lengthZ = 0;

  // Meshsize type
  int meshsizeType = -1;

  // For mesh-stretching
  int stretchX = -1;
  int stretchY = -1;
  int stretchZ = -1;
};

class WallParameters {
public:
  // Scalar value definition. Used to define the pressure, for example.
  RealType scalarLeft;
  RealType scalarRight;
  RealType scalarBottom;
  RealType scalarTop;
  RealType scalarFront;
  RealType scalarBack;

  // Vector values at the boundaries, to define, for example, the velocities.
  RealType vectorLeft[3];
  RealType vectorRight[3];
  RealType vectorBottom[3];
  RealType vectorTop[3];
  RealType vectorFront[3];
  RealType vectorBack[3];

  // Define how will the boundary behave
  BoundaryType typeLeft;
  BoundaryType typeRight;
  BoundaryType typeTop;
  BoundaryType typeBottom;
  BoundaryType typeFront;
  BoundaryType typeBack;
};


class VTKParameters {
public:
  RealType    interval = 0; //! Time interval for file printing
  std::string prefix;       //! Output filename
};

class StdOutParameters {
public:
  RealType interval = 0;
};

class ParallelParameters {
public:
  int rank = -1; //! Rank of the current processor

  int numProcessors[3]; //! Array with the number of processors in each direction

  //@brief Ranks of the neighbours
  //@{
  int leftNb   = MPI_PROC_NULL;
  int rightNb  = MPI_PROC_NULL;
  int bottomNb = MPI_PROC_NULL;
  int topNb    = MPI_PROC_NULL;
  int frontNb  = MPI_PROC_NULL;
  int backNb   = MPI_PROC_NULL;
  //@}

  int indices[3];     //! 3D indices to locate the array
  int localSize[3];   //! Size for the local flow field
  int firstCorner[3]; //! Position of the first element. Used for plotting.

  PetscInt* sizes[3]; //! Arrays with the sizes of the blocks in each direction.
};

class BFStepParameters {
public:
  RealType xRatio = 0;
  RealType yRatio = 0;
};

//@}
class Parameters;

struct ParametersGPUPtrs {
  Parameters*       parametersGPUPtrs_;
  // Wall parameters
  RealType* vectorLeftGPUPtr_;
  RealType* vectorRightGPUPtr_;
  RealType* vectorBottomGPUPtr_;
  RealType* vectorTopGPUPtr_;
  RealType* vectorFrontGPUPtr_;
  RealType* vectorBackGPUPtr_;

  // Parallel parameters
  int*              numProcessorsGPUPtr_;
  int*              indicesGPUPtr_;
  int*              localSizeGPUPtr_;
  int*              firstCornerGPUPtr_;
  PetscInt**        sizesGPUPtr_;
  PetscInt*         sizesGPUPtr1_;
  PetscInt*         sizesGPUPtr2_;
  PetscInt*         sizesGPUPtr3_;
  MeshsizeDelegate* meshsizeDelegateGPUPtr_;

  // Constructor
  ParametersGPUPtrs(
    Parameters*       parametersGPUPtrs,
    RealType*         vectorLeftGPUPtr,
    RealType*         vectorRightGPUPtr,
    RealType*         vectorBottomGPUPtr,
    RealType*         vectorTopGPUPtr,
    RealType*         vectorFrontGPUPtr,
    RealType*         vectorBackGPUPtr,
    int*              numProcessorsGPUPtr,
    int*              indicesGPUPtr,
    int*              localSizeGPUPtr,
    int*              firstCornerGPUPtr,
    PetscInt**        sizesGPUPtr,
    PetscInt*         sizesGPUPtr1,
    PetscInt*         sizesGPUPtr2,
    PetscInt*         sizesGPUPtr3,
    MeshsizeDelegate* meshsizeDelegateGPUPtr
  ):
    parametersGPUPtrs_(parametersGPUPtrs),
    vectorLeftGPUPtr_(vectorLeftGPUPtr),
    vectorRightGPUPtr_(vectorRightGPUPtr),
    vectorBottomGPUPtr_(vectorBottomGPUPtr),
    vectorTopGPUPtr_(vectorTopGPUPtr),
    vectorFrontGPUPtr_(vectorFrontGPUPtr),
    vectorBackGPUPtr_(vectorBackGPUPtr),
    numProcessorsGPUPtr_(numProcessorsGPUPtr),
    indicesGPUPtr_(indicesGPUPtr),
    localSizeGPUPtr_(localSizeGPUPtr),
    firstCornerGPUPtr_(firstCornerGPUPtr),
    sizesGPUPtr_(sizesGPUPtr),
    sizesGPUPtr1_(sizesGPUPtr1),
    sizesGPUPtr2_(sizesGPUPtr2),
    sizesGPUPtr3_(sizesGPUPtr3),
    meshsizeDelegateGPUPtr_(meshsizeDelegateGPUPtr) {}

  // Constructor
};

/** A class to store and pass around the parameters
 */
class Parameters {
public:
  Parameters(const GeometricParameters& geometricParameters, const ParallelParameters& parallelParameters);
  ~Parameters() = default;

  Parameters(const Parameters& parameters):
    simulation(parameters.simulation),
    timestep(parameters.timestep),
    environment(parameters.environment),
    flow(parameters.flow),
    solver(parameters.solver),
    geometry(parameters.geometry),
    walls(parameters.walls),
    vtk(parameters.vtk),
    parallel(parameters.parallel),
    stdOut(parameters.stdOut),
    bfStep(parameters.bfStep),
    turbulence(parameters.turbulence),
    meshsize(MeshsizeDelegate(parameters.geometry, parameters.parallel)) {}

  SimulationParameters    simulation;
  TimestepParameters      timestep;
  EnvironmentalParameters environment;
  FlowParameters          flow;
  SolverParameters        solver;
  GeometricParameters     geometry;
  WallParameters          walls;
  VTKParameters           vtk;
  ParallelParameters      parallel;
  StdOutParameters        stdOut;
  BFStepParameters        bfStep;
  TurbulenceParameters    turbulence;
  MeshsizeDelegate        meshsize;

  static ParametersGPUPtrs mapToGPU(int hostDevice, int targetDevice, Parameters& parameters) {
    size_t      parametersSize = sizeof(parameters);
    size_t      realTypeSize   = sizeof(RealType);
    size_t      intSize        = sizeof(int);
    size_t      vectorDataSize = 3 * realTypeSize;
    Parameters* parametersGPU  = static_cast<Parameters*>(omp_target_alloc(parametersSize, targetDevice));
    if (!parametersGPU) {
      std::cout << "Failed to allocate memory for parameters on GPU" << std::endl;
    }

    bool associatedParameters = omp_target_associate_ptr(&parameters, parametersGPU, parametersSize, 0, targetDevice) == 0;
    if (!associatedParameters) {
      std::cout << "Failed to associate parameters with GPU memory" << std::endl;
    }

    RealType*         vectorLeftGPUPtr;
    RealType*         vectorRightGPUPtr;
    RealType*         vectorBottomGPUPtr;
    RealType*         vectorTopGPUPtr;
    RealType*         vectorFrontGPUPtr;
    RealType*         vectorBackGPUPtr;
    int*              numProcessorsGPUPtr;
    int*              indicesGPUPtr;
    int*              localSizeGPUPtr;
    int*              firstCornerGPUPtr;
    PetscInt**        sizesGPUPtr;
    PetscInt*         sizesGPUPtr1;
    PetscInt*         sizesGPUPtr2;
    PetscInt*         sizesGPUPtr3;
    MeshsizeDelegate* meshsizeGPUPtr;

    // SimulationParameters
    {
      omp_target_memcpy(&parametersGPU->simulation.finalTime, &parameters.simulation.finalTime, realTypeSize, 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&parametersGPU->simulation.scenario, &parameters.simulation.scenario, sizeof(ScenarioType), 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&parametersGPU->simulation.velocityProfile, &parameters.simulation.velocityProfile, sizeof(VelocityProfile), 0, 0, targetDevice, hostDevice);
    }

    // TimestepParameters
    {
      omp_target_memcpy(&parametersGPU->timestep.dt, &parameters.timestep.dt, realTypeSize, 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&parametersGPU->timestep.tau, &parameters.timestep.tau, realTypeSize, 0, 0, targetDevice, hostDevice);
    }

    // EnvironmentalParameters
    {
      omp_target_memcpy(&parametersGPU->environment.gx, &parameters.environment.gx, realTypeSize, 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&parametersGPU->environment.gy, &parameters.environment.gy, realTypeSize, 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&parametersGPU->environment.gz, &parameters.environment.gz, realTypeSize, 0, 0, targetDevice, hostDevice);
    }

    // FlowParameters
    { omp_target_memcpy(&parametersGPU->flow.Re, &parameters.flow.Re, realTypeSize, 0, 0, targetDevice, hostDevice); }

    // SolverParameters
    {
      omp_target_memcpy(&parametersGPU->solver.gamma, &parameters.solver.gamma, realTypeSize, 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&parametersGPU->solver.maxIterations, &parameters.solver.maxIterations, intSize, 0, 0, targetDevice, hostDevice);
    }

    // GeometricParameters
    {
      omp_target_memcpy(&parametersGPU->geometry.dim, &parameters.geometry.dim, intSize, 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&parametersGPU->geometry.sizeX, &parameters.geometry.sizeX, intSize, 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&parametersGPU->geometry.sizeY, &parameters.geometry.sizeY, intSize, 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&parametersGPU->geometry.sizeZ, &parameters.geometry.sizeZ, intSize, 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&parametersGPU->geometry.lengthX, &parameters.geometry.lengthX, realTypeSize, 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&parametersGPU->geometry.lengthY, &parameters.geometry.lengthY, realTypeSize, 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&parametersGPU->geometry.lengthZ, &parameters.geometry.lengthZ, realTypeSize, 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&parametersGPU->geometry.meshsizeType, &parameters.geometry.meshsizeType, intSize, 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&parametersGPU->geometry.stretchX, &parameters.geometry.stretchX, intSize, 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&parametersGPU->geometry.stretchY, &parameters.geometry.stretchY, intSize, 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&parametersGPU->geometry.stretchZ, &parameters.geometry.stretchZ, intSize, 0, 0, targetDevice, hostDevice);
    }

    // WallParameters
    {
      omp_target_memcpy(&parametersGPU->walls.scalarLeft, &parameters.walls.scalarLeft, realTypeSize, 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&parametersGPU->walls.scalarRight, &parameters.walls.scalarRight, realTypeSize, 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&parametersGPU->walls.scalarBottom, &parameters.walls.scalarBottom, realTypeSize, 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&parametersGPU->walls.scalarTop, &parameters.walls.scalarTop, realTypeSize, 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&parametersGPU->walls.scalarFront, &parameters.walls.scalarFront, realTypeSize, 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&parametersGPU->walls.scalarBack, &parameters.walls.scalarBack, realTypeSize, 0, 0, targetDevice, hostDevice);

      omp_target_memcpy(&parametersGPU->walls.typeLeft, &parameters.walls.typeLeft, sizeof(BoundaryType), 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&parametersGPU->walls.typeRight, &parameters.walls.typeRight, sizeof(BoundaryType), 0, 0, targetDevice, hostDevice);

      vectorLeftGPUPtr = static_cast<RealType*>(omp_target_alloc(vectorDataSize, targetDevice));
      if (!vectorLeftGPUPtr) {
        std::cout << "Error: Allocation failed for vectorLeftGPUPtr pointer." << std::endl;
      }
      bool copiedVectorLeftData = omp_target_memcpy(vectorLeftGPUPtr, parameters.walls.vectorLeft, vectorDataSize, 0, 0, targetDevice, hostDevice) == 0;
      if (!copiedVectorLeftData) {
        std::cout << "Error: Copying vectorLeft data to GPU not successful." << std::endl;
      }

      bool copiedVectorLeftPtr = omp_target_memcpy(&parametersGPU->walls.vectorLeft, vectorLeftGPUPtr, sizeof(RealType*), 0, 0, targetDevice, hostDevice) == 0;
      if (!copiedVectorLeftPtr) {
        std::cout << "Error: Copying vectorLeft data pointer to parameters object not successful." << std::endl;
      }

      vectorRightGPUPtr = static_cast<RealType*>(omp_target_alloc(vectorDataSize, targetDevice));
      if (!vectorRightGPUPtr) {
        std::cout << "Error: Allocation failed for vectorRightGPUPtr pointer." << std::endl;
      }
      bool copiedVectorRightData = omp_target_memcpy(vectorRightGPUPtr, parameters.walls.vectorRight, vectorDataSize, 0, 0, targetDevice, hostDevice) == 0;
      if (!copiedVectorRightData) {
        std::cout << "Error: Copying vectorRight data to GPU not successful." << std::endl;
      }
      bool copiedVectorRightPtr = omp_target_memcpy(&parametersGPU->walls.vectorRight, vectorRightGPUPtr, sizeof(RealType*), 0, 0, targetDevice, hostDevice) == 0;
      if (!copiedVectorRightPtr) {
        std::cout << "Error: Copying vectorRight data pointer to parameters object not successful." << std::endl;
      }


      vectorTopGPUPtr = static_cast<RealType*>(omp_target_alloc(vectorDataSize, targetDevice));
      if (!vectorTopGPUPtr) {
        std::cout << "Error: Allocation failed for vectorTopGPUPtr pointer." << std::endl;
      }
      bool copiedVectorTopData = omp_target_memcpy(vectorTopGPUPtr, parameters.walls.vectorTop, vectorDataSize, 0, 0, targetDevice, hostDevice) == 0;
      if (!copiedVectorTopData) {
        std::cout << "Error: Copying vectorTop data to GPU not successful." << std::endl;
      }
      bool copiedVectorTopPtr = omp_target_memcpy(&parametersGPU->walls.vectorTop, vectorTopGPUPtr, sizeof(RealType*), 0, 0, targetDevice, hostDevice) == 0;
      if (!copiedVectorTopPtr) {
        std::cout << "Error: Copying vectorTop data pointer to parameters object not successful." << std::endl;
      }


      vectorBottomGPUPtr = static_cast<RealType*>(omp_target_alloc(vectorDataSize, targetDevice));
      if (!vectorBottomGPUPtr) {
        std::cout << "Error: Allocation failed for vectorBottomGPUPtr pointer." << std::endl;
      }
      bool copiedVectorBottomData = omp_target_memcpy(vectorBottomGPUPtr, parameters.walls.vectorBottom, vectorDataSize, 0, 0, targetDevice, hostDevice) == 0;
      if (!copiedVectorBottomData) {
        std::cout << "Error: Copying vectorBottom data to GPU not successful." << std::endl;
      }
      bool copiedVectorBottomPtr = omp_target_memcpy(&parametersGPU->walls.vectorBottom, vectorBottomGPUPtr, sizeof(RealType*), 0, 0, targetDevice, hostDevice) == 0;
      if (!copiedVectorBottomPtr) {
        std::cout << "Error: Copying vectorBottom data pointer to parameters object not successful." << std::endl;
      }


      vectorFrontGPUPtr = static_cast<RealType*>(omp_target_alloc(vectorDataSize, targetDevice));
      if (!vectorFrontGPUPtr) {
        std::cout << "Error: Allocation failed for vectorFrontGPUPtr pointer." << std::endl;
      }
      bool copiedVectorFrontData = omp_target_memcpy(vectorFrontGPUPtr, parameters.walls.vectorFront, vectorDataSize, 0, 0, targetDevice, hostDevice) == 0;
      if (!copiedVectorFrontData) {
        std::cout << "Error: Copying vectorFront data to GPU not successful." << std::endl;
      }
      bool copiedVectorFrontPtr = omp_target_memcpy(&parametersGPU->walls.vectorFront, vectorFrontGPUPtr, sizeof(RealType*), 0, 0, targetDevice, hostDevice) == 0;
      if (!copiedVectorFrontPtr) {
        std::cout << "Error: Copying vectorFront data pointer to parameters object not successful." << std::endl;
      }


      vectorBackGPUPtr = static_cast<RealType*>(omp_target_alloc(vectorDataSize, targetDevice));
      if (!vectorBackGPUPtr) {
        std::cout << "Error: Allocation failed for vectorBackGPUPtr pointer." << std::endl;
      }
      bool copiedVectorBackData = omp_target_memcpy(vectorBackGPUPtr, parameters.walls.vectorBack, vectorDataSize, 0, 0, targetDevice, hostDevice) == 0;
      if (!copiedVectorBackData) {
        std::cout << "Error: Copying vectorBack data to GPU not successful." << std::endl;
      }
      bool copiedVectorBackPtr = omp_target_memcpy(&parametersGPU->walls.vectorBack, vectorBackGPUPtr, sizeof(RealType*), 0, 0, targetDevice, hostDevice) == 0;
      if (!copiedVectorBackPtr) {
        std::cout << "Error: Copying vectorBack data pointer to parameters object not successful." << std::endl;
      }
    }

    // ParallelParameters
    {
      omp_target_memcpy(&parametersGPU->parallel.rank, &parameters.parallel.rank, intSize, 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&parametersGPU->parallel.leftNb, &parameters.parallel.leftNb, intSize, 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&parametersGPU->parallel.rightNb, &parameters.parallel.rightNb, intSize, 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&parametersGPU->parallel.bottomNb, &parameters.parallel.bottomNb, intSize, 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&parametersGPU->parallel.topNb, &parameters.parallel.topNb, intSize, 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&parametersGPU->parallel.frontNb, &parameters.parallel.frontNb, intSize, 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&parametersGPU->parallel.backNb, &parameters.parallel.backNb, intSize, 0, 0, targetDevice, hostDevice);

      numProcessorsGPUPtr = static_cast<int*>(omp_target_alloc(sizeof(int) * 3, targetDevice));
      if (!numProcessorsGPUPtr) {
        std::cout << "Error: Allocation failed for numProcessorsGPUPtr pointer." << std::endl;
      }
      bool copiedNumProcessorsData = omp_target_memcpy(numProcessorsGPUPtr, parameters.parallel.numProcessors, sizeof(int) * 3, 0, 0, targetDevice, hostDevice) == 0;
      if (!copiedNumProcessorsData) {
        std::cout << "Error: Copying copiedNumProcessorsData data to GPU not successful." << std::endl;
      }
      bool copiedNumProcessorsPtr = omp_target_memcpy(&parametersGPU->parallel.numProcessors, numProcessorsGPUPtr, sizeof(int*), 0, 0, targetDevice, hostDevice) == 0;
      if (!copiedNumProcessorsPtr) {
        std::cout << "Error: Copying numProcessors data pointer to parameters object not successful." << std::endl;
      }


      indicesGPUPtr = static_cast<int*>(omp_target_alloc(sizeof(int) * 3, targetDevice));
      if (!indicesGPUPtr) {
        std::cout << "Error: Allocation failed for indicesGPUPtr pointer." << std::endl;
      }
      bool copiedIndicesData = omp_target_memcpy(indicesGPUPtr, parameters.parallel.indices, sizeof(int) * 3, 0, 0, targetDevice, hostDevice) == 0;
      if (!copiedIndicesData) {
        std::cout << "Error: Copying copiedIndicesData data to GPU not successful." << std::endl;
      }
      bool copiedIndicesPtr = omp_target_memcpy(&parametersGPU->parallel.indices, indicesGPUPtr, sizeof(int*), 0, 0, targetDevice, hostDevice) == 0;
      if (!copiedIndicesPtr) {
        std::cout << "Error: Copying indices data pointer to parameters object not successful." << std::endl;
      }


      localSizeGPUPtr = static_cast<int*>(omp_target_alloc(sizeof(int) * 3, targetDevice));
      if (!localSizeGPUPtr) {
        std::cout << "Error: Allocation failed for localSizeGPUPtr pointer." << std::endl;
      }
      bool copiedLocalSizeData = omp_target_memcpy(localSizeGPUPtr, parameters.parallel.localSize, sizeof(int) * 3, 0, 0, targetDevice, hostDevice) == 0;
      if (!copiedLocalSizeData) {
        std::cout << "Error: Copying copiedLocalSizeData data to GPU not successful." << std::endl;
      }
      bool copiedLocalSizePtr = omp_target_memcpy(&parametersGPU->parallel.localSize, localSizeGPUPtr, sizeof(int*), 0, 0, targetDevice, hostDevice) == 0;
      if (!copiedLocalSizePtr) {
        std::cout << "Error: Copying localSize data pointer to parameters object not successful." << std::endl;
      }


      firstCornerGPUPtr = static_cast<int*>(omp_target_alloc(sizeof(int) * 3, targetDevice));
      if (!firstCornerGPUPtr) {
        std::cout << "Error: Allocation failed for firstCornerGPUPtr pointer." << std::endl;
      }
      bool copiedFirstCornerData = omp_target_memcpy(firstCornerGPUPtr, parameters.parallel.firstCorner, sizeof(int) * 3, 0, 0, targetDevice, hostDevice) == 0;
      if (!copiedFirstCornerData) {
        std::cout << "Error: Copying copiedFirstCornerData data to GPU not successful." << std::endl;
      }
      bool copiedFirstCornerPtr = omp_target_memcpy(&parametersGPU->parallel.firstCorner, firstCornerGPUPtr, sizeof(int*), 0, 0, targetDevice, hostDevice) == 0;
      if (!copiedFirstCornerPtr) {
        std::cout << "Error: Copying firstCorner data pointer to parameters object not successful." << std::endl;
      }

      sizesGPUPtr  = static_cast<PetscInt**>(omp_target_alloc(sizeof(PetscInt*) * 3, targetDevice));
      sizesGPUPtr1 = static_cast<PetscInt*>(omp_target_alloc(sizeof(PetscInt) * 1, targetDevice));
      sizesGPUPtr2 = static_cast<PetscInt*>(omp_target_alloc(sizeof(PetscInt) * 1, targetDevice));
      sizesGPUPtr3 = static_cast<PetscInt*>(omp_target_alloc(sizeof(PetscInt) * 1, targetDevice));

      bool copiedSizesGPUData1 = omp_target_memcpy(sizesGPUPtr1, parameters.parallel.sizes[0], sizeof(PetscInt), 0, 0, targetDevice, hostDevice) == 0;
      if (!copiedSizesGPUData1) {
        std::cout << "Error: Copying copiedSizesGPUData1 data to GPU not successful." << std::endl;
      }

      bool copiedSizesGPUData2 = omp_target_memcpy(sizesGPUPtr2, parameters.parallel.sizes[1], sizeof(PetscInt), 0, 0, targetDevice, hostDevice) == 0;
      if (!copiedSizesGPUData2) {
        std::cout << "Error: Copying copiedSizesGPUData2 data to GPU not successful." << std::endl;
      }

      bool copiedSizesGPUPtr1 = omp_target_memcpy(&sizesGPUPtr[0], &sizesGPUPtr1, sizeof(PetscInt*), 0, 0, targetDevice, hostDevice) == 0;
      if (!copiedSizesGPUPtr1) {
        std::cout << "Error: Copying copiedSizesGPUPtr1 pointer to GPU not successful." << std::endl;
      }

      bool copiedSizesGPUPtr2 = omp_target_memcpy(&sizesGPUPtr[1], &sizesGPUPtr2, sizeof(PetscInt*), 0, 0, targetDevice, hostDevice) == 0;
      if (!copiedSizesGPUPtr2) {
        std::cout << "Error: Copying copiedSizesGPUPtr2 pointer to GPU not successful." << std::endl;
      }

      if (parameters.geometry.dim == 3) {
        bool copiedSizesGPUData3 = omp_target_memcpy(sizesGPUPtr3, parameters.parallel.sizes[2], sizeof(PetscInt), 0, 0, targetDevice, hostDevice) == 0;
        if (!copiedSizesGPUData3) {
          std::cout << "Error: Copying copiedSizesGPUData3 data to GPU not successful." << std::endl;
        }
        bool copiedSizesGPUPtr3 = omp_target_memcpy(&sizesGPUPtr[2], &sizesGPUPtr3, sizeof(PetscInt*), 0, 0, targetDevice, hostDevice) == 0;
        if (!copiedSizesGPUPtr3) {
          std::cout << "Error: Copying copiedSizesGPUPtr3 pointer to GPU not successful." << std::endl;
        }
      }
      bool copiedSizesGPUPtr = omp_target_memcpy(&parametersGPU->parallel.sizes, sizesGPUPtr, sizeof(PetscInt**) * 3, 0, 0, targetDevice, hostDevice) == 0;
      if (!copiedSizesGPUPtr) {
        std::cout << "Error: Copying sizesGPU data pointer to parameters object not successful." << std::endl;
      }
    }

    // StdOutParameters
    { omp_target_memcpy(&parametersGPU->stdOut.interval, &parameters.stdOut.interval, realTypeSize, 0, 0, targetDevice, hostDevice); }

    // BFStepParameters
    {
      omp_target_memcpy(&parametersGPU->bfStep.xRatio, &parameters.bfStep.xRatio, realTypeSize, 0, 0, targetDevice, hostDevice);
      omp_target_memcpy(&parametersGPU->bfStep.yRatio, &parameters.bfStep.yRatio, realTypeSize, 0, 0, targetDevice, hostDevice);
    }

    // TurbulenceParameters
    { omp_target_memcpy(&parametersGPU->turbulence.deltaMixLen, &parameters.turbulence.deltaMixLen, sizeof(TurbulenceType), 0, 0, targetDevice, hostDevice); }

    // MeshsizeDelegate
    {
      meshsizeGPUPtr = MeshsizeDelegate::mapToGPU(hostDevice, targetDevice, parameters.meshsize);
      omp_target_memcpy(&parametersGPU->meshsize, meshsizeGPUPtr, sizeof(MeshsizeDelegate*), 0, 0, targetDevice, hostDevice);
    }

    return ParametersGPUPtrs(
      parametersGPU,
      vectorLeftGPUPtr,
      vectorRightGPUPtr,
      vectorBottomGPUPtr,
      vectorTopGPUPtr,
      vectorFrontGPUPtr,
      vectorBackGPUPtr,
      numProcessorsGPUPtr,
      indicesGPUPtr,
      localSizeGPUPtr,
      firstCornerGPUPtr,
      sizesGPUPtr,
      sizesGPUPtr1,
      sizesGPUPtr2,
      sizesGPUPtr3,
      meshsizeGPUPtr
    );
  }

  static void freeGPU(int hostDevice, int targetDevice, Parameters& parameters, ParametersGPUPtrs& parameterPtrs) {

    bool disassociatedParameters = omp_target_disassociate_ptr(&parameters, targetDevice) == 0;
    if (!disassociatedParameters) {
      std::cout << "Error: Could not disassociate Parameters pointer." << std::endl;
    }

    omp_target_free(parameterPtrs.firstCornerGPUPtr_, targetDevice);
    omp_target_free(parameterPtrs.indicesGPUPtr_, targetDevice);
    omp_target_free(parameterPtrs.localSizeGPUPtr_, targetDevice);
    MeshsizeDelegate::freeGPU(hostDevice, targetDevice, parameterPtrs.meshsizeDelegateGPUPtr_, parameters.meshsize);
    omp_target_free(parameterPtrs.numProcessorsGPUPtr_, targetDevice);
    omp_target_free(parameterPtrs.sizesGPUPtr1_, targetDevice);
    omp_target_free(parameterPtrs.sizesGPUPtr2_, targetDevice);
    omp_target_free(parameterPtrs.sizesGPUPtr3_, targetDevice);
    omp_target_free(parameterPtrs.sizesGPUPtr_, targetDevice);
    omp_target_free(parameterPtrs.vectorBackGPUPtr_, targetDevice);
    omp_target_free(parameterPtrs.vectorBottomGPUPtr_, targetDevice);
    omp_target_free(parameterPtrs.vectorFrontGPUPtr_, targetDevice);
    omp_target_free(parameterPtrs.vectorLeftGPUPtr_, targetDevice);
    omp_target_free(parameterPtrs.vectorRightGPUPtr_, targetDevice);
    omp_target_free(parameterPtrs.vectorTopGPUPtr_, targetDevice);
    omp_target_free(parameterPtrs.parametersGPUPtrs_, targetDevice);
  }
};
/* #pragma omp declare mapper(Parameters p) map(to : p) map(to : p.meshsize, p.meshsize.uniform, p.meshsize.tanh, p.meshsize.tanh.uniformMeshsize_) map(to : p.simulation \
) map(to : p.timestep) map(to : p.environment) map(to : p.flow) map(to : p.solver) map(to : p.geometry) \
  map(to : p.walls, \
        p.walls.vectorLeft[0 : 3], \
        p.walls.vectorRight[0 : 3], \
        p.walls.vectorBottom[0 : 3], \
        p.walls.vectorTop[0 : 3], \
        p.walls.vectorFront[0 : 3], \
        p.walls.vectorBack[0 : 3]) map(to : p.turbulence) */
