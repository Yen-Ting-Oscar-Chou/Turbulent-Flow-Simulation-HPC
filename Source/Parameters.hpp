#pragma once

#include "BoundaryType.hpp"
#include "Definitions.hpp"
#include "Meshsize.hpp"

//! Classes for the parts of the parameters
//@{
class TimestepParameters {
public:
  RealType dt  = 0; //! Timestep
  RealType tau = 0; //! Security factor
};
#pragma omp declare mapper(TimestepParameters p)  \
 map(to: p, p.dt, p.tau)

class SimulationParameters {
public:
  RealType    finalTime = 0; //! Final time for the simulation
  std::string type;          //! Type of the simulation (DNS vs. Turbulence)
  std::string scenario;      //! If channel or cavity, for example
  std::string velocityProfile = "uniform"; //! block or parabolic, block by default
};
#pragma omp declare mapper(SimulationParameters p)  \
 map(to: p, p.finalTime, p.type, p.scenario, p.velocityProfile)

class EnvironmentalParameters {
public:
  // Gravity components
  RealType gx = 0;
  RealType gy = 0;
  RealType gz = 0;
};
#pragma omp declare mapper(EnvironmentalParameters p)  \
 map(to: p, p.gx, p.gy, p.gz)

class FlowParameters {
public:
  RealType Re = 0; //! Reynolds number
};
#pragma omp declare mapper(FlowParameters p)  \
 map(to: p, p.Re)

class TurbulenceParameters {
  public: 
    std::function<RealType (const RealType, const RealType)> deltaMixLen = NULL;
};
#pragma omp declare mapper(TurbulenceParameters p)  \
 map(to: p, p.deltaMixLen)

class SolverParameters {
public:
  RealType gamma         = 0;  //! Donor cell balance coefficient
  int      maxIterations = -1; //! Maximum number of iterations in the linear solver
};
#pragma omp declare mapper(SolverParameters p)  \
 map(to: p, p.gamma, p.maxIterations)

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
#pragma omp declare mapper(GeometricParameters p)  \
 map(to: p, p.dim, p.sizeX, p.sizeY, p.sizeZ,p.lengthX,p.lengthY,p.lengthZ,p.meshsizeType,p.stretchX ,p.stretchY,p.stretchZ)

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
#pragma omp declare mapper(WallParameters p)  \
 map(to: p, p.scalarLeft, p.scalarRight, p.scalarBottom, p.scalarTop, p.scalarFront, p.scalarBack, p.typeLeft, p.typeRight, p.typeTop, p.typeBottom, p.typeFront, p.typeBack, p.vectorLeft[0:3], p.vectorRight[0:3], p.vectorBottom[0:3], p.vectorTop[0:3], p.vectorFront[0:3], p.vectorBack[0:3])


class VTKParameters {
public:
  RealType    interval = 0; //! Time interval for file printing
  std::string prefix;       //! Output filename
};
#pragma omp declare mapper(VTKParameters v) \
  map(to: v, v.interval, v.prefix)

class StdOutParameters {
public:
  RealType interval = 0;
};
#pragma omp declare mapper(StdOutParameters s) \
  map(to: s, s.interval)

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

#ifdef ENABLE_PETSC
  PetscInt* sizes[3]; //! Arrays with the sizes of the blocks in each direction.
#else
  int* sizes[3];
#endif
};
// TODO map sizes[3]
#pragma omp declare mapper(ParallelParameters p) \
  map(to: p, p.rank, p.numProcessors[0:3], p.leftNb, p.rightNb, p.bottomNb, p.topNb, p.frontNb, p.backNb, p.indices[0:3], p.localSize[0:3], p.firstCorner[0:3])

class BFStepParameters {
public:
  RealType xRatio = 0;
  RealType yRatio = 0;
};
#pragma omp declare mapper(BFStepParameters b) \
  map(to: b, b.xRatio, b.yRatio)

//@}

/** A class to store and pass around the parameters
 */
class Parameters {
public:
  Parameters();
  ~Parameters();

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
  Meshsize* meshsize;
};
#pragma omp declare mapper(Parameters p) \
  map(to: p.simulation, p.timestep, p.environment, p.flow, p.solver, \
          p.geometry, p.walls, p.vtk, p.parallel, p.stdOut, \
          p.bfStep, p.turbulence, p.meshsize) 
