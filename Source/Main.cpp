#include "StdAfx.hpp"

#include <fenv.h>

#include "Clock.hpp"
#include "Configuration.hpp"
#include "ParallelManagers/PetscParallelConfiguration.hpp"
#include "Simulation.hpp"
#include "TurbulentSimulation.hpp"

int main(int argc, char* argv[]) {
  spdlog::set_level(spdlog::level::info);

  // Parallelisation related. Initialise and identify.
  // ---------------------------------------------------
  int rank  = -1; // This is the processor's identifier
  int nproc = 0;  // Number of processors in the group
#ifdef ENABLE_PETSC
  PetscInitialize(&argc, &argv, "petsc_commandline_arg", PETSC_NULLPTR);
#else
  MPI_Init(&argc, &argv);
#endif

  MPI_Comm_size(PETSC_COMM_WORLD, &nproc);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  spdlog::info("Rank: {}, Nproc: {}", rank, nproc);
  //----------------------------------------------------


#ifdef ENABLE_SINGLE_PRECISION
  spdlog::warn(
    "Using single floating-point precision for benchmarking only; make sure to switch to double precision "
    "for production runs."
  );
#else
  spdlog::info("Using double floating-point precision");
#endif

#ifndef NDEBUG
  spdlog::warn("Running in Debug mode; make sure to switch to Release mode for production/benchmark runs.");
#ifdef __linux__
  feenableexcept(FE_DIVBYZERO | FE_INVALID);
#endif
#else
  spdlog::info("Running in Release mode");
#endif

  if (!argv[1]) {
    spdlog::error("You need to pass a configuration file: mpirun -np 1 ./NS-EOF ExampleCases/Cavity2D.xml.");
    throw std::runtime_error("Argument parsing error");
  }

  // Read configuration and store information in parameters object
  Configuration                                configuration(argv[1]);
  Parameters                                   parameters = configuration.loadParameters();
  ParallelManagers::PetscParallelConfiguration parallelConfiguration(parameters);
  FlowField*                                   flowField  = NULL;
  Simulation*                                  simulation = NULL;

  spdlog::debug(
    "Processor {} with index {}, {}, {} is computing the size of its subdomain and obtains {}, {} and {}.",
    parameters.parallel.rank,
    parameters.parallel.indices[0],
    parameters.parallel.indices[1],
    parameters.parallel.indices[2],
    parameters.parallel.localSize[0],
    parameters.parallel.localSize[1],
    parameters.parallel.localSize[2]
  );
  spdlog::debug("Left neighbour: {}, right neighbour: {}", parameters.parallel.leftNb, parameters.parallel.rightNb);
  spdlog::debug("Top neighbour: {}, bottom neighbour: {}", parameters.parallel.topNb, parameters.parallel.bottomNb);
  spdlog::debug("Front neighbour: {}, back neighbour: {}", parameters.parallel.frontNb, parameters.parallel.backNb);
  spdlog::debug("Min. meshsizes: {}, {}, {}", parameters.meshsize->getDxMin(), parameters.meshsize->getDyMin(), parameters.meshsize->getDzMin());

  // Initialise simulation
  if (parameters.simulation.type == "turbulence") {
    if (rank == 0) {
      spdlog::info("Start turbulence simulation in {}D", parameters.geometry.dim);
    }
    TurbulentFlowField* turbulentFlowField = new TurbulentFlowField(parameters);
    if (turbulentFlowField == NULL) {
      throw std::runtime_error("turbulentFlowField == NULL!");
    }
    TurbulentSimulation* turbulentSimulation = new TurbulentSimulation(&parameters, turbulentFlowField);
    flowField                                = turbulentFlowField;
    simulation                               = turbulentSimulation;
  } else if (parameters.simulation.type == "dns") {
    if (rank == 0) {
      spdlog::info("Start DNS simulation in {}D", parameters.geometry.dim);
    }
    flowField = new FlowField(parameters);
    if (flowField == NULL) {
      throw std::runtime_error("flowField == NULL!");
    }
    simulation = new Simulation(&parameters, flowField);
  } else {
    throw std::runtime_error("Unknown simulation type! Currently supported: dns, turbulence");
  }

  // Call initialisation of simulation (initialise flow field)
  if (simulation == NULL) {
    throw std::runtime_error("simulation == NULL!");
  }
  simulation->initializeFlowField();

  RealType time       = 0.0;
  RealType timeVtk    = parameters.vtk.interval;
  RealType timeStdOut = parameters.stdOut.interval;
  int      timeSteps  = 0;

  // Plot initial state
  simulation->plotVTK(timeSteps, time);

  Clock clock;
  int   hostDevice   = omp_get_initial_device();
  int   targetDevice = omp_get_default_device();

  // Time loop
  if (parameters.simulation.type == "dns") {
    auto ptrs = Simulation::mapToGPU(hostDevice, targetDevice, *simulation);
    while (time < parameters.simulation.finalTime) {
      simulation->solveTimestep(hostDevice, targetDevice, ptrs.flowFieldGPUptrs_);
      Parameters::mapToCPU(hostDevice, targetDevice, *simulation->parameters_, ptrs.parameterPtrs_);
      timeSteps++;
      time += parameters.timestep.dt;

      if ((rank == 0) && (timeStdOut <= time)) {
        spdlog::info("Current time: {}\tTimestep: {}", time, parameters.timestep.dt);
        timeStdOut += parameters.stdOut.interval;
      }

      if (timeVtk <= time) {
        FlowField::mapToCPUVTK(hostDevice, targetDevice, *simulation->flowField_, ptrs.flowFieldGPUptrs_);
        simulation->plotVTK(timeSteps, time);
        timeVtk += parameters.vtk.interval;
      }
    }
    Simulation::mapToCPUAndFree(hostDevice, targetDevice, *simulation, ptrs);
  } else if (parameters.simulation.type == "turbulence") {
    TurbulentSimulation* turbulentSimulation = static_cast<TurbulentSimulation*>(simulation);
    auto                 ptrs                = TurbulentSimulation::mapToGPU(hostDevice, targetDevice, *turbulentSimulation);
    while (time < parameters.simulation.finalTime) {
      turbulentSimulation->solveTimestep(hostDevice, targetDevice, ptrs.flowFieldGPUptrs_);
      Parameters::mapToCPU(hostDevice, targetDevice, *simulation->parameters_, ptrs.parameterPtrs_);
      timeSteps++;
      time += parameters.timestep.dt;

      if ((rank == 0) && (timeStdOut <= time)) {
        spdlog::info("Current time: {}\tTimestep: {}", time, parameters.timestep.dt);
        timeStdOut += parameters.stdOut.interval;
      }

      if (timeVtk <= time) {
        TurbulentFlowField::mapToCPUVTK(hostDevice, targetDevice, *turbulentSimulation->turbulentField_, ptrs.flowFieldGPUptrs_);
        turbulentSimulation->plotVTK(timeSteps, time);
        timeVtk += parameters.vtk.interval;
      }
    }
    TurbulentSimulation::mapToCPUAndFree(hostDevice, targetDevice, *turbulentSimulation, ptrs);
  }
  spdlog::info("Finished simulation with a duration of {}ns", clock.getTime());

  // Plot final solution
  simulation->plotVTK(timeSteps, time);

  delete simulation;
  simulation = NULL;

  delete flowField;
  flowField = NULL;

#ifdef ENABLE_PETSC
  PetscFinalize();
#else
  MPI_Finalize();
#endif

  return EXIT_SUCCESS;
}
