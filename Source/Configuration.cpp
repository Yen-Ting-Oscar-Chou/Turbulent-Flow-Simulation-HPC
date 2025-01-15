#include "StdAfx.hpp"

#include "Configuration.hpp"

#include <tinyxml2.h>

void readFloatMandatory(RealType& storage, tinyxml2::XMLElement* node, const char* tag) {
  double value; // Use to be able to select precision
  if (node->QueryDoubleAttribute(tag, &value) != tinyxml2::XML_SUCCESS) {
    throw std::runtime_error("Error while reading mandatory argument");
  } else {
    storage = static_cast<RealType>(value);
  }
}

void readFloatOptional(RealType& storage, tinyxml2::XMLElement* node, const char* tag, RealType defaultValue = 0) {
  double value; // Use to be able to select precision
  int    result = node->QueryDoubleAttribute(tag, &value);
  if (result == tinyxml2::XML_NO_ATTRIBUTE) {
    storage = defaultValue;
  } else if (result == tinyxml2::XML_WRONG_ATTRIBUTE_TYPE) {
    throw std::runtime_error("Error while reading optional argument");
  } else {
    storage = static_cast<RealType>(value);
  }
}

void readIntMandatory(int& storage, tinyxml2::XMLElement* node, const char* tag) {
  int value;
  if (node->QueryIntAttribute(tag, &value) != tinyxml2::XML_SUCCESS) {
    throw std::runtime_error("Error while reading mandatory argument");
  } else {
    storage = value;
  }
}

void readIntOptional(int& storage, tinyxml2::XMLElement* node, const char* tag, int defaultValue = 0) {
  int result = node->QueryIntAttribute(tag, &storage);
  if (result == tinyxml2::XML_NO_ATTRIBUTE) {
    storage = defaultValue;
  } else if (result == tinyxml2::XML_WRONG_ATTRIBUTE_TYPE) {
    throw std::runtime_error("Error while reading optional argument");
  }
}

void readBoolMandatory(bool& storage, tinyxml2::XMLElement* node, const char* tag) {
  bool value;
  if (node->QueryBoolAttribute(tag, &value) != tinyxml2::XML_SUCCESS) {
    throw std::runtime_error("Error while reading mandatory argument");
  } else {
    storage = value;
  }
}

void readBoolOptional(bool& storage, tinyxml2::XMLElement* node, const char* tag, bool defaultValue = false) {
  int result = node->QueryBoolAttribute(tag, &storage);
  if (result == tinyxml2::XML_NO_ATTRIBUTE) {
    storage = defaultValue;
  } else if (result == tinyxml2::XML_WRONG_ATTRIBUTE_TYPE) {
    throw std::runtime_error("Error while reading optional argument");
  }
}

void readStringMandatory(std::string& storage, tinyxml2::XMLElement* node) {
  const char* myText = node->GetText();
  if (myText == NULL) {
    const std::string nodename = node->Name();
    spdlog::error("No string specified for this node: {}", nodename);
    throw std::runtime_error("Error while reading mandatory string");
  } else {
    storage = node->GetText();
    if (!storage.compare("")) {
      throw std::runtime_error("Missing mandatory string!");
    }
  }
}

void readWall(tinyxml2::XMLElement* wall, RealType* vector, RealType& scalar) {
  tinyxml2::XMLElement* quantity = wall->FirstChildElement("vector");
  if (quantity != NULL) {
    readFloatOptional(vector[0], quantity, "x");
    readFloatOptional(vector[1], quantity, "y");
    readFloatOptional(vector[2], quantity, "z");
  }
  quantity = wall->FirstChildElement("scalar");
  if (quantity != NULL) {
    readFloatOptional(scalar, quantity, "value");
  }
}

void broadcastString(std::string& target, const MPI_Comm& communicator, int root = 0) {
  int stringSize = 0, rank = -1;
  MPI_Comm_rank(communicator, &rank);
  if (rank == root) {
    stringSize = static_cast<int>(target.size());
  }
  MPI_Bcast(&stringSize, 1, MPI_INT, 0, communicator);
  char* name = new char[stringSize + 1]; // One more for the null character
  if (rank == root) {
    target.copy(name, stringSize, 0);
  }
  name[stringSize] = '\0';
  MPI_Bcast(name, stringSize + 1, MPI_CHAR, 0, communicator);
  if (rank != root) {
    target = name;
  }
  delete[] name;
  name = NULL;
}

Configuration::Configuration() { filename_ = ""; }

Configuration::Configuration(const std::string& filename) { filename_ = filename; }

void Configuration::setFileName(const std::string& filename) { filename_ = filename; }

Parameters Configuration::loadParameters(const MPI_Comm& communicator) {
  tinyxml2::XMLDocument confFile;
  tinyxml2::XMLElement* node;
  tinyxml2::XMLElement* subNode;
  std::string           deltaMixLen = "";
  GeometricParameters   geometricParameters;
  ParallelParameters    parallelParameters;

  int rank = -1;
  MPI_Comm_rank(communicator, &rank);

  // We only read on rank 0; afterwards, all configuration parameters are broadcasted to all processes.
  // So, if you add new parameters in the configuration, make sure to broadcast them to the other processes!
  if (rank == 0) {
    // Parse the configuration file and check validity
    confFile.LoadFile(filename_.c_str());
    if (confFile.FirstChildElement() == NULL) {
      throw std::runtime_error("Error parsing the configuration file");
    }

    //--------------------------------------------------
    // Load geometric parameters
    //--------------------------------------------------
    node = confFile.FirstChildElement()->FirstChildElement("geometry");

    if (node == NULL) {
      throw std::runtime_error("Error loading geometry properties");
    }

    readIntMandatory(geometricParameters.sizeX, node, "sizeX");
    readIntMandatory(geometricParameters.sizeY, node, "sizeY");
    readIntOptional(geometricParameters.sizeZ, node, "sizeZ");

    if (geometricParameters.sizeX < 2 || geometricParameters.sizeY < 2 || geometricParameters.sizeZ < 0) {
      throw std::runtime_error("Invalid size specified in configuration file");
    }

    geometricParameters.dim = 0;
    if (node->QueryIntAttribute("dim", &(geometricParameters.dim)) != tinyxml2::XML_WRONG_ATTRIBUTE_TYPE) {
      if (geometricParameters.dim == 0) {
        if (geometricParameters.sizeZ == 0) {
          geometricParameters.sizeZ = 1;
          geometricParameters.dim   = 2;
        } else {
          geometricParameters.dim = 3;
        }
      }
    }

    if (geometricParameters.dim == 3 && geometricParameters.sizeZ == 1) {
      throw std::runtime_error("Inconsistent data: 3D geometry specified with Z size zero");
    }

    if (geometricParameters.dim == 2 && geometricParameters.sizeZ != 1) {
      geometricParameters.sizeZ = 1;
    }

    // Determine the sizes of the cells
    readFloatMandatory(geometricParameters.lengthX, node, "lengthX");
    readFloatMandatory(geometricParameters.lengthY, node, "lengthY");
    readFloatMandatory(geometricParameters.lengthZ, node, "lengthZ");
    // Read geometry->meshsize parameters
    std::string meshsizeType = "";
    subNode                  = node->FirstChildElement("mesh");
    readStringMandatory(meshsizeType, subNode);
    if (meshsizeType == "uniform") {
      geometricParameters.meshsizeType = Uniform;
    } else if (meshsizeType == "stretched") {
      geometricParameters.meshsizeType = TanhStretching;
      bool buffer                      = false;
      readBoolMandatory(buffer, node, "stretchX");
      geometricParameters.stretchX = static_cast<int>(buffer);
      readBoolMandatory(buffer, node, "stretchY");
      geometricParameters.stretchY = static_cast<int>(buffer);
      if (geometricParameters.dim == 3) {
        readBoolMandatory(buffer, node, "stretchZ");
        geometricParameters.stretchZ = static_cast<int>(buffer);
      } else {
        geometricParameters.stretchZ = false;
      }
    } else {
      throw std::runtime_error("Unknown 'mesh'!");
    }

    // Now, the size of the elements should be set

    dim_ = geometricParameters.dim;

    //--------------------------------------------------
    // Parallel parameters
    //--------------------------------------------------
    node = confFile.FirstChildElement()->FirstChildElement("parallel");

    if (node == NULL) {
      throw std::runtime_error("Error loading parallel parameters");
    }

    readIntOptional(parallelParameters.numProcessors[0], node, "numProcessorsX", 1);
    readIntOptional(parallelParameters.numProcessors[1], node, "numProcessorsY", 1);
    readIntOptional(parallelParameters.numProcessors[2], node, "numProcessorsZ", 1);

    // Start neighbors on null in case that no parallel configuration is used later.
    parallelParameters.leftNb   = MPI_PROC_NULL;
    parallelParameters.rightNb  = MPI_PROC_NULL;
    parallelParameters.bottomNb = MPI_PROC_NULL;
    parallelParameters.topNb    = MPI_PROC_NULL;
    parallelParameters.frontNb  = MPI_PROC_NULL;
    parallelParameters.backNb   = MPI_PROC_NULL;

    // Yet more parameters initialized in case that no parallel configuration is applied
    parallelParameters.localSize[0] = geometricParameters.sizeX;
    parallelParameters.localSize[1] = geometricParameters.sizeY;
    parallelParameters.localSize[2] = geometricParameters.sizeZ;

    parallelParameters.firstCorner[0] = 0;
    parallelParameters.firstCorner[1] = 0;
    parallelParameters.firstCorner[2] = 0;

    // VTK output is named after the rank, so we define it here, again, in case that it's not
    // initialized anywhere else.
    parallelParameters.rank = rank;
  }

  MPI_Bcast(&(geometricParameters.sizeX), 1, MPI_INT, 0, communicator);
  MPI_Bcast(&(geometricParameters.sizeY), 1, MPI_INT, 0, communicator);
  MPI_Bcast(&(geometricParameters.sizeZ), 1, MPI_INT, 0, communicator);

  MPI_Bcast(&(geometricParameters.dim), 1, MPI_INT, 0, communicator);

  MPI_Bcast(&(geometricParameters.meshsizeType), 1, MPI_INT, 0, communicator);
  MPI_Bcast(&(geometricParameters.stretchX), 1, MPI_INT, 0, communicator);
  MPI_Bcast(&(geometricParameters.stretchY), 1, MPI_INT, 0, communicator);
  MPI_Bcast(&(geometricParameters.stretchZ), 1, MPI_INT, 0, communicator);
  MPI_Bcast(&(geometricParameters.lengthX), 1, MY_MPI_FLOAT, 0, communicator);
  MPI_Bcast(&(geometricParameters.lengthY), 1, MY_MPI_FLOAT, 0, communicator);
  MPI_Bcast(&(geometricParameters.lengthZ), 1, MY_MPI_FLOAT, 0, communicator);

  MPI_Bcast(parallelParameters.numProcessors, 3, MPI_INT, 0, communicator);

  Parameters parameters(geometricParameters, parallelParameters);

  if (rank == 0) {
    //--------------------------------------------------
    // Timestep parameters
    //--------------------------------------------------
    node = confFile.FirstChildElement()->FirstChildElement("timestep");

    if (node == NULL) {
      throw std::runtime_error("Error loading timestep parameters");
    }

    readFloatOptional(parameters.timestep.dt, node, "dt", 1);
    readFloatOptional(parameters.timestep.tau, node, "tau", 0.5);

    //--------------------------------------------------
    // Flow parameters
    //--------------------------------------------------
    node = confFile.FirstChildElement()->FirstChildElement("flow");

    if (node == NULL) {
      throw std::runtime_error("Error loading flow parameters");
    }

    readFloatMandatory(parameters.flow.Re, node, "Re");

    //--------------------------------------------------
    // Solver parameters
    //--------------------------------------------------
    node = confFile.FirstChildElement()->FirstChildElement("solver");

    if (node == NULL) {
      throw std::runtime_error("Error loading solver parameters");
    }

    readFloatMandatory(parameters.solver.gamma, node, "gamma");
    readIntOptional(parameters.solver.maxIterations, node, "maxIterations");

    //--------------------------------------------------
    // Environmental parameters
    //--------------------------------------------------
    node = confFile.FirstChildElement()->FirstChildElement("environment");

    if (node == NULL) {
      throw std::runtime_error("Error loading environmental parameters");
    }

    readFloatOptional(parameters.environment.gx, node, "gx");
    readFloatOptional(parameters.environment.gy, node, "gy");
    readFloatOptional(parameters.environment.gz, node, "gz");

    //--------------------------------------------------
    // Simulation parameters
    //--------------------------------------------------
    node = confFile.FirstChildElement()->FirstChildElement("simulation");

    if (node == NULL) {
      throw std::runtime_error("Error loading simulation parameters");
    }

    readFloatMandatory(parameters.simulation.finalTime, node, "finalTime");

    subNode = node->FirstChildElement("type");
    if (subNode != NULL) {
      readStringMandatory(parameters.simulation.type, subNode);
    } else {
      throw std::runtime_error("Missing type in simulation parameters");
    }

    subNode = node->FirstChildElement("scenario");
    if (subNode != NULL) {
      std::string scenario;
      readStringMandatory(scenario, subNode);
      if (scenario == "cavity") {
        parameters.simulation.scenario = CAVITY;
      } else if (scenario == "channel") {
        parameters.simulation.scenario = CHANNEL;
      } else {
        throw std::runtime_error("Unkown scenario");
      }
    } else {
      throw std::runtime_error("Missing scenario in simulation parameters");
    }

    subNode = node->FirstChildElement("velocityProfile");
    if (subNode != NULL) {
      readStringMandatory(parameters.simulation.velocityProfile, subNode);
      if (parameters.simulation.velocityProfile != "parabolic" && parameters.simulation.velocityProfile != "uniform") {
        throw std::runtime_error("Unsupported velocity profile in simulation parameters");
      }
    }

    //--------------------------------------------------
    // VTK parameters
    //--------------------------------------------------
    node = confFile.FirstChildElement()->FirstChildElement("vtk");

    if (node == NULL) {
      throw std::runtime_error("Error loading VTK parameters");
    }

    readFloatOptional(parameters.vtk.interval, node, "interval");
    readStringMandatory(parameters.vtk.prefix, node);

    //--------------------------------------------------
    // StdOut parameters
    //--------------------------------------------------
    node = confFile.FirstChildElement()->FirstChildElement("stdOut");

    if (node == NULL) {
      throw std::runtime_error("Error loading StdOut parameters");
    }

    // If no value given, print every step
    readFloatOptional(parameters.stdOut.interval, node, "interval", 1);

    //--------------------------------------------------
    // Walls
    //--------------------------------------------------
    node = confFile.FirstChildElement()->FirstChildElement("walls");

    if (node == NULL) {
      throw std::runtime_error("Error loading wall parameters");
    }

    tinyxml2::XMLElement* wall;
    wall = node->FirstChildElement("left");
    if (wall != NULL) {
      readWall(wall, parameters.walls.vectorLeft, parameters.walls.scalarLeft);
    }

    wall = node->FirstChildElement("right");
    if (wall != NULL) {
      readWall(wall, parameters.walls.vectorRight, parameters.walls.scalarRight);
    }

    wall = node->FirstChildElement("bottom");
    if (wall != NULL) {
      readWall(wall, parameters.walls.vectorBottom, parameters.walls.scalarBottom);
    }

    wall = node->FirstChildElement("top");
    if (wall != NULL) {
      readWall(wall, parameters.walls.vectorTop, parameters.walls.scalarTop);
    }

    wall = node->FirstChildElement("front");
    if (wall != NULL) {
      readWall(wall, parameters.walls.vectorFront, parameters.walls.scalarFront);
    }

    wall = node->FirstChildElement("back");
    if (wall != NULL) {
      readWall(wall, parameters.walls.vectorBack, parameters.walls.scalarBack);
    }

    // Set the scalar values to zero;
    parameters.walls.scalarRight  = 0.0;
    parameters.walls.scalarBottom = 0.0;
    parameters.walls.scalarTop    = 0.0;
    parameters.walls.scalarFront  = 0.0;
    parameters.walls.scalarBack   = 0.0;

    //--------------------------------------------------
    // Backward facing step
    //--------------------------------------------------
    parameters.bfStep.xRatio = -1.0;
    parameters.bfStep.yRatio = -1.0;
    node                     = confFile.FirstChildElement()->FirstChildElement("backwardFacingStep");
    if (node != NULL) {
      readFloatMandatory(parameters.bfStep.xRatio, node, "xRatio");
      readFloatMandatory(parameters.bfStep.yRatio, node, "yRatio");
    }

    //------------------------------------------------------
    // TODO WS2: Turbulence
    //------------------------------------------------------
    if (parameters.simulation.type == "turbulence") {
      node = confFile.FirstChildElement()->FirstChildElement("deltaMixLen");
      readStringMandatory(deltaMixLen, node);
      if (deltaMixLen == "turbulence") {
        parameters.turbulence.deltaMixLen = TURBULENT;
      } else if (deltaMixLen == "laminar") {
        parameters.turbulence.deltaMixLen = LAMINAR;
      } else if (deltaMixLen == "zero") {
        parameters.turbulence.deltaMixLen = ZERO;
      } else {
        throw std::runtime_error("Error loading delta for mixing lengths");
      }
    }
  }

  // Broadcasting of the values
  MPI_Bcast(&(parameters.timestep.dt), 1, MY_MPI_FLOAT, 0, communicator);
  MPI_Bcast(&(parameters.timestep.tau), 1, MY_MPI_FLOAT, 0, communicator);

  MPI_Bcast(&(parameters.flow.Re), 1, MY_MPI_FLOAT, 0, communicator);

  MPI_Bcast(&(parameters.solver.gamma), 1, MY_MPI_FLOAT, 0, communicator);
  MPI_Bcast(&(parameters.solver.maxIterations), 1, MY_MPI_FLOAT, 0, communicator);

  MPI_Bcast(&(parameters.environment.gx), 1, MY_MPI_FLOAT, 0, communicator);
  MPI_Bcast(&(parameters.environment.gy), 1, MY_MPI_FLOAT, 0, communicator);
  MPI_Bcast(&(parameters.environment.gz), 1, MY_MPI_FLOAT, 0, communicator);

  MPI_Bcast(&(parameters.simulation.finalTime), 1, MY_MPI_FLOAT, 0, communicator);

  MPI_Bcast(&(parameters.vtk.interval), 1, MY_MPI_FLOAT, 0, communicator);
  MPI_Bcast(&(parameters.stdOut.interval), 1, MPI_INT, 0, communicator);

  broadcastString(parameters.vtk.prefix, communicator);
  broadcastString(parameters.simulation.type, communicator);
  MPI_Bcast(&(parameters.simulation.scenario), 1, MPI_INT, 0, communicator);

  MPI_Bcast(&(parameters.bfStep.xRatio), 1, MY_MPI_FLOAT, 0, communicator);
  MPI_Bcast(&(parameters.bfStep.yRatio), 1, MY_MPI_FLOAT, 0, communicator);

  MPI_Bcast(&(parameters.walls.scalarLeft), 1, MY_MPI_FLOAT, 0, communicator);
  MPI_Bcast(&(parameters.walls.scalarRight), 1, MY_MPI_FLOAT, 0, communicator);
  MPI_Bcast(&(parameters.walls.scalarBottom), 1, MY_MPI_FLOAT, 0, communicator);
  MPI_Bcast(&(parameters.walls.scalarTop), 1, MY_MPI_FLOAT, 0, communicator);
  MPI_Bcast(&(parameters.walls.scalarFront), 1, MY_MPI_FLOAT, 0, communicator);
  MPI_Bcast(&(parameters.walls.scalarBack), 1, MY_MPI_FLOAT, 0, communicator);

  MPI_Bcast(parameters.walls.vectorLeft, 3, MY_MPI_FLOAT, 0, communicator);
  MPI_Bcast(parameters.walls.vectorRight, 3, MY_MPI_FLOAT, 0, communicator);
  MPI_Bcast(parameters.walls.vectorBottom, 3, MY_MPI_FLOAT, 0, communicator);
  MPI_Bcast(parameters.walls.vectorTop, 3, MY_MPI_FLOAT, 0, communicator);
  MPI_Bcast(parameters.walls.vectorFront, 3, MY_MPI_FLOAT, 0, communicator);
  MPI_Bcast(parameters.walls.vectorBack, 3, MY_MPI_FLOAT, 0, communicator);

  // TODO WS2: broadcast turbulence parameters
  if (parameters.simulation.type == "turbulence") {
    broadcastString(parameters.simulation.velocityProfile, communicator);
    broadcastString(deltaMixLen, communicator);
    if (rank != 0) {
      if (deltaMixLen == "turbulence") {
        parameters.turbulence.deltaMixLen = TURBULENT;
      } else if (deltaMixLen == "laminar") {
        parameters.turbulence.deltaMixLen = LAMINAR;
      } else if (deltaMixLen == "zero") {
        parameters.turbulence.deltaMixLen = ZERO;
      } else {
        throw std::runtime_error("Error loading delta for mixing lengths");
      }
    }
  }

  return parameters;
}
