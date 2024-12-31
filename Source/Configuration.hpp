#pragma once

#include "Definitions.hpp"
#include "Parameters.hpp"
#include "MixingLengths.hpp"

class Configuration {
private:
  int         dim_;
  std::string filename_;

public:
  Configuration();
  Configuration(const std::string& filename);
  ~Configuration() = default;

  void setFileName(const std::string& filename);
  Parameters loadParameters(const MPI_Comm& communicator = PETSC_COMM_WORLD);
};
