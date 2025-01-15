#pragma once

#include "Definitions.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"

namespace Solvers {

  class SORSolver {
  public:
    SORSolver()  = default;
    ~SORSolver() = default;

#pragma omp declare target
    void solve(FlowField& flowField, const Parameters& parameters);
#pragma omp end declare target

    inline void reInitMatrix() {};
  };

} // namespace Solvers
