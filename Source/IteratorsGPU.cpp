#include "StdAfx.hpp"

#include "IteratorsGPU.hpp"

GPUFieldIterator::GPUFieldIterator(FlowField& flowField, const Parameters& parameters, StencilDelegate& stencil, int lowOffset, int highOffset):
  GPUIterator<FlowField>(flowField, parameters),
  stencil_(stencil),
  lowOffset_(lowOffset),
  highOffset_(highOffset) {}

void GPUFieldIterator::iterate() {
  const int cellsX = GPUIterator<FlowField>::flowField_.getCellsX();
  const int cellsY = GPUIterator<FlowField>::flowField_.getCellsY();
  const int cellsZ = GPUIterator<FlowField>::flowField_.getCellsZ();

  Parameters      localParameters = Parameters(GPUIterator<FlowField>::parameters_);
  StencilDelegate localStencil    = stencil_;

  // The index k can be used for the 2D and 3D cases.
  if (GPUIterator<FlowField>::parameters_.geometry.dim == 2) {
    // Loop without lower boundaries. These will be dealt with by the global boundary stencils
    // or by the subdomain boundary GPUIterator.
#pragma omp target data      map(to : lowOffset_, highOffset_, cellsX, cellsY, localParameters, localStencil)
    {
      if (localStencil.getType() == RHS) {
        ScalarField rhsLocal = GPUIterator<FlowField>::flowField_.RHS_;
        VectorField fghLocal = GPUIterator<FlowField>::flowField_.FGH_;
#pragma omp target map(to : fghLocal, fghLocal.data_[0 : fghLocal.size_]) map(tofrom : rhsLocal, rhsLocal.data_[0 : rhsLocal.size_])
#pragma omp teams distribute parallel for
        for (int j = 1 + lowOffset_; j < cellsY - 1 + highOffset_; j++) {
          for (int i = 1 + lowOffset_; i < cellsX - 1 + highOffset_; i++) {
            localStencil.applyRHS(localParameters, rhsLocal, fghLocal, i, j);
          }
        }
        for (int j = 1 + lowOffset_; j < cellsY - 1 + highOffset_; j++) {
          for (int i = 1 + lowOffset_; i < cellsX - 1 + highOffset_; i++) {
            GPUIterator<FlowField>::flowField_.getRHS().getScalar(i, j) = rhsLocal.getScalar(i, j);
          }
        }
      } else if (localStencil.getType() == FGH) {
        VectorField    velocityLocal = GPUIterator<FlowField>::flowField_.velocity_;
        IntScalarField flagsLocal    = GPUIterator<FlowField>::flowField_.flags_;
        VectorField    fghLocal      = GPUIterator<FlowField>::flowField_.FGH_;
#pragma omp target map(to : velocityLocal, velocityLocal.data_[0 : velocityLocal.size_]) map(to : flagsLocal, flagsLocal.data_[0 : flagsLocal.size_]) \
  map(tofrom : fghLocal, fghLocal.data_[0 : fghLocal.size_])
#pragma omp teams distribute parallel for
        for (int j = 1 + lowOffset_; j < cellsY - 1 + highOffset_; j++) {
          for (int i = 1 + lowOffset_; i < cellsX - 1 + highOffset_; i++) {
            localStencil.applyFGH(localParameters, velocityLocal, fghLocal, flagsLocal, i, j);
          }
        }
        for (int j = 1 + lowOffset_; j < cellsY - 1 + highOffset_; j++) {
          for (int i = 1 + lowOffset_; i < cellsX - 1 + highOffset_; i++) {
            GPUIterator<FlowField>::flowField_.getFGH().getVectorElement(i, j, 0) = fghLocal.getVectorElement(i, j, 0);
            GPUIterator<FlowField>::flowField_.getFGH().getVectorElement(i, j, 1) = fghLocal.getVectorElement(i, j, 1);
          }
        }
      } else if (localStencil.getType() == VELOCITY) {
        VectorField    velocityLocal = GPUIterator<FlowField>::flowField_.velocity_;
        VectorField    fghLocal      = GPUIterator<FlowField>::flowField_.FGH_;
        ScalarField    pressureLocal = GPUIterator<FlowField>::flowField_.pressure_;
        IntScalarField flagsLocal    = GPUIterator<FlowField>::flowField_.flags_;
#pragma omp target map(tofrom : velocityLocal, velocityLocal.data_[0 : velocityLocal.size_]) map(to : flagsLocal, flagsLocal.data_[0 : flagsLocal.size_]) \
  map(to : fghLocal, fghLocal.data_[0 : fghLocal.size_]) map(to : pressureLocal, pressureLocal.data_[0 : pressureLocal.size_])
#pragma omp teams distribute parallel for
        for (int j = 1 + lowOffset_; j < cellsY - 1 + highOffset_; j++) {
          for (int i = 1 + lowOffset_; i < cellsX - 1 + highOffset_; i++) {
            localStencil.applyVelocity(localParameters, velocityLocal, fghLocal, pressureLocal, flagsLocal, i, j);
          }
        }
        for (int j = 1 + lowOffset_; j < cellsY - 1 + highOffset_; j++) {
          for (int i = 1 + lowOffset_; i < cellsX - 1 + highOffset_; i++) {
            GPUIterator<FlowField>::flowField_.getVelocity().getVectorElement(i, j, 0) = velocityLocal.getVectorElement(i, j, 0);
            GPUIterator<FlowField>::flowField_.getVelocity().getVectorElement(i, j, 1) = velocityLocal.getVectorElement(i, j, 1);
          }
        }
      } else {
        assert(false);
      }
    }
  }

  if (GPUIterator<FlowField>::parameters_.geometry.dim == 3) {
#pragma omp target data      map(to : lowOffset_, highOffset_, cellsX, cellsY, cellsZ, localParameters, localStencil)
    {
      if (localStencil.getType() == RHS) {
        ScalarField rhsLocal = GPUIterator<FlowField>::flowField_.RHS_;
        VectorField fghLocal = GPUIterator<FlowField>::flowField_.FGH_;
#pragma omp target map(to : fghLocal, fghLocal.data_[0 : fghLocal.size_]) map(tofrom : rhsLocal, rhsLocal.data_[0 : rhsLocal.size_])
#pragma omp teams distribute parallel for
        for (int k = 1 + lowOffset_; k < cellsZ - 1 + highOffset_; k++) {
          for (int j = 1 + lowOffset_; j < cellsY - 1 + highOffset_; j++) {
            for (int i = 1 + lowOffset_; i < cellsX - 1 + highOffset_; i++) {
              localStencil.applyRHS(localParameters, rhsLocal, fghLocal, i, j, k);
            }
          }
        }
        for (int k = 1 + lowOffset_; k < cellsZ - 1 + highOffset_; k++) {
          for (int j = 1 + lowOffset_; j < cellsY - 1 + highOffset_; j++) {
            for (int i = 1 + lowOffset_; i < cellsX - 1 + highOffset_; i++) {
              GPUIterator<FlowField>::flowField_.getRHS().getScalar(i, j, k) = rhsLocal.getScalar(i, j, k);
            }
          }
        }
      } else if (localStencil.getType() == FGH) {
        VectorField    velocityLocal = GPUIterator<FlowField>::flowField_.velocity_;
        IntScalarField flagsLocal    = GPUIterator<FlowField>::flowField_.flags_;
        VectorField    fghLocal      = GPUIterator<FlowField>::flowField_.FGH_;
#pragma omp target map(to : velocityLocal, velocityLocal.data_[0 : velocityLocal.size_]) map(to : flagsLocal, flagsLocal.data_[0 : flagsLocal.size_]) \
  map(tofrom : fghLocal, fghLocal.data_[0 : fghLocal.size_])
#pragma omp teams distribute parallel for
        for (int k = 1 + lowOffset_; k < cellsZ - 1 + highOffset_; k++) {
          for (int j = 1 + lowOffset_; j < cellsY - 1 + highOffset_; j++) {
            for (int i = 1 + lowOffset_; i < cellsX - 1 + highOffset_; i++) {
              localStencil.applyFGH(localParameters, velocityLocal, fghLocal, flagsLocal, i, j, k);
            }
          }
        }
        for (int k = 1 + lowOffset_; k < cellsZ - 1 + highOffset_; k++) {
          for (int j = 1 + lowOffset_; j < cellsY - 1 + highOffset_; j++) {
            for (int i = 1 + lowOffset_; i < cellsX - 1 + highOffset_; i++) {
              GPUIterator<FlowField>::flowField_.getFGH().getVectorElement(i, j, k, 0) = fghLocal.getVectorElement(i, j, k, 0);
              GPUIterator<FlowField>::flowField_.getFGH().getVectorElement(i, j, k, 1) = fghLocal.getVectorElement(i, j, k, 1);
              GPUIterator<FlowField>::flowField_.getFGH().getVectorElement(i, j, k, 2) = fghLocal.getVectorElement(i, j, k, 2);
            }
          }
        }
      } else if (localStencil.getType() == VELOCITY) {
        VectorField    velocityLocal = GPUIterator<FlowField>::flowField_.velocity_;
        VectorField    fghLocal      = GPUIterator<FlowField>::flowField_.FGH_;
        ScalarField    pressureLocal = GPUIterator<FlowField>::flowField_.pressure_;
        IntScalarField flagsLocal    = GPUIterator<FlowField>::flowField_.flags_;
#pragma omp target map(tofrom : velocityLocal, velocityLocal.data_[0 : velocityLocal.size_]) map(to : flagsLocal, flagsLocal.data_[0 : flagsLocal.size_]) \
  map(to : fghLocal, fghLocal.data_[0 : fghLocal.size_]) map(to : pressureLocal, pressureLocal.data_[0 : pressureLocal.size_])
#pragma omp teams distribute parallel for
        for (int k = 1 + lowOffset_; k < cellsZ - 1 + highOffset_; k++) {
          for (int j = 1 + lowOffset_; j < cellsY - 1 + highOffset_; j++) {
            for (int i = 1 + lowOffset_; i < cellsX - 1 + highOffset_; i++) {
              localStencil.applyVelocity(localParameters, velocityLocal, fghLocal, pressureLocal, flagsLocal, i, j, k);
            }
          }
        }
        for (int k = 1 + lowOffset_; k < cellsZ - 1 + highOffset_; k++) {
          for (int j = 1 + lowOffset_; j < cellsY - 1 + highOffset_; j++) {
            for (int i = 1 + lowOffset_; i < cellsX - 1 + highOffset_; i++) {
              GPUIterator<FlowField>::flowField_.getVelocity().getVectorElement(i, j, k, 0) = velocityLocal.getVectorElement(i, j, k, 0);
              GPUIterator<FlowField>::flowField_.getVelocity().getVectorElement(i, j, k, 1) = velocityLocal.getVectorElement(i, j, k, 1);
              GPUIterator<FlowField>::flowField_.getVelocity().getVectorElement(i, j, k, 2) = velocityLocal.getVectorElement(i, j, k, 2);
            }
          }
        }
      } else {
        assert(false);
      }
    }
  }
}

GPUFieldIteratorTurbulent::GPUFieldIteratorTurbulent(TurbulentFlowField& flowField, const Parameters& parameters, StencilDelegate& stencil, int lowOffset, int highOffset):
  GPUIterator<TurbulentFlowField>(flowField, parameters),
  stencil_(stencil),
  lowOffset_(lowOffset),
  highOffset_(highOffset) {}

void GPUFieldIteratorTurbulent::iterate() {
  const int cellsX = GPUIterator<TurbulentFlowField>::flowField_.getCellsX();
  const int cellsY = GPUIterator<TurbulentFlowField>::flowField_.getCellsY();
  const int cellsZ = GPUIterator<TurbulentFlowField>::flowField_.getCellsZ();

  // Parameters           localParameters = GPUIterator<FlowFieldType>::parameters_;
  Parameters      localParameters = Parameters(GPUIterator<TurbulentFlowField>::parameters_);
  StencilDelegate localStencil    = stencil_;

  // The index k can be used for the 2D and 3D cases.
  if (GPUIterator<TurbulentFlowField>::parameters_.geometry.dim == 2) {
// Loop without lower boundaries. These will be dealt with by the global boundary stencils
// or by the subdomain boundary GPUIterator.
#pragma omp target data map(to : lowOffset_, highOffset_, cellsX, cellsY, localParameters, localStencil)
    {
      if (localStencil.getType() == VISCOSITY) {
        IntScalarField flagsLocal     = GPUIterator<TurbulentFlowField>::flowField_.flags_;
        ScalarField    viscosityLocal = GPUIterator<TurbulentFlowField>::flowField_.viscosity_;
        ScalarField    distanceLocal  = GPUIterator<TurbulentFlowField>::flowField_.distance_;
        VectorField    velocityLocal  = GPUIterator<TurbulentFlowField>::flowField_.velocity_;
#pragma omp target map(to : velocityLocal, velocityLocal.data_[0 : velocityLocal.size_]) map(to : flagsLocal, flagsLocal.data_[0 : flagsLocal.size_]) \
  map(to : distanceLocal, distanceLocal.data_[0 : distanceLocal.size_]) map(tofrom : viscosityLocal, viscosityLocal.data_[0 : viscosityLocal.size_])
#pragma omp teams distribute parallel for
        for (int j = 1 + lowOffset_; j < cellsY - 1 + highOffset_; j++) {
          for (int i = 1 + lowOffset_; i < cellsX - 1 + highOffset_; i++) {
            localStencil.applyViscosity(localParameters, velocityLocal, viscosityLocal, distanceLocal, flagsLocal, i, j);
          }
        }
        for (int j = 1 + lowOffset_; j < cellsY - 1 + highOffset_; j++) {
          for (int i = 1 + lowOffset_; i < cellsX - 1 + highOffset_; i++) {
            GPUIterator<TurbulentFlowField>::flowField_.getViscosity().getScalar(i, j) = viscosityLocal.getScalar(i, j);
          }
        }
      } else if (stencil_.getType() == TURBFGH) {
        VectorField    velocityLocal  = GPUIterator<TurbulentFlowField>::flowField_.velocity_;
        IntScalarField flagsLocal     = GPUIterator<TurbulentFlowField>::flowField_.flags_;
        VectorField    fghLocal       = GPUIterator<TurbulentFlowField>::flowField_.FGH_;
        ScalarField    viscosityLocal = GPUIterator<TurbulentFlowField>::flowField_.viscosity_;
#pragma omp target map(to : velocityLocal, velocityLocal.data_[0 : velocityLocal.size_]) map(to : flagsLocal, flagsLocal.data_[0 : flagsLocal.size_]) \
  map(tofrom : fghLocal, fghLocal.data_[0 : fghLocal.size_]) map(to : viscosityLocal, viscosityLocal.data_[0 : viscosityLocal.size_])
#pragma omp teams distribute parallel for
        for (int j = 1 + lowOffset_; j < cellsY - 1 + highOffset_; j++) {
          for (int i = 1 + lowOffset_; i < cellsX - 1 + highOffset_; i++) {
            localStencil.applyTurbulentFGH(localParameters, fghLocal, viscosityLocal, velocityLocal, flagsLocal, i, j);
          }
        }
        for (int j = 1 + lowOffset_; j < cellsY - 1 + highOffset_; j++) {
          for (int i = 1 + lowOffset_; i < cellsX - 1 + highOffset_; i++) {
            GPUIterator<TurbulentFlowField>::flowField_.getFGH().getVectorElement(i, j, 0) = fghLocal.getVectorElement(i, j, 0);
            GPUIterator<TurbulentFlowField>::flowField_.getFGH().getVectorElement(i, j, 1) = fghLocal.getVectorElement(i, j, 1);
            GPUIterator<TurbulentFlowField>::flowField_.getFGH().getVectorElement(i, j, 2) = fghLocal.getVectorElement(i, j, 2);
          }
        }
      } else {
        assert(false);
      }
    }
  }

  if (GPUIterator<TurbulentFlowField>::parameters_.geometry.dim == 3) {
#pragma omp target data map(to : lowOffset_, highOffset_, cellsX, cellsY, cellsZ, localParameters, localStencil)
    {
      if (localStencil.getType() == VISCOSITY) {
        IntScalarField flagsLocal     = GPUIterator<TurbulentFlowField>::flowField_.flags_;
        ScalarField    viscosityLocal = GPUIterator<TurbulentFlowField>::flowField_.viscosity_;
        ScalarField    distanceLocal  = GPUIterator<TurbulentFlowField>::flowField_.distance_;
        VectorField    velocityLocal  = GPUIterator<TurbulentFlowField>::flowField_.velocity_;
#pragma omp target map(to : velocityLocal, velocityLocal.data_[0 : velocityLocal.size_]) map(to : flagsLocal, flagsLocal.data_[0 : flagsLocal.size_]) \
  map(to : distanceLocal, distanceLocal.data_[0 : distanceLocal.size_]) map(tofrom : viscosityLocal, viscosityLocal.data_[0 : viscosityLocal.size_])
#pragma omp teams distribute parallel for
        for (int k = 1 + lowOffset_; k < cellsZ - 1 + highOffset_; k++) {
          for (int j = 1 + lowOffset_; j < cellsY - 1 + highOffset_; j++) {
            for (int i = 1 + lowOffset_; i < cellsX - 1 + highOffset_; i++) {
              localStencil.applyViscosity(localParameters, velocityLocal, viscosityLocal, distanceLocal, flagsLocal, i, j, k);
            }
          }
        }
        for (int k = 1 + lowOffset_; k < cellsZ - 1 + highOffset_; k++) {
          for (int j = 1 + lowOffset_; j < cellsY - 1 + highOffset_; j++) {
            for (int i = 1 + lowOffset_; i < cellsX - 1 + highOffset_; i++) {
              GPUIterator<TurbulentFlowField>::flowField_.getViscosity().getScalar(i, j, k) = viscosityLocal.getScalar(i, j, k);
            }
          }
        }
      } else if (localStencil.getType() == TURBFGH) {
        VectorField    velocityLocal  = GPUIterator<TurbulentFlowField>::flowField_.velocity_;
        IntScalarField flagsLocal     = GPUIterator<TurbulentFlowField>::flowField_.flags_;
        VectorField    fghLocal       = GPUIterator<TurbulentFlowField>::flowField_.FGH_;
        ScalarField    viscosityLocal = GPUIterator<TurbulentFlowField>::flowField_.viscosity_;
#pragma omp target map(to : velocityLocal, velocityLocal.data_[0 : velocityLocal.size_]) map(to : flagsLocal, flagsLocal.data_[0 : flagsLocal.size_]) \
  map(tofrom : fghLocal, fghLocal.data_[0 : fghLocal.size_]) map(to : viscosityLocal, viscosityLocal.data_[0 : viscosityLocal.size_])
#pragma omp teams distribute parallel for
        for (int k = 1 + lowOffset_; k < cellsZ - 1 + highOffset_; k++) {
          for (int j = 1 + lowOffset_; j < cellsY - 1 + highOffset_; j++) {
            for (int i = 1 + lowOffset_; i < cellsX - 1 + highOffset_; i++) {
              localStencil.applyTurbulentFGH(localParameters, fghLocal, viscosityLocal, velocityLocal, flagsLocal, i, j, k);
            }
          }
        }
        for (int k = 1 + lowOffset_; k < cellsZ - 1 + highOffset_; k++) {
          for (int j = 1 + lowOffset_; j < cellsY - 1 + highOffset_; j++) {
            for (int i = 1 + lowOffset_; i < cellsX - 1 + highOffset_; i++) {
              GPUIterator<TurbulentFlowField>::flowField_.getFGH().getVectorElement(i, j, k, 0) = fghLocal.getVectorElement(i, j, k, 0);
              GPUIterator<TurbulentFlowField>::flowField_.getFGH().getVectorElement(i, j, k, 1) = fghLocal.getVectorElement(i, j, k, 1);
              GPUIterator<TurbulentFlowField>::flowField_.getFGH().getVectorElement(i, j, k, 2) = fghLocal.getVectorElement(i, j, k, 2);
            }
          }
        }
      } else {
        assert(false);
      }
    }
  }
}
