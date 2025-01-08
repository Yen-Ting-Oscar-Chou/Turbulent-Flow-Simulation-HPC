#include "StdAfx.hpp"

#include "IteratorsGPU.hpp"

template <class FlowFieldType>
GPUFieldIterator<FlowFieldType>::GPUFieldIterator(FlowFieldType& flowField, const Parameters& parameters, FieldStencilDelegate& stencil, int lowOffset, int highOffset):
  GPUIterator<FlowFieldType>(flowField, parameters),
  stencil_(stencil),
  lowOffset_(lowOffset),
  highOffset_(highOffset) {}

template <class FlowFieldType>
void GPUFieldIterator<FlowFieldType>::iterate() {
  const int cellsX = GPUIterator<FlowFieldType>::flowField_.getCellsX();
  const int cellsY = GPUIterator<FlowFieldType>::flowField_.getCellsY();
  const int cellsZ = GPUIterator<FlowFieldType>::flowField_.getCellsZ();

  //Parameters           localParameters = GPUIterator<FlowFieldType>::parameters_;
  Parameters           localParameters = Parameters(GPUIterator<FlowFieldType>::parameters_);
  FieldStencilDelegate localStencil    = stencil_;
  VectorField    fghLocal      = GPUIterator<FlowFieldType>::flowField_.FGH_;
  ScalarField rhsLocal = GPUIterator<FlowFieldType>::flowField_.RHS_;

  // The index k can be used for the 2D and 3D cases.
  if (GPUIterator<FlowFieldType>::parameters_.geometry.dim == 2) {
    // Loop without lower boundaries. These will be dealt with by the global boundary stencils
    // or by the subdomain boundary GPUIterator.

    #pragma omp target data map(to : cellsX, cellsY, localParameters, localStencil)
    switch (localStencil.getType()) {
    case RHS:
      #pragma omp target map(to : fghLocal, fghLocal.data_[0 : fghLocal.size_]) map(tofrom : rhsLocal, rhsLocal.data_[0 : rhsLocal.size_])
      { 
        for (int j = 1; j < cellsY - 1; j++) {
          for (int i = 1; i < cellsX - 1; i++) {
            localStencil.applyGPU(localParameters, rhsLocal, fghLocal, i, j);
          }
        }
      }
      for (int j = 1 + lowOffset_; j < cellsY - 1 + highOffset_; j++) {
        for (int i = 1 + lowOffset_; i < cellsX - 1 + highOffset_; i++) {
          GPUIterator<FlowFieldType>::flowField_.getRHS().getScalar(i, j) = rhsLocal.getScalar(i, j);
        }
      }
      break;

    case FGH:
      VectorField    velocityLocal = GPUIterator<FlowFieldType>::flowField_.velocity_;
      IntScalarField flagsLocal    = GPUIterator<FlowFieldType>::flowField_.flags_;
      #pragma omp target map(to : velocityLocal, velocityLocal.data_[0 : velocityLocal.size_]) \
        map(to : flagsLocal, flagsLocal.data_[0 : flagsLocal.size_]) \
        map(tofrom : fghLocal, fghLocal.data_[0 : fghLocal.size_])
      {
        for (int j = 1; j < cellsY - 1; j++) {
          for (int i = 1; i < cellsX - 1; i++) {
            localStencil.applyGPU(localParameters, velocityLocal, fghLocal, flagsLocal, i, j);
          }
        }
      }
      for (int j = 1 + lowOffset_; j < cellsY - 1 + highOffset_; j++) {
        for (int i = 1 + lowOffset_; i < cellsX - 1 + highOffset_; i++) {
          GPUIterator<FlowFieldType>::flowField_.getFGH().getVectorElement(i, j, 0) = fghLocal.getVectorElement(i, j, 0);
          GPUIterator<FlowFieldType>::flowField_.getFGH().getVectorElement(i, j, 1) = fghLocal.getVectorElement(i, j, 1);
        }
      }
      break;
    }
  }

  if (GPUIterator<FlowFieldType>::parameters_.geometry.dim == 3) {
    for (int k = 1 + lowOffset_; k < cellsZ - 1 + highOffset_; k++) {
      for (int j = 1 + lowOffset_; j < cellsY - 1 + highOffset_; j++) {
        for (int i = 1 + lowOffset_; i < cellsX - 1 + highOffset_; i++) {
          stencil_.apply(GPUIterator<FlowFieldType>::parameters_, GPUIterator<FlowFieldType>::flowField_, i, j, k);
        }
      }
    }
  }
}

template class GPUFieldIterator<FlowField>;
template class GPUFieldIterator<TurbulentFlowField>;
