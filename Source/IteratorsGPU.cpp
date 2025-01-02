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
  std::cout << stencil_.getType() << std::endl;
  const int cellsX = GPUIterator<FlowFieldType>::flowField_.getCellsX();
  const int cellsY = GPUIterator<FlowFieldType>::flowField_.getCellsY();
  const int cellsZ = GPUIterator<FlowFieldType>::flowField_.getCellsZ();

  auto localFlowField = GPUIterator<FlowFieldType>::flowField_;
  auto localParameters = GPUIterator<FlowFieldType>::parameters_;
  auto localStencil = stencil_;

  // The index k can be used for the 2D and 3D cases.
  if (GPUIterator<FlowFieldType>::parameters_.geometry.dim == 2) {
    // Loop without lower boundaries. These will be dealt with by the global boundary stencils
    // or by the subdomain boundary GPUIterator.

  #pragma omp target parallel for collapse(2) map(to : cellsX, cellsY, lowOffset_, highOffset_, localParameters) map(tofrom : localFlowField, localStencil) 
  {
    for (int j = 1 + lowOffset_; j < cellsY - 1 + highOffset_; j++) {
      for (int i = 1 + lowOffset_; i < cellsX - 1 + highOffset_; i++) {
        localStencil.apply(localParameters, localFlowField, i, j);
      }
    }
  }
    
  }

  if (GPUIterator<FlowFieldType>::parameters_.geometry.dim == 3) {
    for (int k = 1 + lowOffset_; k < cellsZ - 1 + highOffset_; k++) {
      for (int j = 1 + lowOffset_; j < cellsY - 1 + highOffset_; j++) {
        for (int i = 1 + lowOffset_; i < cellsX - 1 + highOffset_; i++) {
          stencil_.apply(localParameters, GPUIterator<FlowFieldType>::flowField_, i, j, k);
        }
      }
    }
  }
}

template class GPUFieldIterator<FlowField>;
template class GPUFieldIterator<TurbulentFlowField>;
