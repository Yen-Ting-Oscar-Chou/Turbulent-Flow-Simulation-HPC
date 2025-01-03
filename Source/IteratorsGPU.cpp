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

  FlowField            localFlowField  = GPUIterator<FlowFieldType>::flowField_;
  Parameters           localParameters = GPUIterator<FlowFieldType>::parameters_;
  FieldStencilDelegate localStencil    = stencil_;

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

    switch (stencil_.getType()) {
    case RHS:
      for (int j = 1 + lowOffset_; j < cellsY - 1 + highOffset_; j++) {
        for (int i = 1 + lowOffset_; i < cellsX - 1 + highOffset_; i++) {
          GPUIterator<FlowFieldType>::flowField_.getRHS().getScalar(i, j) = localFlowField.getRHS().getScalar(i, j);
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

  //delete localFlowField;
}

template class GPUFieldIterator<FlowField>;
template class GPUFieldIterator<TurbulentFlowField>;
