#include "StdAfx.hpp"

#include "Iterators.hpp"

template <class FlowFieldType>
FieldIterator<FlowFieldType>::FieldIterator(FlowFieldType& flowField, const Parameters& parameters, Stencils::FieldStencil<FlowFieldType>& stencil, int lowOffset, int highOffset):
  Iterator<FlowFieldType>(flowField, parameters),
  stencil_(stencil),
  lowOffset_(lowOffset),
  highOffset_(highOffset) {}

template <class FlowFieldType>
void FieldIterator<FlowFieldType>::iterate() {
  const int cellsX = Iterator<FlowFieldType>::flowField_.getCellsX();
  const int cellsY = Iterator<FlowFieldType>::flowField_.getCellsY();
  const int cellsZ = Iterator<FlowFieldType>::flowField_.getCellsZ();
  // The index k can be used for the 2D and 3D cases.
  if (Iterator<FlowFieldType>::parameters_.geometry.dim == 2) {
    // Loop without lower boundaries. These will be dealt with by the global boundary stencils
    // or by the subdomain boundary iterators.
    for (int j = 1 + lowOffset_; j < cellsY - 1 + highOffset_; j++) {
      for (int i = 1 + lowOffset_; i < cellsX - 1 + highOffset_; i++) {
        stencil_.apply(Iterator<FlowFieldType>::parameters_, Iterator<FlowFieldType>::flowField_, i, j);
      }
    }
  }

  if (Iterator<FlowFieldType>::parameters_.geometry.dim == 3) {
    for (int k = 1 + lowOffset_; k < cellsZ - 1 + highOffset_; k++) {
      for (int j = 1 + lowOffset_; j < cellsY - 1 + highOffset_; j++) {
        for (int i = 1 + lowOffset_; i < cellsX - 1 + highOffset_; i++) {
          stencil_.apply(Iterator<FlowFieldType>::parameters_, Iterator<FlowFieldType>::flowField_, i, j, k);
        }
      }
    }
  }
}

template <class FlowFieldType>
FieldIteratorGPU<FlowFieldType>::FieldIteratorGPU(int lowOffset, int highOffset):
  lowOffset_(lowOffset),
  highOffset_(highOffset) {}

template <class FlowFieldType>
void FieldIteratorGPU<FlowFieldType>::iterate(StencilType type, const Parameters& parameters, FlowFieldType& flowField, StencilDelegate& stencil) {
  const int cellsX = flowField.getCellsX();
  const int cellsY = flowField.getCellsY();
  const int cellsZ = flowField.getCellsZ();
  stencil.setType(type);
  // The index k can be used for the 2D and 3D cases.
  if (parameters.geometry.dim == 2) {
    // Loop without lower boundaries. These will be dealt with by the global boundary stencils
    // or by the subdomain boundary iterators.
    for (int j = 1 + lowOffset_; j < cellsY - 1 + highOffset_; j++) {
      for (int i = 1 + lowOffset_; i < cellsX - 1 + highOffset_; i++) {
        stencil.apply(parameters, flowField, i, j);
      }
    }
  }

  if (parameters.geometry.dim == 3) {
    for (int k = 1 + lowOffset_; k < cellsZ - 1 + highOffset_; k++) {
      for (int j = 1 + lowOffset_; j < cellsY - 1 + highOffset_; j++) {
        for (int i = 1 + lowOffset_; i < cellsX - 1 + highOffset_; i++) {
          stencil.apply(parameters, flowField, i, j, k);
        }
      }
    }
  }
}

template <class FlowFieldType>
GlobalBoundaryIteratorGPU<FlowFieldType>::GlobalBoundaryIteratorGPU(
  FlowFieldType& flowField, const Parameters& parameters, StencilDelegate& stencil, int lowOffset, int highOffset
):
  flowField_(flowField),
  parameters_(parameters),
  stencil_(stencil),
  lowOffset_(lowOffset),
  highOffset_(highOffset) {}


template <class FlowFieldType>
void GlobalBoundaryIteratorGPU<FlowFieldType>::iterate(StencilType type) {
  stencil_.setType(type);
  if (parameters_.geometry.dim == 2) {
    if (parameters_.parallel.leftNb < 0) {
      for (int j = lowOffset_; j < flowField_.getCellsY() + highOffset_; j++) {
        stencil_.applyLeftWall(parameters_, flowField_, lowOffset_, j);
      }
    }

    if (parameters_.parallel.rightNb < 0) {
      for (int j = lowOffset_; j < flowField_.getCellsY() + highOffset_; j++) {
        stencil_.applyRightWall(parameters_, flowField_, flowField_.getCellsX() + highOffset_ - 1, j);
      }
    }

    if (parameters_.parallel.bottomNb < 0) {
      for (int i = lowOffset_; i < flowField_.getCellsX() + highOffset_; i++) {
        stencil_.applyBottomWall(parameters_, flowField_, i, lowOffset_);
      }
    }

    if (parameters_.parallel.topNb < 0) {
      for (int i = lowOffset_; i < flowField_.getCellsX() + highOffset_; i++) {
        stencil_.applyTopWall(parameters_, flowField_, i, flowField_.getCellsY() + highOffset_ - 1);
      }
    }
  }

  if (parameters_.geometry.dim == 3) {
    if (parameters_.parallel.leftNb < 0) {
      for (int j = lowOffset_; j < flowField_.getCellsY() + highOffset_; j++) {
        for (int k = lowOffset_; k < flowField_.getCellsZ() + highOffset_; k++) {
          stencil_.applyLeftWall(parameters_, flowField_, lowOffset_, j, k);
        }
      }
    }

    if (parameters_.parallel.rightNb < 0) {
      for (int j = lowOffset_; j < flowField_.getCellsY() + highOffset_; j++) {
        for (int k = lowOffset_; k < flowField_.getCellsZ() + highOffset_; k++) {
          stencil_.applyRightWall(parameters_, flowField_, flowField_.getCellsX() + highOffset_ - 1, j, k);
        }
      }
    }

    if (parameters_.parallel.bottomNb < 0) {
      for (int i = lowOffset_; i < flowField_.getCellsX() + highOffset_; i++) {
        for (int k = lowOffset_; k < flowField_.getCellsZ() + highOffset_; k++) {
          stencil_.applyBottomWall(parameters_, flowField_, i, lowOffset_, k);
        }
      }
    }

    if (parameters_.parallel.topNb < 0) {
      for (int i = lowOffset_; i < flowField_.getCellsX() + highOffset_; i++) {
        for (int k = lowOffset_; k < flowField_.getCellsZ() + highOffset_; k++) {
          stencil_.applyTopWall(parameters_, flowField_, i, flowField_.getCellsY() + highOffset_ - 1, k);
        }
      }
    }

    if (parameters_.parallel.frontNb < 0) {
      for (int i = lowOffset_; i < flowField_.getCellsX() + highOffset_; i++) {
        for (int j = lowOffset_; j < flowField_.getCellsY() + highOffset_; j++) {
          stencil_.applyFrontWall(parameters_, flowField_, i, j, lowOffset_);
        }
      }
    }

    if (parameters_.parallel.backNb < 0) {
      for (int i = lowOffset_; i < flowField_.getCellsX() + highOffset_; i++) {
        for (int j = lowOffset_; j < flowField_.getCellsY() + highOffset_; j++) {
          stencil_.applyBackWall(parameters_, flowField_, i, j, flowField_.getCellsZ() + highOffset_ - 1);
        }
      }
    }
  }
}

template <class FlowFieldType>
GlobalBoundaryIterator<FlowFieldType>::GlobalBoundaryIterator(
  FlowFieldType& flowField, const Parameters& parameters, Stencils::BoundaryStencil<FlowFieldType>& stencil, int lowOffset, int highOffset
):
  Iterator<FlowFieldType>(flowField, parameters),
  lowOffset_(lowOffset),
  highOffset_(highOffset),
  leftWallStencil_(stencil),
  rightWallStencil_(stencil),
  bottomWallStencil_(stencil),
  topWallStencil_(stencil),
  frontWallStencil_(stencil),
  backWallStencil_(stencil) {}

template <class FlowFieldType>
GlobalBoundaryIterator<FlowFieldType>::GlobalBoundaryIterator(
  FlowFieldType&                            flowField,
  const Parameters&                         parameters,
  Stencils::BoundaryStencil<FlowFieldType>& leftWallStencil,
  Stencils::BoundaryStencil<FlowFieldType>& rightWallStencil,
  Stencils::BoundaryStencil<FlowFieldType>& bottomWallStencil,
  Stencils::BoundaryStencil<FlowFieldType>& topWallStencil,
  int                                       lowOffset,
  int                                       highOffset
):
  Iterator<FlowFieldType>(flowField, parameters),
  lowOffset_(lowOffset),
  highOffset_(highOffset),
  leftWallStencil_(leftWallStencil),
  rightWallStencil_(rightWallStencil),
  bottomWallStencil_(bottomWallStencil),
  topWallStencil_(topWallStencil)
  // This is plain bad, but it will work. The references had to be initialized somehow.
  ,
  frontWallStencil_(leftWallStencil),
  backWallStencil_(leftWallStencil) {

  if (Iterator<FlowFieldType>::parameters_.geometry.dim == 3) {
    throw std::runtime_error("Trying to use 2D constructor for a 3D field");
  }
}

template <class FlowFieldType>
GlobalBoundaryIterator<FlowFieldType>::GlobalBoundaryIterator(
  FlowFieldType&                            flowField,
  const Parameters&                         parameters,
  Stencils::BoundaryStencil<FlowFieldType>& leftWallStencil,
  Stencils::BoundaryStencil<FlowFieldType>& rightWallStencil,
  Stencils::BoundaryStencil<FlowFieldType>& bottomWallStencil,
  Stencils::BoundaryStencil<FlowFieldType>& topWallStencil,
  Stencils::BoundaryStencil<FlowFieldType>& frontWallStencil,
  Stencils::BoundaryStencil<FlowFieldType>& backWallStencil,
  int                                       lowOffset,
  int                                       highOffset
):
  Iterator<FlowFieldType>(flowField, parameters),
  lowOffset_(lowOffset),
  highOffset_(highOffset),
  leftWallStencil_(leftWallStencil),
  rightWallStencil_(rightWallStencil),
  bottomWallStencil_(bottomWallStencil),
  topWallStencil_(topWallStencil),
  frontWallStencil_(frontWallStencil),
  backWallStencil_(backWallStencil) {}

template <class FlowFieldType>
void GlobalBoundaryIterator<FlowFieldType>::iterate() {
  if (Iterator<FlowFieldType>::parameters_.geometry.dim == 2) {
    if (Iterator<FlowFieldType>::parameters_.parallel.leftNb < 0) {
      for (int j = lowOffset_; j < Iterator<FlowFieldType>::flowField_.getCellsY() + highOffset_; j++) {
        leftWallStencil_.applyLeftWall(Iterator<FlowFieldType>::parameters_, Iterator<FlowFieldType>::flowField_, lowOffset_, j);
      }
    }

    if (Iterator<FlowFieldType>::parameters_.parallel.rightNb < 0) {
      for (int j = lowOffset_; j < Iterator<FlowFieldType>::flowField_.getCellsY() + highOffset_; j++) {
        rightWallStencil_.applyRightWall(
          Iterator<FlowFieldType>::parameters_, Iterator<FlowFieldType>::flowField_, Iterator<FlowFieldType>::flowField_.getCellsX() + highOffset_ - 1, j
        );
      }
    }

    if (Iterator<FlowFieldType>::parameters_.parallel.bottomNb < 0) {
      for (int i = lowOffset_; i < Iterator<FlowFieldType>::flowField_.getCellsX() + highOffset_; i++) {
        bottomWallStencil_.applyBottomWall(Iterator<FlowFieldType>::parameters_, Iterator<FlowFieldType>::flowField_, i, lowOffset_);
      }
    }

    if (Iterator<FlowFieldType>::parameters_.parallel.topNb < 0) {
      for (int i = lowOffset_; i < Iterator<FlowFieldType>::flowField_.getCellsX() + highOffset_; i++) {
        topWallStencil_.applyTopWall(
          Iterator<FlowFieldType>::parameters_, Iterator<FlowFieldType>::flowField_, i, Iterator<FlowFieldType>::flowField_.getCellsY() + highOffset_ - 1
        );
      }
    }
  }

  if (Iterator<FlowFieldType>::parameters_.geometry.dim == 3) {
    if (Iterator<FlowFieldType>::parameters_.parallel.leftNb < 0) {
      for (int j = lowOffset_; j < Iterator<FlowFieldType>::flowField_.getCellsY() + highOffset_; j++) {
        for (int k = lowOffset_; k < Iterator<FlowFieldType>::flowField_.getCellsZ() + highOffset_; k++) {
          leftWallStencil_.applyLeftWall(Iterator<FlowFieldType>::parameters_, Iterator<FlowFieldType>::flowField_, lowOffset_, j, k);
        }
      }
    }

    if (Iterator<FlowFieldType>::parameters_.parallel.rightNb < 0) {
      for (int j = lowOffset_; j < Iterator<FlowFieldType>::flowField_.getCellsY() + highOffset_; j++) {
        for (int k = lowOffset_; k < Iterator<FlowFieldType>::flowField_.getCellsZ() + highOffset_; k++) {
          rightWallStencil_.applyRightWall(
            Iterator<FlowFieldType>::parameters_, Iterator<FlowFieldType>::flowField_, Iterator<FlowFieldType>::flowField_.getCellsX() + highOffset_ - 1, j, k
          );
        }
      }
    }

    if (Iterator<FlowFieldType>::parameters_.parallel.bottomNb < 0) {
      for (int i = lowOffset_; i < Iterator<FlowFieldType>::flowField_.getCellsX() + highOffset_; i++) {
        for (int k = lowOffset_; k < Iterator<FlowFieldType>::flowField_.getCellsZ() + highOffset_; k++) {
          bottomWallStencil_.applyBottomWall(Iterator<FlowFieldType>::parameters_, Iterator<FlowFieldType>::flowField_, i, lowOffset_, k);
        }
      }
    }

    if (Iterator<FlowFieldType>::parameters_.parallel.topNb < 0) {
      for (int i = lowOffset_; i < Iterator<FlowFieldType>::flowField_.getCellsX() + highOffset_; i++) {
        for (int k = lowOffset_; k < Iterator<FlowFieldType>::flowField_.getCellsZ() + highOffset_; k++) {
          topWallStencil_.applyTopWall(
            Iterator<FlowFieldType>::parameters_, Iterator<FlowFieldType>::flowField_, i, Iterator<FlowFieldType>::flowField_.getCellsY() + highOffset_ - 1, k
          );
        }
      }
    }

    if (Iterator<FlowFieldType>::parameters_.parallel.frontNb < 0) {
      for (int i = lowOffset_; i < Iterator<FlowFieldType>::flowField_.getCellsX() + highOffset_; i++) {
        for (int j = lowOffset_; j < Iterator<FlowFieldType>::flowField_.getCellsY() + highOffset_; j++) {
          frontWallStencil_.applyFrontWall(Iterator<FlowFieldType>::parameters_, Iterator<FlowFieldType>::flowField_, i, j, lowOffset_);
        }
      }
    }

    if (Iterator<FlowFieldType>::parameters_.parallel.backNb < 0) {
      for (int i = lowOffset_; i < Iterator<FlowFieldType>::flowField_.getCellsX() + highOffset_; i++) {
        for (int j = lowOffset_; j < Iterator<FlowFieldType>::flowField_.getCellsY() + highOffset_; j++) {
          backWallStencil_.applyBackWall(
            Iterator<FlowFieldType>::parameters_, Iterator<FlowFieldType>::flowField_, i, j, Iterator<FlowFieldType>::flowField_.getCellsZ() + highOffset_ - 1
          );
        }
      }
    }
  }
}

template <class FlowFieldType>
ParallelBoundaryIterator<FlowFieldType>::ParallelBoundaryIterator(
  FlowFieldType& flowField, const Parameters& parameters, Stencils::BoundaryStencil<FlowFieldType>& stencil, int lowOffset, int highOffset
):
  Iterator<FlowFieldType>(flowField, parameters),
  stencil_(stencil),
  lowOffset_(lowOffset),
  highOffset_(highOffset) {}

template <class FlowFieldType>
void ParallelBoundaryIterator<FlowFieldType>::iterate() {
  if (Iterator<FlowFieldType>::parameters_.geometry.dim == 2) {
    if (Iterator<FlowFieldType>::parameters_.parallel.leftNb >= 0) {
      for (int j = lowOffset_; j < Iterator<FlowFieldType>::flowField_.getCellsY() + highOffset_; j++) {
        stencil_.applyLeftWall(Iterator<FlowFieldType>::parameters_, Iterator<FlowFieldType>::flowField_, lowOffset_, j);
      }
    }

    if (Iterator<FlowFieldType>::parameters_.parallel.rightNb >= 0) {
      for (int j = lowOffset_; j < Iterator<FlowFieldType>::flowField_.getCellsY() + highOffset_; j++) {
        stencil_.applyRightWall(Iterator<FlowFieldType>::parameters_, Iterator<FlowFieldType>::flowField_, Iterator<FlowFieldType>::flowField_.getCellsX() + highOffset_ - 1, j);
      }
    }

    if (Iterator<FlowFieldType>::parameters_.parallel.bottomNb >= 0) {
      for (int i = lowOffset_; i < Iterator<FlowFieldType>::flowField_.getCellsX() + highOffset_; i++) {
        stencil_.applyBottomWall(Iterator<FlowFieldType>::parameters_, Iterator<FlowFieldType>::flowField_, i, lowOffset_);
      }
    }

    if (Iterator<FlowFieldType>::parameters_.parallel.topNb >= 0) {
      for (int i = lowOffset_; i < Iterator<FlowFieldType>::flowField_.getCellsX() + highOffset_; i++) {
        stencil_.applyTopWall(Iterator<FlowFieldType>::parameters_, Iterator<FlowFieldType>::flowField_, i, Iterator<FlowFieldType>::flowField_.getCellsY() + highOffset_ - 1);
      }
    }
  }

  if (Iterator<FlowFieldType>::parameters_.geometry.dim == 3) {
    if (Iterator<FlowFieldType>::parameters_.parallel.leftNb >= 0) {
      for (int j = lowOffset_; j < Iterator<FlowFieldType>::flowField_.getCellsY() + highOffset_; j++) {
        for (int k = lowOffset_; k < Iterator<FlowFieldType>::flowField_.getCellsZ() + highOffset_; k++) {
          stencil_.applyLeftWall(Iterator<FlowFieldType>::parameters_, Iterator<FlowFieldType>::flowField_, lowOffset_, j, k);
        }
      }
    }

    if (Iterator<FlowFieldType>::parameters_.parallel.rightNb >= 0) {
      for (int j = lowOffset_; j < Iterator<FlowFieldType>::flowField_.getCellsY() + highOffset_; j++) {
        for (int k = lowOffset_; k < Iterator<FlowFieldType>::flowField_.getCellsZ() + highOffset_; k++) {
          stencil_.applyRightWall(
            Iterator<FlowFieldType>::parameters_, Iterator<FlowFieldType>::flowField_, Iterator<FlowFieldType>::flowField_.getCellsX() + highOffset_ - 1, j, k
          );
        }
      }
    }

    if (Iterator<FlowFieldType>::parameters_.parallel.bottomNb >= 0) {
      for (int i = lowOffset_; i < Iterator<FlowFieldType>::flowField_.getCellsX() + highOffset_; i++) {
        for (int k = lowOffset_; k < Iterator<FlowFieldType>::flowField_.getCellsZ() + highOffset_; k++) {
          stencil_.applyBottomWall(Iterator<FlowFieldType>::parameters_, Iterator<FlowFieldType>::flowField_, i, lowOffset_, k);
        }
      }
    }

    if (Iterator<FlowFieldType>::parameters_.parallel.topNb >= 0) {
      for (int i = lowOffset_; i < Iterator<FlowFieldType>::flowField_.getCellsX() + highOffset_; i++) {
        for (int k = lowOffset_; k < Iterator<FlowFieldType>::flowField_.getCellsZ() + highOffset_; k++) {
          stencil_.applyTopWall(Iterator<FlowFieldType>::parameters_, Iterator<FlowFieldType>::flowField_, i, Iterator<FlowFieldType>::flowField_.getCellsY() + highOffset_ - 1, k);
        }
      }
    }

    if (Iterator<FlowFieldType>::parameters_.parallel.frontNb >= 0) {
      for (int i = lowOffset_; i < Iterator<FlowFieldType>::flowField_.getCellsX() + highOffset_; i++) {
        for (int j = lowOffset_; j < Iterator<FlowFieldType>::flowField_.getCellsY() + highOffset_; j++) {
          stencil_.applyFrontWall(Iterator<FlowFieldType>::parameters_, Iterator<FlowFieldType>::flowField_, i, j, lowOffset_);
        }
      }
    }

    if (Iterator<FlowFieldType>::parameters_.parallel.backNb >= 0) {
      for (int i = lowOffset_; i < Iterator<FlowFieldType>::flowField_.getCellsX() + highOffset_; i++) {
        for (int j = lowOffset_; j < Iterator<FlowFieldType>::flowField_.getCellsY() + highOffset_; j++) {
          stencil_.applyBackWall(
            Iterator<FlowFieldType>::parameters_, Iterator<FlowFieldType>::flowField_, i, j, Iterator<FlowFieldType>::flowField_.getCellsZ() + highOffset_ - 1
          );
        }
      }
    }
  }
}

template class FieldIterator<FlowField>;
template class FieldIterator<TurbulentFlowField>;
template class GlobalBoundaryIterator<FlowField>;
template class GlobalBoundaryIterator<TurbulentFlowField>;
template class FieldIteratorGPU<FlowField>;
template class FieldIteratorGPU<TurbulentFlowField>;
template class GlobalBoundaryIteratorGPU<FlowField>;
template class GlobalBoundaryIteratorGPU<TurbulentFlowField>;
template class ParallelBoundaryIterator<FlowField>;
template class ParallelBoundaryIterator<TurbulentFlowField>;
