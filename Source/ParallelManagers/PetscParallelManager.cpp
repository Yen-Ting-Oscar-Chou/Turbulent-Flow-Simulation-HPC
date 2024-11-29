#include "StdAfx.hpp"

#include "PetscParallelManager.hpp"
//#include <petsclog.h>

ParallelManagers::PetscParallelManager::PetscParallelManager(Parameters& parameters, FlowField & flowField):
  _parameters(parameters), _flowField(flowField) {
  std::cout << "current rank is:" << _parameters.parallel.rank << std::endl;

   // Handles the dimensions internally.
  _pressureBufferFillStencil = new Stencils::PressureBufferFillStencil(parameters);
  _velocityBufferFillStencil = new Stencils::VelocityBufferFillStencil(parameters);
  _pressureBufferReadStencil = new Stencils::PressureBufferReadStencil(parameters);
  _velocityBufferReadStencil = new Stencils::VelocityBufferReadStencil(parameters);

  int lowOffset = 2;
  int highOffset = -1;

  // Construct iterators.
  _parallelBoundaryPressureFillIterator = new ParallelBoundaryIterator<FlowField>(_flowField, _parameters, *_pressureBufferFillStencil, lowOffset, highOffset);
  _parallelBoundaryPressureReadIterator = new ParallelBoundaryIterator<FlowField>(_flowField, _parameters, *_pressureBufferReadStencil, lowOffset, highOffset);
  
  _parallelBoundaryVelocityFillIterator = new ParallelBoundaryIterator<FlowField>(_flowField, _parameters, *_velocityBufferFillStencil, lowOffset, highOffset);
  _parallelBoundaryVelocityReadIterator = new ParallelBoundaryIterator<FlowField>(_flowField, _parameters, *_velocityBufferReadStencil, lowOffset, highOffset);
}

ParallelManagers::PetscParallelManager::~PetscParallelManager() {
  delete _velocityBufferReadStencil;
  delete _pressureBufferReadStencil;
  delete _velocityBufferFillStencil;
  delete _velocityBufferReadStencil;
  delete _parallelBoundaryPressureFillIterator;
  delete _parallelBoundaryPressureReadIterator;
  delete _parallelBoundaryVelocityFillIterator;
  delete _parallelBoundaryVelocityReadIterator;
}

int ParallelManagers::PetscParallelManager::computeRankFromIndices(int i, int j, int k) const {
  if (i < 0 || i >= _parameters.parallel.numProcessors[0] ||
      j < 0 || j >= _parameters.parallel.numProcessors[1] ||
      k < 0 || k >= _parameters.parallel.numProcessors[2]) {
    return MPI_PROC_NULL;
  }
  int nrank = i + j * _parameters.parallel.numProcessors[0];
  if (_parameters.geometry.dim == 3) {
    nrank += k * _parameters.parallel.numProcessors[0] * _parameters.parallel.numProcessors[1];
  }
  return nrank;
}

void ParallelManagers::PetscParallelManager::communicateVelocities() {
  int counter = 0;
  std::array<MPI_Request, 12> requests{};

  int index0 = _parameters.parallel.indices[0];
  int index1 = _parameters.parallel.indices[1];
  int index2 = 0;
  int left = computeRankFromIndices(index0 - 1, index1, index2);
  int right = computeRankFromIndices(index0 + 1, index1, index2);
  int top = computeRankFromIndices(index0, index1 + 1, index2);
  int bottom = computeRankFromIndices(index0, index1 - 1, index2);
  int front = MPI_PROC_NULL;
  int back = MPI_PROC_NULL;
  if (_parameters.geometry.dim == 3) {
    index2 = _parameters.parallel.indices[2];
    front = computeRankFromIndices(index0, index1, index2 - 1);
    back = computeRankFromIndices(index0, index1, index2 + 1);
  }

  _parallelBoundaryVelocityFillIterator->iterate();
  // Put recv before send to avoid exceptions
  MPI_Irecv(_velocityBufferReadStencil->getvelocityLeft().data(), _velocityBufferReadStencil->getvelocityLeft().size(), MY_MPI_FLOAT, left, 1, PETSC_COMM_WORLD, &requests.at(counter++));
  MPI_Irecv(_velocityBufferReadStencil->getvelocityRight().data(), _velocityBufferReadStencil->getvelocityRight().size(), MY_MPI_FLOAT, right, 1, PETSC_COMM_WORLD, &requests.at(counter++));
  MPI_Irecv(_velocityBufferReadStencil->getvelocityBottom().data(), _velocityBufferReadStencil->getvelocityBottom().size(), MY_MPI_FLOAT, bottom, 1, PETSC_COMM_WORLD, &requests.at(counter++));
  MPI_Irecv(_velocityBufferReadStencil->getvelocityTop().data(), _velocityBufferReadStencil->getvelocityTop().size(), MY_MPI_FLOAT, top, 1, PETSC_COMM_WORLD, &requests.at(counter++));
  MPI_Irecv(_velocityBufferReadStencil->getvelocityFront().data(), _velocityBufferReadStencil->getvelocityFront().size(), MY_MPI_FLOAT, front, 1, PETSC_COMM_WORLD,&requests.at(counter++));
  MPI_Irecv(_velocityBufferReadStencil->getvelocityBack().data(), _velocityBufferReadStencil->getvelocityBack().size(), MY_MPI_FLOAT, back, 1, PETSC_COMM_WORLD,&requests.at(counter++));

  MPI_Isend(_velocityBufferFillStencil->getvelocityLeft().data(), _velocityBufferFillStencil->getvelocityLeft().size(), MY_MPI_FLOAT, left, 1, PETSC_COMM_WORLD, &requests.at(counter++));
  MPI_Isend(_velocityBufferFillStencil->getvelocityRight().data(), _velocityBufferFillStencil->getvelocityRight().size(), MY_MPI_FLOAT, right, 1, PETSC_COMM_WORLD, &requests.at(counter++));
  MPI_Isend(_velocityBufferFillStencil->getvelocityBottom().data(), _velocityBufferFillStencil->getvelocityBottom().size(), MY_MPI_FLOAT, bottom, 1, PETSC_COMM_WORLD, &requests.at(counter++));
  MPI_Isend(_velocityBufferFillStencil->getvelocityTop().data(), _velocityBufferFillStencil->getvelocityTop().size(), MY_MPI_FLOAT, top, 1, PETSC_COMM_WORLD, &requests.at(counter++));
  MPI_Isend(_velocityBufferFillStencil->getvelocityFront().data(), _velocityBufferFillStencil->getvelocityFront().size(), MY_MPI_FLOAT, front, 1, PETSC_COMM_WORLD, &requests.at(counter++));
  MPI_Isend(_velocityBufferFillStencil->getvelocityBack().data(), _velocityBufferFillStencil->getvelocityBack().size(), MY_MPI_FLOAT, back, 1, PETSC_COMM_WORLD, &requests.at(counter++));
  
  MPI_Waitall(counter, requests.data(), MPI_STATUSES_IGNORE);

  /* MPI_Barrier(PETSC_COMM_WORLD);

  if(_parameters.parallel.rank == computeRankFromIndices(0, 0, 0)) {
    std::cout << "bottom left rank is:" << _parameters.parallel.rank << std::endl;
    std::cout << "right send buffer = ";
    for (int i = 0; i < _velocityBufferFillStencil->getvelocityRight().size(); i++) {
      std::cout << _velocityBufferFillStencil->getvelocityRight()[i] << " ";
    }
    std::cout << std::endl;
  }
  
  if(_parameters.parallel.rank == computeRankFromIndices(1, 0, 0)) {
    std::cout << "bottom right rank is:" << _parameters.parallel.rank << std::endl;
    std::cout << "left recv buffer = ";
    for (int i = 0; i < _velocityBufferFillStencil->getvelocityLeft().size(); i++) {
      std::cout << _velocityBufferFillStencil->getvelocityLeft()[i] << " ";
    }
    std::cout << std::endl;
  }

  if(_parameters.parallel.rank == computeRankFromIndices(0, 0, 0)) {
    std::cout << "bottom left rank is:" << _parameters.parallel.rank << std::endl;
    std::cout << "top send buffer = ";
    for (int i = 0; i < _velocityBufferFillStencil->getvelocityTop().size(); i++) {
      std::cout << _velocityBufferFillStencil->getvelocityTop()[i] << " ";
    }
  }
    std::cout << std::endl;
  if(_parameters.parallel.rank == computeRankFromIndices(0, 1, 0)) {
    std::cout << "top left rank is:" << _parameters.parallel.rank << std::endl;
    std::cout << "bottom recv buffer = ";
    for (int i = 0; i < _velocityBufferFillStencil->getvelocityBottom().size(); i++) {
      std::cout << _velocityBufferFillStencil->getvelocityBottom()[i] << " ";
    }
    std::cout << std::endl;
  }

  if(_parameters.parallel.rank == computeRankFromIndices(1, 1, 0)) {
    std::cout << "top right rank is:" << _parameters.parallel.rank << std::endl;
    std::cout << "left send buffer = ";
    for (int i = 0; i < _velocityBufferFillStencil->getvelocityLeft().size(); i++) {
      std::cout << _velocityBufferFillStencil->getvelocityLeft()[i] << " ";
    }
    std::cout << std::endl;
  }
  if(_parameters.parallel.rank == computeRankFromIndices(0, 1, 0)) {
    std::cout << "top left rank is:" << _parameters.parallel.rank << std::endl;
    std::cout << "right recv buffer = ";
    for (int i = 0; i < _velocityBufferFillStencil->getvelocityRight().size(); i++) {
      std::cout << _velocityBufferFillStencil->getvelocityRight()[i] << " ";
    }
    std::cout << std::endl;
  }

  if(_parameters.parallel.rank == computeRankFromIndices(1, 1, 0)) {
    std::cout << "top right rank is:" << _parameters.parallel.rank << std::endl;
    std::cout << "bottom send buffer = ";
    for (int i = 0; i < _velocityBufferFillStencil->getvelocityBottom().size(); i++) {
      std::cout << _velocityBufferFillStencil->getvelocityBottom()[i] << " ";
    }
  }
    std::cout << std::endl;
  if(_parameters.parallel.rank == computeRankFromIndices(1, 0, 0)) {
    std::cout << "bottom right rank is:" << _parameters.parallel.rank << std::endl;
    std::cout << "top recv buffer = ";
    for (int i = 0; i < _velocityBufferFillStencil->getvelocityTop().size(); i++) {
      std::cout << _velocityBufferFillStencil->getvelocityTop()[i] << " ";
    }
    std::cout << std::endl;
  }

  std::cout << "Done printing" << std::endl;
 */
  _parallelBoundaryVelocityReadIterator->iterate();
  
  MPI_Barrier(PETSC_COMM_WORLD);

  // print the pressure values for all ranks, including ghost cells
  if (_parameters.parallel.rank == computeRankFromIndices(0, 0, 0)) {
    std::cout << "current rank is:" << _parameters.parallel.rank << std::endl;
    std::cout << "Velocity u values including ghost cells: " << std::endl;
    
      for (int j = _flowField.getCellsY() - 1; j >= 0; j--) {
          for (int i = 0; i < _flowField.getCellsX(); i++) {
            std::cout << _flowField.getVelocity().getVector(i, j)[0] << " ";
        }
      std::cout << std::endl;
      }
    std::cout << std::endl;

    std::cout << "current rank is:" << _parameters.parallel.rank << std::endl;
    std::cout << "Velocity v values including ghost cells: " << std::endl;
    for (int j = _flowField.getCellsY() - 1; j >= 0; j--) {
      for (int i = 0; i < _flowField.getCellsX(); i++) {
        std::cout << _flowField.getVelocity().getVector(i, j)[1] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }

  MPI_Barrier(PETSC_COMM_WORLD);
  
  if (_parameters.parallel.rank == computeRankFromIndices(1, 0, 0)) {
    std::cout << "current rank is:" << _parameters.parallel.rank << std::endl;
    std::cout << "Velocity u values including ghost cells: " << std::endl;
    
      for (int j = _flowField.getCellsY() - 1; j >= 0; j--) {
          for (int i = 0; i < _flowField.getCellsX(); i++) {
            std::cout << _flowField.getVelocity().getVector(i, j)[0] << " ";
        }
      std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "current rank is:" << _parameters.parallel.rank << std::endl;
    std::cout << "Velocity v values including ghost cells: " << std::endl;
    for (int j = _flowField.getCellsY() - 1; j >= 0; j--) {
      for (int i = 0; i < _flowField.getCellsX(); i++) {
        std::cout << _flowField.getVelocity().getVector(i, j)[1] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
  MPI_Barrier(PETSC_COMM_WORLD);

  if (_parameters.parallel.rank == computeRankFromIndices(0, 1, 0)) {
    std::cout << "current rank is:" << _parameters.parallel.rank << std::endl;
    std::cout << "Velocity values for u including ghost cells: " << std::endl;
    
      for (int j = _flowField.getCellsY() - 1; j >= 0; j--) {
          for (int i = 0; i < _flowField.getCellsX(); i++) {
            std::cout << _flowField.getVelocity().getVector(i, j)[0] << " ";
        }
      std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "current rank is:" << _parameters.parallel.rank << std::endl;
    std::cout << "Velocity v values including ghost cells: " << std::endl;
        for (int j = _flowField.getCellsY() - 1; j >= 0; j--) {
      for (int i = 0; i < _flowField.getCellsX(); i++) {
        std::cout << _flowField.getVelocity().getVector(i, j)[1] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
  MPI_Barrier(PETSC_COMM_WORLD);

  if (_parameters.parallel.rank == computeRankFromIndices(1, 1, 0)) {
    std::cout << "current rank is:" << _parameters.parallel.rank << std::endl;
    std::cout << "Velocity u values including ghost cells: " << std::endl;
    
      for (int j = _flowField.getCellsY() - 1; j >= 0; j--) {
          for (int i = 0; i < _flowField.getCellsX(); i++) {
            std::cout << _flowField.getVelocity().getVector(i, j)[0] << " ";
        }
      std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "current rank is:" << _parameters.parallel.rank << std::endl;
    std::cout << "Velocity v values including ghost cells: " << std::endl;
        for (int j = _flowField.getCellsY() - 1; j >= 0; j--) {
      for (int i = 0; i < _flowField.getCellsX(); i++) {
        std::cout << _flowField.getVelocity().getVector(i, j)[1] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }

}

void ParallelManagers::PetscParallelManager::communicatePressure() {
  int counter = 0;
  std::array<MPI_Request, 12> requests{};
  int index0 = _parameters.parallel.indices[0];
  int index1 = _parameters.parallel.indices[1];
  int index2 = _parameters.parallel.indices[2];
  int left = computeRankFromIndices(index0 - 1, index1, index2);
  int right = computeRankFromIndices(index0 + 1, index1, index2);
  int top = computeRankFromIndices(index0, index1 + 1, index2);
  int bottom = computeRankFromIndices(index0, index1 - 1, index2);
  int front = MPI_PROC_NULL;
  int back = MPI_PROC_NULL;

  if (_parameters.geometry.dim == 3) {
    front = computeRankFromIndices(index0, index1, index2 - 1);
    back = computeRankFromIndices(index0, index1, index2 + 1);
  }

  _parallelBoundaryPressureFillIterator->iterate();
  
/* 
   if(_parameters.parallel.rank == 0) {
      std::cout << "current rank is:" << _parameters.parallel.rank << std::endl;
      
      std::cout << "left fill stencil = ";
      for (int i = 0; i < _pressureBufferFillStencil..size(); i++) {
        std::cout << _velocityBufferFillStencil->getvelocityLeft()[i] << " ";
      }
      std::cout << std::endl;
      std::cout << "right fille buffer = ";
      for (int i = 0; i < _pressureBufferStencil.rightSend.size(); i++) {
        std::cout << pressureBuffers.rightSend[i] << " ";
      }
      std::cout << std::endl;
      std::cout << "right fill stencil = ";
      for (int i = 0; i < pressureBuffers.leftSend.size(); i++) {
        std::cout << _velocityBufferReadStencil->getvelocityRight()[i] << " ";
      }
      std::cout << std::endl;
      std::cout << "top send buffer = ";
      for (int i = 0; i < pressureBuffers.topSend.size(); i++) {
        std::cout << pressureBuffers.topSend[i] << " ";
      }
      std::cout << std::endl;
      std::cout << "top fill stencil = ";
      for (int i = 0; i < pressureBuffers.topSend.size(); i++) {
        std::cout << _velocityBufferReadStencil->getvelocityTop()[i] << " ";
      }
      std::cout << std::endl;
      std::cout << "bottom send buffer = ";
      for (int i = 0; i < pressureBuffers.bottomSend.size(); i++) {
        std::cout << pressureBuffers.bottomSend[i] << " ";
      }
      std::cout << std::endl;
      std::cout << "bottom fill stencil = ";
      for (int i = 0; i < pressureBuffers.bottomSend.size(); i++) {
        std::cout << _velocityBufferReadStencil->getvelocityBottom()[i] << " ";
      }
      std::cout << std::endl;
      std::cout << "Done printing" << std::endl;
  }  */
  
  // Put recv before send to avoid exceptions
  MPI_Irecv(_pressureBufferReadStencil->getpressureLeft().data(), _pressureBufferReadStencil->getpressureLeft().size(), MY_MPI_FLOAT, left, 1, PETSC_COMM_WORLD, &requests.at(counter++));
  MPI_Irecv(_pressureBufferReadStencil->getpressureRight().data(), _pressureBufferReadStencil->getpressureRight().size(), MY_MPI_FLOAT, right, 1, PETSC_COMM_WORLD, &requests.at(counter++));
  MPI_Irecv(_pressureBufferReadStencil->getpressureBottom().data(), _pressureBufferReadStencil->getpressureBottom().size(), MY_MPI_FLOAT, bottom, 1, PETSC_COMM_WORLD, &requests.at(counter++));
  MPI_Irecv(_pressureBufferReadStencil->getpressureTop().data(), _pressureBufferReadStencil->getpressureTop().size(), MY_MPI_FLOAT, top, 1, PETSC_COMM_WORLD, &requests.at(counter++));
  MPI_Irecv(_pressureBufferReadStencil->getpressureFront().data(), _pressureBufferReadStencil->getpressureFront().size(), MY_MPI_FLOAT, front, 1, PETSC_COMM_WORLD, &requests.at(counter++));
  MPI_Irecv(_pressureBufferReadStencil->getpressureBack().data(), _pressureBufferReadStencil->getpressureBack().size(), MY_MPI_FLOAT, back, 1, PETSC_COMM_WORLD, &requests.at(counter++));
  
  MPI_Isend(_pressureBufferFillStencil->getpressureLeft().data(), _pressureBufferFillStencil->getpressureLeft().size(), MY_MPI_FLOAT, left, 1, PETSC_COMM_WORLD, &requests.at(counter++));
  MPI_Isend(_pressureBufferFillStencil->getpressureRight().data(), _pressureBufferFillStencil->getpressureRight().size(), MY_MPI_FLOAT, right, 1, PETSC_COMM_WORLD, &requests.at(counter++));
  MPI_Isend(_pressureBufferFillStencil->getpressureBottom().data(), _pressureBufferFillStencil->getpressureBottom().size(), MY_MPI_FLOAT, bottom, 1, PETSC_COMM_WORLD, &requests.at(counter++));
  MPI_Isend(_pressureBufferFillStencil->getpressureTop().data(), _pressureBufferFillStencil->getpressureTop().size(), MY_MPI_FLOAT, top, 1, PETSC_COMM_WORLD, &requests.at(counter++));
  MPI_Isend(_pressureBufferFillStencil->getpressureFront().data(), _pressureBufferFillStencil->getpressureFront().size(), MY_MPI_FLOAT, front, 1, PETSC_COMM_WORLD, &requests.at(counter++));
  MPI_Isend(_pressureBufferFillStencil->getpressureBack().data(), _pressureBufferFillStencil->getpressureBack().size(), MY_MPI_FLOAT, back, 1, PETSC_COMM_WORLD, &requests.at(counter++));

  MPI_Waitall(counter, requests.data(), MPI_STATUSES_IGNORE);
  /*
  // print the pressure values for all ranks, including ghost cells
  if (_parameters.parallel.rank == computeRankFromIndices(0, 0, 0)) {
    std::cout << "current rank is:" << _parameters.parallel.rank << std::endl;
    std::cout << "Pressure values including ghost cells: " << std::endl;
    
      for (int j = _flowField.getCellsY() - 1; j >= 0; j--) {
          for (int i = 0; i < _flowField.getCellsX(); i++) {
            std::cout << _flowField.getPressure().getScalar(i, j) << " ";
        }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
 
  MPI_Barrier(PETSC_COMM_WORLD);
  
  if (_parameters.parallel.rank == computeRankFromIndices(1, 0, 0)) {
    std::cout << "current rank is:" << _parameters.parallel.rank << std::endl;
    std::cout << "Pressure values including ghost cells: " << std::endl;
    
      for (int j = _flowField.getCellsY() - 1; j >= 0; j--) {
          for (int i = 0; i < _flowField.getCellsX(); i++) {
            std::cout << _flowField.getPressure().getScalar(i, j) << " ";
        }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
  MPI_Barrier(PETSC_COMM_WORLD);

  if (_parameters.parallel.rank == computeRankFromIndices(0, 1, 0)) {
    std::cout << "current rank is:" << _parameters.parallel.rank << std::endl;
    std::cout << "Pressure values including ghost cells: " << std::endl;
    
      for (int j = _flowField.getCellsY() - 1; j >= 0; j--) {
          for (int i = 0; i < _flowField.getCellsX(); i++) {
            std::cout << _flowField.getPressure().getScalar(i, j) << " ";
        }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
  MPI_Barrier(PETSC_COMM_WORLD);

  if (_parameters.parallel.rank == computeRankFromIndices(1, 1, 0)) {
    std::cout << "current rank is:" << _parameters.parallel.rank << std::endl;
    std::cout << "Pressure values including ghost cells: " << std::endl;
    
      for (int j = _flowField.getCellsY() - 1; j >= 0; j--) {
          for (int i = 0; i < _flowField.getCellsX(); i++) {
            std::cout << _flowField.getPressure().getScalar(i, j) << " ";
        }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
  MPI_Barrier(PETSC_COMM_WORLD);
  */
  // 
  /*
  if(_parameters.parallel.rank == computeRankFromIndices(0, 0, 0)) {
    std::cout << "bottom left rank is:" << _parameters.parallel.rank << std::endl;
    std::cout << "right send buffer = ";
    for (int i = 0; i < _pressureBufferFillStencil->getpressureRight().size(); i++) {
      std::cout << _pressureBufferFillStencil->getpressureRight()[i] << " ";
    }
    std::cout << std::endl;
  }
  
  if(_parameters.parallel.rank == computeRankFromIndices(1, 0, 0)) {
    std::cout << "bottom right rank is:" << _parameters.parallel.rank << std::endl;
    std::cout << "left recv buffer = ";
    for (int i = 0; i < _pressureBufferReadStencil->getpressureLeft().size(); i++) {
      std::cout << _pressureBufferReadStencil->getpressureLeft()[i] << " ";
    }
    std::cout << std::endl;
  }

  if(_parameters.parallel.rank == computeRankFromIndices(0, 0, 0)) {
    std::cout << "bottom left rank is:" << _parameters.parallel.rank << std::endl;
    std::cout << "top send buffer = ";
    for (int i = 0; i < _pressureBufferFillStencil->getpressureTop().size(); i++) {
      std::cout << _pressureBufferFillStencil->getpressureTop()[i] << " ";
    }
  }
    std::cout << std::endl;
  if(_parameters.parallel.rank == computeRankFromIndices(0, 1, 0)) {
    std::cout << "top left rank is:" << _parameters.parallel.rank << std::endl;
    std::cout << "bottom recv buffer = ";
    for (int i = 0; i < _pressureBufferReadStencil->getpressureBottom().size(); i++) {
      std::cout << _pressureBufferReadStencil->getpressureBottom()[i] << " ";
    }
    std::cout << std::endl;
  }

  if(_parameters.parallel.rank == computeRankFromIndices(1, 1, 0)) {
    std::cout << "top right rank is:" << _parameters.parallel.rank << std::endl;
    std::cout << "left send buffer = ";
    for (int i = 0; i < _pressureBufferFillStencil->getpressureLeft().size(); i++) {
      std::cout << _pressureBufferFillStencil->getpressureLeft()[i] << " ";
    }
    std::cout << std::endl;
  }
  if(_parameters.parallel.rank == computeRankFromIndices(0, 1, 0)) {
    std::cout << "top left rank is:" << _parameters.parallel.rank << std::endl;
    std::cout << "right recv buffer = ";
    for (int i = 0; i < _pressureBufferReadStencil->getpressureRight().size(); i++) {
      std::cout << _pressureBufferReadStencil->getpressureRight()[i] << " ";
    }
    std::cout << std::endl;
  }

  if(_parameters.parallel.rank == computeRankFromIndices(1, 1, 0)) {
    std::cout << "top right rank is:" << _parameters.parallel.rank << std::endl;
    std::cout << "bottom send buffer = ";
    for (int i = 0; i < _pressureBufferFillStencil->getpressureBottom().size(); i++) {
      std::cout << _pressureBufferFillStencil->getpressureBottom()[i] << " ";
    }
  }
    std::cout << std::endl;
  if(_parameters.parallel.rank == computeRankFromIndices(1, 0, 0)) {
    std::cout << "bottom right rank is:" << _parameters.parallel.rank << std::endl;
    std::cout << "top recv buffer = ";
    for (int i = 0; i < _pressureBufferReadStencil->getpressureTop().size(); i++) {
      std::cout << _pressureBufferReadStencil->getpressureTop()[i] << " ";
    }
    std::cout << std::endl;
  }
  std::cout << "Done printing" << std::endl;
  */
  _parallelBoundaryPressureReadIterator->iterate();
  /*
  MPI_Barrier(PETSC_COMM_WORLD);
  // print the pressure values for all ranks, including ghost cells
  if (_parameters.parallel.rank == computeRankFromIndices(0, 0, 0)) {
    std::cout << "current rank is:" << _parameters.parallel.rank << std::endl;
    std::cout << "Pressure values including ghost cells: " << std::endl;
    
      for (int j = _flowField.getCellsY() - 1; j >= 0; j--) {
          for (int i = 0; i < _flowField.getCellsX(); i++) {
            std::cout << _flowField.getPressure().getScalar(i, j) << " ";
        }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
 
  MPI_Barrier(PETSC_COMM_WORLD);
  
  if (_parameters.parallel.rank == computeRankFromIndices(1, 0, 0)) {
    std::cout << "current rank is:" << _parameters.parallel.rank << std::endl;
    std::cout << "Pressure values including ghost cells: " << std::endl;
    
      for (int j = _flowField.getCellsY() - 1; j >= 0; j--) {
          for (int i = 0; i < _flowField.getCellsX(); i++) {
            std::cout << _flowField.getPressure().getScalar(i, j) << " ";
        }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
  MPI_Barrier(PETSC_COMM_WORLD);

  if (_parameters.parallel.rank == computeRankFromIndices(0, 1, 0)) {
    std::cout << "current rank is:" << _parameters.parallel.rank << std::endl;
    std::cout << "Pressure values including ghost cells: " << std::endl;
    
      for (int j = _flowField.getCellsY() - 1; j >= 0; j--) {
          for (int i = 0; i < _flowField.getCellsX(); i++) {
            std::cout << _flowField.getPressure().getScalar(i, j) << " ";
        }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
  MPI_Barrier(PETSC_COMM_WORLD);

  if (_parameters.parallel.rank == computeRankFromIndices(1, 1, 0)) {
    std::cout << "current rank is:" << _parameters.parallel.rank << std::endl;
    std::cout << "Pressure values including ghost cells: " << std::endl;
    
      for (int j = _flowField.getCellsY() - 1; j >= 0; j--) {
          for (int i = 0; i < _flowField.getCellsX(); i++) {
            std::cout << _flowField.getPressure().getScalar(i, j) << " ";
        }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
  MPI_Barrier(PETSC_COMM_WORLD);
  */
}
