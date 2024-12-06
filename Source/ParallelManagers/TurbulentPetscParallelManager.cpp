#include "StdAfx.hpp"

#include "TurbulentPetscParallelManager.hpp"

ParallelManagers::TurbulentPetscParallelManager::TurbulentPetscParallelManager(Parameters& parameters, TurbulentFlowField& flowField):
  PetscParallelManager(parameters, flowField),
  _turbulentFlowField(flowField),
  _viscosityBufferFillStencil(parameters),
  _viscosityBufferReadStencil(parameters),
  _parallelBoundaryViscosityFillIterator(_turbulentFlowField, parameters, _viscosityBufferFillStencil, 1, -1),
  _parallelBoundaryViscosityReadIterator(_turbulentFlowField, parameters, _viscosityBufferReadStencil, 1, -1) {}

void ParallelManagers::TurbulentPetscParallelManager::communicateObstacleCoordinates(
  std::vector<std::tuple<RealType, RealType>>& coordinatesList2D, std::vector<std::tuple<RealType, RealType, RealType>>& coordinatesList3D
) {
  int max_vector_size_2d;
  int max_vector_size_3d;
  int local_size_2d = coordinatesList2D.size() * 2;
  int local_size_3d = coordinatesList3D.size() * 3;

  MPI_Allreduce(&local_size_2d, &max_vector_size_2d, 1, MPI_INT, MPI_MAX, PETSC_COMM_WORLD);
  MPI_Allreduce(&local_size_3d, &max_vector_size_3d, 1, MPI_INT, MPI_MAX, PETSC_COMM_WORLD);

  int num_ranks = _parameters.parallel.numProcessors[0] * _parameters.parallel.numProcessors[1] * _parameters.parallel.numProcessors[2];

  std::vector<RealType> list2d;
  std::vector<RealType> list3d;

  if (_parameters.parallel.rank == 0) {
    std::vector<RealType> recv_buf_2d;
    recv_buf_2d.resize(max_vector_size_2d);
    std::vector<RealType> recv_buf_3d;
    recv_buf_3d.resize(max_vector_size_3d);

    // flatten coordinate lists
    for (const auto& coordsObst : coordinatesList2D) {
      list2d.push_back(std::get<0>(coordsObst));
      list2d.push_back(std::get<1>(coordsObst));
    }
    for (const auto& coordsObst : coordinatesList3D) {
      list3d.push_back(std::get<0>(coordsObst));
      list3d.push_back(std::get<1>(coordsObst));
      list3d.push_back(std::get<2>(coordsObst));
    }

    for (int i = 1; i < num_ranks; i++) {
      // gather lists from other ranks
      MPI_Recv(&recv_buf_2d[0], max_vector_size_2d, MY_MPI_FLOAT, i, 2, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(&recv_buf_3d[0], max_vector_size_3d, MY_MPI_FLOAT, i, 3, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
      
      // insert items into rank 0 collector
      for (int j = 0; j < max_vector_size_2d; j += 2) {
        RealType first_element = recv_buf_2d[j];
        if (first_element == MY_FLOAT_MAX) {
          break;
        }
        list2d.push_back(first_element);
        list2d.push_back(recv_buf_2d[j + 1]);
      }

      for (int j = 0; j < max_vector_size_3d; j += 3) {
        RealType first_element = recv_buf_3d[j];
        if (first_element == MY_FLOAT_MAX) {
          break;
        }
        list3d.push_back(first_element);
        list3d.push_back(recv_buf_3d[j + 1]);
        list3d.push_back(recv_buf_3d[j + 2]);
      }
    }
    
    int list2d_size = list2d.size();
    int list3d_size = list3d.size();
    
    // communicate new coordinate list sizes to other ranks
    MPI_Bcast(&list2d_size, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&list3d_size, 1, MPI_INT, 0, PETSC_COMM_WORLD);

    // communicate new coordinate lists to other ranks
    MPI_Bcast(&list2d[0], list2d_size, MY_MPI_FLOAT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&list3d[0], list3d_size, MY_MPI_FLOAT, 0, PETSC_COMM_WORLD);
  } else {
    std::vector<RealType> send_buf_2d;
    std::vector<RealType> send_buf_3d;
    
    // flatten own obstacle list
    for (const auto& coordsObst2D : coordinatesList2D) {
      send_buf_2d.push_back(std::get<0>(coordsObst2D));
      send_buf_2d.push_back(std::get<1>(coordsObst2D));
    }
    for (const auto& coordsObst3D : coordinatesList3D) {
      send_buf_3d.push_back(std::get<0>(coordsObst3D));
      send_buf_3d.push_back(std::get<1>(coordsObst3D));
      send_buf_3d.push_back(std::get<2>(coordsObst3D));
    }
    
    // resize to fixed length
    send_buf_2d.resize(max_vector_size_2d, MY_FLOAT_MAX);
    send_buf_3d.resize(max_vector_size_3d, MY_FLOAT_MAX);

    // send flattened obstacle coordinates
    MPI_Send(&send_buf_2d[0], max_vector_size_2d, MY_MPI_FLOAT, 0, 2, PETSC_COMM_WORLD);
    MPI_Send(&send_buf_3d[0], max_vector_size_3d, MY_MPI_FLOAT, 0, 3, PETSC_COMM_WORLD);

    int list2d_size;
    int list3d_size;
    // receive new list sizes from rank 0
    MPI_Bcast(&list2d_size, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&list3d_size, 1, MPI_INT, 0, PETSC_COMM_WORLD);

    // resize to fixed size
    list2d.resize(list2d_size);
    list3d.resize(list3d_size);

    // receive new coordinate lists from rank 0
    MPI_Bcast(&list2d[0], list2d_size, MY_MPI_FLOAT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&list3d[0], list3d_size, MY_MPI_FLOAT, 0, PETSC_COMM_WORLD);
  }

  coordinatesList2D.clear();
  coordinatesList3D.clear();
  
  // reassemble flattened list into list of tuples
  for (long unsigned int i = 0; i < list2d.size(); i += 2) {
    coordinatesList2D.push_back(std::make_tuple(list2d[i], list2d[i + 1]));
  }

  for (long unsigned int i = 0; i < list3d.size(); i += 3) {
    coordinatesList3D.push_back(std::make_tuple(list3d[i], list3d[i + 1], list3d[i + 2]));
  }
}

void ParallelManagers::TurbulentPetscParallelManager::communicateViscosity() {
  int                         counter = 0;
  std::array<MPI_Request, 12> requests{};

  _parallelBoundaryViscosityFillIterator.iterate();

  MPI_Irecv(
    _viscosityBufferReadStencil.getviscosityLeft().data(), _viscosityBufferReadStencil.getviscosityLeft().size(), MY_MPI_FLOAT, _left, 1, PETSC_COMM_WORLD, &requests.at(counter++)
  );
  MPI_Irecv(
    _viscosityBufferReadStencil.getviscosityRight().data(),
    _viscosityBufferReadStencil.getviscosityRight().size(),
    MY_MPI_FLOAT,
    _right,
    1,
    PETSC_COMM_WORLD,
    &requests.at(counter++)
  );
  MPI_Irecv(
    _viscosityBufferReadStencil.getviscosityBottom().data(),
    _viscosityBufferReadStencil.getviscosityBottom().size(),
    MY_MPI_FLOAT,
    _bottom,
    1,
    PETSC_COMM_WORLD,
    &requests.at(counter++)
  );
  MPI_Irecv(
    _viscosityBufferReadStencil.getviscosityTop().data(), _viscosityBufferReadStencil.getviscosityTop().size(), MY_MPI_FLOAT, _top, 1, PETSC_COMM_WORLD, &requests.at(counter++)
  );
  MPI_Irecv(
    _viscosityBufferReadStencil.getviscosityFront().data(),
    _viscosityBufferReadStencil.getviscosityFront().size(),
    MY_MPI_FLOAT,
    _front,
    1,
    PETSC_COMM_WORLD,
    &requests.at(counter++)
  );
  MPI_Irecv(
    _viscosityBufferReadStencil.getviscosityBack().data(), _viscosityBufferReadStencil.getviscosityBack().size(), MY_MPI_FLOAT, _back, 1, PETSC_COMM_WORLD, &requests.at(counter++)
  );

  MPI_Isend(
    _viscosityBufferFillStencil.getviscosityLeft().data(), _viscosityBufferFillStencil.getviscosityLeft().size(), MY_MPI_FLOAT, _left, 1, PETSC_COMM_WORLD, &requests.at(counter++)
  );
  MPI_Isend(
    _viscosityBufferFillStencil.getviscosityRight().data(),
    _viscosityBufferFillStencil.getviscosityRight().size(),
    MY_MPI_FLOAT,
    _right,
    1,
    PETSC_COMM_WORLD,
    &requests.at(counter++)
  );
  MPI_Isend(
    _viscosityBufferFillStencil.getviscosityBottom().data(),
    _viscosityBufferFillStencil.getviscosityBottom().size(),
    MY_MPI_FLOAT,
    _bottom,
    1,
    PETSC_COMM_WORLD,
    &requests.at(counter++)
  );
  MPI_Isend(
    _viscosityBufferFillStencil.getviscosityTop().data(), _viscosityBufferFillStencil.getviscosityTop().size(), MY_MPI_FLOAT, _top, 1, PETSC_COMM_WORLD, &requests.at(counter++)
  );
  MPI_Isend(
    _viscosityBufferFillStencil.getviscosityFront().data(),
    _viscosityBufferFillStencil.getviscosityFront().size(),
    MY_MPI_FLOAT,
    _front,
    1,
    PETSC_COMM_WORLD,
    &requests.at(counter++)
  );
  MPI_Isend(
    _viscosityBufferFillStencil.getviscosityBack().data(), _viscosityBufferFillStencil.getviscosityBack().size(), MY_MPI_FLOAT, _back, 1, PETSC_COMM_WORLD, &requests.at(counter++)
  );

  MPI_Waitall(counter, requests.data(), MPI_STATUSES_IGNORE);

  _parallelBoundaryViscosityReadIterator.iterate();
}
