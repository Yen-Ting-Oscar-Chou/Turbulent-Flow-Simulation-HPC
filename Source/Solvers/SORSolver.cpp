#include "StdAfx.hpp"

#include "SORSolver.hpp"

void Solvers::SORSolver::solve(FlowField& flowField, const Parameters& parameters) {
  RealType resnorm = DBL_MAX, tol = 1e-4;

  double omg        = 1.7;
  int    iterations = -1;
  int    it         = 0;

  int          nx = flowField.getNx(), ny = flowField.getNy(), nz = flowField.getNz();
  ScalarField& P = flowField.getPressure();
  if (parameters.geometry.dim == 3) {
    do {
#pragma omp target teams distribute parallel for collapse(3) num_teams(NUM_TEAMS) num_threads(NUM_THREADS)
      for (int k = 2; k < nz + 2; k++) {
        for (int j = 2; j < ny + 2; j++) {
          for (int i = 2; i < nx + 2; i++) {
            if ((i + j) % 2 == 0) {
              SOR3D(flowField, parameters, omg, i, j, k);
            }
          }
        }
      }

#pragma omp target teams distribute parallel for collapse(3) num_teams(NUM_TEAMS) num_threads(NUM_THREADS)
      for (int k = 2; k < nz + 2; k++) {
        for (int j = 2; j < ny + 2; j++) {
          for (int i = 2; i < nx + 2; i++) {
            if ((i + j) % 2 == 1) {
              SOR3D(flowField, parameters, omg, i, j, k);
            }
          }
        }
      }

#pragma omp target teams distribute parallel for collapse(2) num_teams(NUM_TEAMS) num_threads(NUM_THREADS)
      for (int j = 2; j < ny + 2; j++) {
        for (int k = 2; k < nz + 2; k++) {
          P.getScalar(1, j, k)      = P.getScalar(2, j, k);
          P.getScalar(nx + 2, j, k) = P.getScalar(nx + 1, j, k);
        }
      }

#pragma omp target teams distribute parallel for collapse(2) num_teams(NUM_TEAMS) num_threads(NUM_THREADS)
      for (int i = 2; i < nx + 2; i++) {
        for (int k = 2; k < nz + 2; k++) {
          P.getScalar(i, 1, k)      = P.getScalar(i, 2, k);
          P.getScalar(i, ny + 2, k) = P.getScalar(i, ny + 1, k);
        }
      }

#pragma omp target teams distribute parallel for collapse(2) num_teams(NUM_TEAMS) num_threads(NUM_THREADS)
      for (int i = 2; i < nx + 2; i++) {
        for (int j = 2; j < ny + 2; j++) {
          P.getScalar(i, j, 1)      = P.getScalar(i, j, 2);
          P.getScalar(i, j, nz + 2) = P.getScalar(i, j, nz + 1);
        }
      }

      resnorm = 0;
#pragma omp target teams distribute parallel for collapse(3) num_teams(NUM_TEAMS) num_threads(NUM_THREADS) reduction(+ : resnorm)
      for (int k = 2; k < nz + 2; k++) {
        for (int j = 2; j < ny + 2; j++) {
          for (int i = 2; i < nx + 2; i++) {
            const RealType dx_0  = parameters.meshsize->getDx(i, j, k);
            const RealType dx_M1 = parameters.meshsize->getDx(i - 1, j, k);
            const RealType dx_P1 = parameters.meshsize->getDx(i + 1, j, k);
            const RealType dy_0  = parameters.meshsize->getDy(i, j, k);
            const RealType dy_M1 = parameters.meshsize->getDy(i, j - 1, k);
            const RealType dy_P1 = parameters.meshsize->getDy(i, j + 1, k);
            const RealType dz_0  = parameters.meshsize->getDz(i, j, k);
            const RealType dz_M1 = parameters.meshsize->getDz(i, j, k - 1);
            const RealType dz_P1 = parameters.meshsize->getDz(i, j, k + 1);

            const RealType dx_W = 0.5 * (dx_0 + dx_M1);
            const RealType dx_E = 0.5 * (dx_0 + dx_P1);
            const RealType dx_S = 0.5 * (dy_0 + dy_M1);
            const RealType dx_N = 0.5 * (dy_0 + dy_P1);
            const RealType dx_B = 0.5 * (dz_0 + dz_M1);
            const RealType dx_T = 0.5 * (dz_0 + dz_P1);

            const RealType a_W = 2.0 / (dx_W * (dx_W + dx_E));
            const RealType a_E = 2.0 / (dx_E * (dx_W + dx_E));
            const RealType a_N = 2.0 / (dx_N * (dx_N + dx_S));
            const RealType a_S = 2.0 / (dx_S * (dx_N + dx_S));
            const RealType a_T = 2.0 / (dx_T * (dx_T + dx_B));
            const RealType a_B = 2.0 / (dx_B * (dx_T + dx_B));
            const RealType a_C = -2.0 / (dx_E * dx_W) - 2.0 / (dx_N * dx_S) - 2.0 / (dx_B * dx_T);

            resnorm += pow(
              (flowField.getRHS().getScalar(i, j, k) - a_W * P.getScalar(i - 1, j, k) - a_E * P.getScalar(i + 1, j, k) - a_S * P.getScalar(i, j - 1, k)
               - a_N * P.getScalar(i, j + 1, k) - a_B * P.getScalar(i, j, k - 1) - a_T * P.getScalar(i, j, k + 1) - a_C * P.getScalar(i, j, k)),
              2
            );
          }
        }
      }
      resnorm = sqrt(resnorm / (nx * ny * nz));
#ifndef NDEBUG
      spdlog::debug("Residual norm : {}", resnorm);
#endif

      it++;
      iterations--;
    } while (resnorm > tol && iterations);
  }
  if (parameters.geometry.dim == 2) {
    do {
#pragma omp target teams distribute parallel for collapse(2) num_teams(NUM_TEAMS) num_threads(NUM_THREADS)
      for (int j = 2; j < ny + 2; j++) {
        for (int i = 2; i < nx + 2; i++) {
          if ((i + j) % 2 == 0) {
            SOR2D(flowField, parameters, omg, i, j);
          }
        }
      }

#pragma omp target teams distribute parallel for collapse(2) num_teams(NUM_TEAMS) num_threads(NUM_THREADS)
      for (int j = 2; j < ny + 2; j++) {
        for (int i = 2; i < nx + 2; i++) {
          if ((i + j) % 2 == 1) {
            SOR2D(flowField, parameters, omg, i, j);
          }
        }
      }

      resnorm = 0.0;
#pragma omp target teams distribute parallel for collapse(2) num_teams(NUM_TEAMS) num_threads(NUM_THREADS) reduction(+ : resnorm)
      for (int j = 2; j < ny + 2; j++) {
        for (int i = 2; i < nx + 2; i++) {
          const RealType dx_0  = parameters.meshsize->getDx(i, j);
          const RealType dx_M1 = parameters.meshsize->getDx(i - 1, j);
          const RealType dx_P1 = parameters.meshsize->getDx(i + 1, j);
          const RealType dy_0  = parameters.meshsize->getDy(i, j);
          const RealType dy_M1 = parameters.meshsize->getDy(i, j - 1);
          const RealType dy_P1 = parameters.meshsize->getDy(i, j + 1);

          const RealType dx_W = 0.5 * (dx_0 + dx_M1);
          const RealType dx_E = 0.5 * (dx_0 + dx_P1);
          const RealType dx_S = 0.5 * (dy_0 + dy_M1);
          const RealType dx_N = 0.5 * (dy_0 + dy_P1);

          const RealType a_W      = 2.0 / (dx_W * (dx_W + dx_E));
          const RealType a_E      = 2.0 / (dx_E * (dx_W + dx_E));
          const RealType a_N      = 2.0 / (dx_N * (dx_N + dx_S));
          const RealType a_S      = 2.0 / (dx_S * (dx_N + dx_S));
          const RealType a_C      = -2.0 / (dx_E * dx_W) - 2.0 / (dx_N * dx_S);
          const RealType residual = flowField.getRHS().getScalar(i, j) - a_W * P.getScalar(i - 1, j) - a_E * P.getScalar(i + 1, j) - a_S * P.getScalar(i, j - 1)
                                    - a_N * P.getScalar(i, j + 1) - a_C * P.getScalar(i, j);
          resnorm += residual * residual;
        }
      }
      resnorm = sqrt(resnorm / (nx * ny));

#pragma omp target teams distribute parallel for num_teams(NUM_TEAMS) num_threads(NUM_THREADS)
      for (int j = 2; j < ny + 2; j++) {
        P.getScalar(1, j)      = P.getScalar(2, j);
        P.getScalar(nx + 2, j) = P.getScalar(nx + 1, j);
      }

#pragma omp target teams distribute parallel for num_teams(NUM_TEAMS) num_threads(NUM_THREADS)
      for (int i = 2; i < nx + 2; i++) {
        P.getScalar(i, 1)      = P.getScalar(i, 2);
        P.getScalar(i, ny + 2) = P.getScalar(i, ny + 1);
      }

#ifndef NDEBUG
      spdlog::debug("Residual norm : {}", resnorm);
#endif
      iterations--;
      it++;
    } while (resnorm > tol && iterations);
  }

#ifndef NDEBUG
  spdlog::debug("SORSolver needed {} iterations", it);
#endif
}
