// #pragma once

// #include "Definitions.hpp"
// #include "FlowField.hpp"
// #include "Parameters.hpp"

// namespace Solvers {

//   class SORSolver {
//   public:
//     SORSolver()  = default;
//     ~SORSolver() = default;

// #pragma omp declare target
//     void        solve(FlowField& flowField, const Parameters& parameters);
//     inline void SOR2D(FlowField& flowField, const Parameters& parameters, double omg, int i, int j) {
//       ScalarField&   P     = flowField.getPressure();
//       const RealType dx_0  = parameters.meshsize->getDx(i, j);
//       const RealType dx_M1 = parameters.meshsize->getDx(i - 1, j);
//       const RealType dx_P1 = parameters.meshsize->getDx(i + 1, j);
//       const RealType dy_0  = parameters.meshsize->getDy(i, j);
//       const RealType dy_M1 = parameters.meshsize->getDy(i, j - 1);
//       const RealType dy_P1 = parameters.meshsize->getDy(i, j + 1);

//       const RealType dx_W = 0.5 * (dx_0 + dx_M1);
//       const RealType dx_E = 0.5 * (dx_0 + dx_P1);
//       const RealType dx_S = 0.5 * (dy_0 + dy_M1);
//       const RealType dx_N = 0.5 * (dy_0 + dy_P1);

//       const RealType a_W = 2.0 / (dx_W * (dx_W + dx_E));
//       const RealType a_E = 2.0 / (dx_E * (dx_W + dx_E));
//       const RealType a_N = 2.0 / (dx_N * (dx_N + dx_S));
//       const RealType a_S = 2.0 / (dx_S * (dx_N + dx_S));
//       const RealType a_C = -2.0 / (dx_E * dx_W) - 2.0 / (dx_N * dx_S);

//       const RealType gaussSeidel
//         = 1.0 / a_C * (flowField.getRHS().getScalar(i, j) - a_W * P.getScalar(i - 1, j) - a_E * P.getScalar(i + 1, j) - a_S * P.getScalar(i, j - 1) - a_N * P.getScalar(i, j + 1));
//       P.getScalar(i, j) = omg * gaussSeidel + (1.0 - omg) * P.getScalar(i, j);
//     }

//     inline void SOR3D(FlowField& flowField, const Parameters& parameters, double omg, int i, int j, int k) {
//       ScalarField& P = flowField.getPressure();
//       const RealType dx_0  = parameters.meshsize->getDx(i, j, k);
//       const RealType dx_M1 = parameters.meshsize->getDx(i - 1, j, k);
//       const RealType dx_P1 = parameters.meshsize->getDx(i + 1, j, k);
//       const RealType dy_0  = parameters.meshsize->getDy(i, j, k);
//       const RealType dy_M1 = parameters.meshsize->getDy(i, j - 1, k);
//       const RealType dy_P1 = parameters.meshsize->getDy(i, j + 1, k);
//       const RealType dz_0  = parameters.meshsize->getDz(i, j, k);
//       const RealType dz_M1 = parameters.meshsize->getDz(i, j, k - 1);
//       const RealType dz_P1 = parameters.meshsize->getDz(i, j, k + 1);

//       const RealType dx_W = 0.5 * (dx_0 + dx_M1);
//       const RealType dx_E = 0.5 * (dx_0 + dx_P1);
//       const RealType dx_S = 0.5 * (dy_0 + dy_M1);
//       const RealType dx_N = 0.5 * (dy_0 + dy_P1);
//       const RealType dx_B = 0.5 * (dz_0 + dz_M1);
//       const RealType dx_T = 0.5 * (dz_0 + dz_P1);

//       const RealType a_W = 2.0 / (dx_W * (dx_W + dx_E));
//       const RealType a_E = 2.0 / (dx_E * (dx_W + dx_E));
//       const RealType a_N = 2.0 / (dx_N * (dx_N + dx_S));
//       const RealType a_S = 2.0 / (dx_S * (dx_N + dx_S));
//       const RealType a_T = 2.0 / (dx_T * (dx_T + dx_B));
//       const RealType a_B = 2.0 / (dx_B * (dx_T + dx_B));
//       const RealType a_C = -2.0 / (dx_E * dx_W) - 2.0 / (dx_N * dx_S) - 2.0 / (dx_B * dx_T);
//       P.getScalar(
//         i, j, k
//       ) = omg / a_C
//             * (flowField.getRHS().getScalar(i, j, k) - a_W * P.getScalar(i - 1, j, k) - a_E * P.getScalar(i + 1, j, k) - a_S * P.getScalar(i, j - 1, k) - a_N * P.getScalar(i, j + 1, k) - a_B * P.getScalar(i, j, k - 1) - a_T * P.getScalar(i, j, k + 1))
//           + (1.0 - omg) * P.getScalar(i, j, k);
//     }
// #pragma omp end declare target

//     inline void reInitMatrix() {};
//   };

// } // namespace Solvers
