#pragma once

#include "Definitions.hpp"
#include <cmath>

inline RealType ReX(float x, float rho, float U, float mu) {
    return (rho * U * x ) / mu;
}

inline RealType deltaZero([[maybe_unused]] float x, [[maybe_unused]] float rho, [[maybe_unused]] float U, [[maybe_unused]] float mu) {
    return 0.0;
};

inline RealType deltaLaminar(float x, float rho, float U, float mu) {
    return 4.91 * x / sqrt(ReX(x, rho, U, mu));
};

inline RealType deltaTurbulent(float x, float rho, float U, float mu) {
    RealType rex = ReX(x, rho, U, mu);
    return 0.382 * x / (std::pow(rex, 1 / 5));
};