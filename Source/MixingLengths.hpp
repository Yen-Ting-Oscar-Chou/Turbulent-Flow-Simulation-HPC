#pragma once

#include "Definitions.hpp"
#include <cmath>

inline RealType ReX(RealType x, RealType Re) {
    return Re * x;
}

inline RealType deltaZero([[maybe_unused]] RealType x, [[maybe_unused]] RealType Re) {
    return 0.0;
};

inline RealType deltaLaminar(RealType x, RealType Re) {
    return 4.91 * x / sqrt(ReX(x, Re));
};

inline RealType deltaTurbulent(RealType x, RealType Re) {
    const RealType rex = ReX(x, Re);
    return 0.382 * x / (std::pow(rex, 1 / 5));
};