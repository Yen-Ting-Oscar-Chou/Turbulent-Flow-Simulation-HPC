#pragma once

#include "Definitions.hpp"
#include <cmath>

inline RealType ReX(const RealType x, const RealType Re) {
    return Re * x;
}

inline RealType deltaZero([[maybe_unused]] const RealType x, [[maybe_unused]] const RealType Re) {
    return 0.0;
}

inline RealType deltaLaminar(const RealType x, const RealType Re) {
    return 4.91 * x / sqrt(ReX(x, Re));
}

inline RealType deltaTurbulent(const RealType x, const RealType Re) {
    const RealType rex = ReX(x, Re);
    return 0.382 * x / (std::pow(rex, 1 / 5));
}