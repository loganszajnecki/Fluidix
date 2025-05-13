#pragma once
#include <algorithm>

inline double minmod(double a, double b) {
    if (a * b <= 0.0) return 0.0;
    return (std::abs(a) < std::abs(b)) ? a : b;
}

inline double mc(double a, double b) {
    if (a * b <= 0.0) return 0.0;
    return (std::abs(a + b) < 2.0 * std::min(std::abs(a), std::abs(b))) ? 0.5 * (a + b) : 2.0 * std::min(std::abs(a), std::abs(b)) * ((a > 0) ? 1.0 : -1.0);
}

inline double vanLeer(double a, double b) {
    return (a * b > 0.0) ? (2.0 * a * b) / (a + b) : 0.0;
}
