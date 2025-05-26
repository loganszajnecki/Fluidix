#include "RusanovFlux.h"
#include <cmath>

ConservedVariables RusanovFlux::compute(const ConservedVariables& UL, 
                                        const ConservedVariables& UR) const {
    auto VL = gas_.conservedToPrimitive(UL);
    auto VR = gas_.conservedToPrimitive(UR);

    double aL = gas_.soundSpeed(VL);
    double aR = gas_.soundSpeed(VR);

    double s_max = std::max(std::abs(VL.u) + aL, std::abs(VR.u) + aR);

    auto FL = physicalFlux(UL);
    auto FR = physicalFlux(UR);

    return (FL + FR) * 0.5 - (UR - UL) * (0.5 * s_max);
}

ConservedVariables RusanovFlux::physicalFlux(const ConservedVariables& U) const {
    auto V = gas_.conservedToPrimitive(U);
    double rho = V.rho;
    double u = V.u;
    double p = V.p;

    double rho_u = rho * u;
    double E = U.E;

    return {
        rho_u,
        rho_u * u + p,
        u * (E + p)
    };
}
