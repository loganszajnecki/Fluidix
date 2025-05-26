#include "GasModel.h"
#include <cmath>

PrimitiveVariables GasModel::conservedToPrimitive(const ConservedVariables& U) const {
    double rho = U.rho;
    double u = U.rho_u / rho;
    double E = U.E;
    double p = (gamma_ - 1.0) * (E - 0.5 * rho * u * u);
    return {rho, u, p};
}

double GasModel::soundSpeed(const PrimitiveVariables& V) const {
    return std::sqrt(gamma_ * V.p / V.rho);
}