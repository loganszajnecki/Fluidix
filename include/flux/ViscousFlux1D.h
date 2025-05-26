#pragma once
#include "GasModel.h"
#include "ConservedVariables.h"
#include <vector>

class ViscousFlux1D {
public:
    ViscousFlux1D(const GasModel& gas, double mu, double k)
        : gas_(gas), mu_(mu), k_(k) {}

    // Compute viscous fluxes for all interior points
    void compute(const std::vector<PrimitiveVariables>& V,
                 std::vector<ConservedVariables>& viscous_flux,
                 double dx) const;

private:
    const GasModel& gas_;
    double mu_;  // dynamic viscosity
    double k_;   // thermal conductivity
};
