#pragma once

struct ConservedVariables {
    double rho;   // density , assume kg/m^3 for now
    double rho_u; // momentum, SI units for now
    double E;     // energy, SI units for now

    ConservedVariables() : rho(0.0), rho_u(0.0), E(0.0) {}
    ConservedVariables(double r, double ru, double e) : rho(r), rho_u(ru), E(e) {}

    ConservedVariables operator+(const ConservedVariables& other) const;
    ConservedVariables operator-(const ConservedVariables& other) const;
    ConservedVariables operator*(double scalar) const;
};