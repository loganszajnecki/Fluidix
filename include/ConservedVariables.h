#pragma once

struct ConservedVariables {
    double rho;    // density
    double rho_u;  // momentum
    double E;      // total energy

    // Constructors
    ConservedVariables() : rho(0.0), rho_u(0.0), E(0.0) {}
    ConservedVariables(double r, double ru, double e) : rho(r), rho_u(ru), E(e) {}

    // Member operators
    inline ConservedVariables operator+(const ConservedVariables& other) const {
        return {rho + other.rho, rho_u + other.rho_u, E + other.E};
    }

    inline ConservedVariables operator-(const ConservedVariables& other) const {
        return {rho - other.rho, rho_u - other.rho_u, E - other.E};
    }

    inline ConservedVariables operator*(double scalar) const {
        return {rho * scalar, rho_u * scalar, E * scalar};
    }

    inline ConservedVariables operator/(double scalar) const {
        return {rho / scalar, rho_u / scalar, E / scalar};
    }

    inline ConservedVariables operator-() const {
        return {-rho, -rho_u, -E};
    }
};

// Non-member symmetric scalar * vector operator
inline ConservedVariables operator*(double scalar, const ConservedVariables& U) {
    return U * scalar;
}
