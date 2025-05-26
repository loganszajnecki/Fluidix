#include "ConservedVariables.h"

ConservedVariables ConservedVariables::operator+(const ConservedVariables& other) const {
    return {rho + other.rho, rho_u + other.rho_u, E + other.E};
}
ConservedVariables ConservedVariables::operator-(const ConservedVariables& other) const {
    return {rho - other.rho, rho_u - other.rho_u, E - other.E};
}
ConservedVariables ConservedVariables::operator*(double scalar) const {
    return {rho * scalar, rho_u * scalar, E * scalar};
}
