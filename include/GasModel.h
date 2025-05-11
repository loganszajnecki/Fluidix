#pragma once
#include "ConservedVariables.h"

struct PrimitiveVariables {
    double rho;
    double u;
    double p;
};

class GasModel {
public:
    GasModel(double gamma) : gamma_(gamma) {}

    PrimitiveVariables conservedToPrimitive(const ConservedVariables& U) const;
    double soundSpeed(const PrimitiveVariables& V) const;

private:
    double gamma_;
};
