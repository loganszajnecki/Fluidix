#pragma once

#include "ConservedVariables.h"
#include "GasModel.h"

class RusanovFlux {
public:
    RusanovFlux(const GasModel& gas) : gas_(gas) {}

    ConservedVariables compute(const ConservedVariables& UL,
                                const ConservedVariables& UR) const;

private:
    const GasModel& gas_;

    ConservedVariables physicalFlux(const ConservedVariables& U) const;
};
