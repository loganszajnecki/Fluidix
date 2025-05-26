#pragma once 
#include "GasModel.h"
#include "ConservedVariables.h"
#include <cmath>

class HLLCFlux {
public:
    HLLCFlux(const GasModel& gas) : gas_(gas) {}

    ConservedVariables compute(const ConservedVariables& UL, const ConservedVariables& UR) const;

private:
    const GasModel& gas_;
    ConservedVariables physicalFlux(const ConservedVariables& U) const;
};