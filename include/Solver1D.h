#pragma once
#include "Grid1D.h"
#include "GasModel.h"
#include "RusanovFlux.h"
#include "HLLCFlux.h"

class Solver1D {
public:
    Solver1D(Grid1D& grid, const GasModel& gas);

    void advance(double dt);

private:
    Grid1D& grid_;
    const GasModel& gas_; 
    // RusanovFlux flux_;
    HLLCFlux flux_;
};
