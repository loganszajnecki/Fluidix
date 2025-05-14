#pragma once
#include "Grid1D.h"
#include "GasModel.h"
#include "RusanovFlux.h"
#include "HLLCFlux.h"

class Solver1D {
public:
    Solver1D(Grid1D& grid, const GasModel& gas);

    void advance(double dt);

    void computeRHS(const Grid1D& grid, std::vector<ConservedVariables>& rhs);


private:
    Grid1D& grid_;
    const GasModel& gas_; 
    // RusanovFlux flux_;
    HLLCFlux flux_;
};
