#pragma once
#include "Grid1D.h"
#include "GasModel.h"
#include "RusanovFlux.h"

class Solver1D {
public:
    Solver1D(Grid1D& grid, const GasModel& gas);

    void advance(double dt);

private:
    Grid1D& grid_;
    RusanovFlux flux_;
};
