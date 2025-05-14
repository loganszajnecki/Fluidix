#pragma once
#include "TimeIntegrator.h"
#include "Solver1D.h"

class RK3 : public TimeIntegrator {
public:
    RK3(Solver1D& solver);

    void advance(Grid1D& grid, double dt) override;
private:
    Solver1D& solver_;
};