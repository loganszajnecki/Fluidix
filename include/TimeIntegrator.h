#pragma once
#include "Grid1D.h"

class TimeIntegrator {
public:
    virtual ~TimeIntegrator() = default;

    virtual void advance(Grid1D& grid, double dt) = 0;
};
