#pragma once
#include <vector>
#include "ConservedVariables.h"

class Grid1D {
public:
    Grid1D(size_t N, double x_start, double x_end);

    size_t num_cells() const { return cells.size(); }
    double dx() const { return dx_; }

    ConservedVariables& operator[](size_t i) { return cells[i]; }
    const ConservedVariables& operator[](size_t i) const { return cells[i]; }
    const std::vector<ConservedVariables>& data() const { return cells; }

private:
    std::vector<ConservedVariables> cells;
    double x0, x1, dx_;
};