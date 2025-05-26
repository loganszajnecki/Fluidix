#include "Grid1D.h"

Grid1D::Grid1D(size_t N, double x_start, double x_end)
    : x0(x_start), x1(x_end), dx_((x_end - x_start) / N), cells(N) {}
