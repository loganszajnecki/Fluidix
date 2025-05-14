#include "Grid1D.h"
#include <stdexcept>

Grid1D::Grid1D(size_t N, double x_start, double x_end)
    : x0(x_start), x1(x_end), dx_((x_end - x_start) / N), cells(N) {}

void Grid1D::set(const std::vector<ConservedVariables>& new_cells) {
    if (new_cells.size() != cells.size()) {
        throw std::runtime_error("Grid1D::set() - size mismatch");
    }
    cells = new_cells;
}
