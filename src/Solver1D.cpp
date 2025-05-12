#include "Solver1D.h"
#include "MUSCLReconstructor1D.h"

Solver1D::Solver1D(Grid1D& grid, const GasModel& gas)
    : grid_(grid), flux_(gas) {}

void Solver1D::advance(double dt) {
    size_t N = grid_.num_cells();
    double dx = grid_.dx();

    std::vector<ConservedVariables> fluxes(N + 1); // N+1 interfaces

    // Compute numerical fluxes at interfaces, Rusanov Flux 
    // for (size_t i = 0; i < N - 1; ++i) {
    //     fluxes[i + 1] = flux_.compute(grid_[i], grid_[i + 1]);
    // }

    MUSCLReconstructor1D reconstructor;
    std::vector<ConservedVariables> UL, UR;
    reconstructor.reconstruct(grid_.data(), UL, UR);

    for (size_t i = 0; i < UL.size(); ++i) {
        fluxes[i + 1] = flux_.compute(UL[i], UR[i]);
    }

    // Update conserved variables using FV update
    std::vector<ConservedVariables> new_state(N);
    for (size_t i = 1; i < N - 1; ++i) {
        auto dF = fluxes[i + 1] - fluxes[i];
        new_state[i] = grid_[i] - dF * (dt / dx);
    }

    // Copy new state into the grid
    for (size_t i = 1; i < N - 1; ++i) {
        grid_[i] = new_state[i];
    }
}
