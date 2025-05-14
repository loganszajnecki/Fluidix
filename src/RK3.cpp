#include "RK3.h"

RK3::RK3(Solver1D& solver) : solver_(solver) {}

void RK3::advance(Grid1D& grid, double dt) {
    size_t N = grid.num_cells();

     // Stage 0: store U^n
    std::vector<ConservedVariables> U0 = grid.data();
    std::vector<ConservedVariables> U1(N), U2(N), rhs(N);

    // === Stage 1 ===
    solver_.computeRHS(grid, rhs);
    for (size_t i = 0; i < N; ++i)
        U1[i] = U0[i] + dt * rhs[i];

    // === Stage 2 ===
    grid.set(U1);                     // temporary update
    solver_.computeRHS(grid, rhs);   // RHS at U1
    for (size_t i = 0; i < N; ++i)
        U2[i] = 0.75 * U0[i] + 0.25 * (U1[i] + dt * rhs[i]);

    // === Stage 3 ===
    grid.set(U2);
    solver_.computeRHS(grid, rhs);   // RHS at U2
    for (size_t i = 0; i < N; ++i)
        U1[i] = (1.0/3.0) * U0[i] + (2.0/3.0) * (U2[i] + dt * rhs[i]);

    // Final update
    grid.set(U1);
}