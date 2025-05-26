#include "Solver1D.h"
#include "VTKWriter1D.h"
#include "GnuplotWriter1D.h"
#include <cmath>
#include <iostream>
#include <string>
#include <cstdlib>

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: ./Fluidix <num_cells> <output_file>\n";
        return 1;
    }

    int N = std::atoi(argv[1]);
    std::string output_path = argv[2];

    Grid1D grid(N, 0.0, 1.0);
    GasModel gas(1.4);

    // Initial condition: Sod shock tube
    for (size_t i = 0; i < N; ++i) {
        double x = i * grid.dx();
        if (x < 0.5) {
            grid[i] = ConservedVariables(1.0, 0.0, 2.5);
        } else {
            grid[i] = ConservedVariables(0.125, 0.0, 0.25);
        }
    }

    Solver1D solver(grid, gas);

    const double CFL = 0.5;
    const double final_time = 0.2;
    double time = 0.0;

    while (time < final_time) {
        // Estimate max wave speed for dt
        double max_speed = 0.0;
        for (size_t i = 0; i < grid.num_cells(); ++i) {
            auto prim = gas.conservedToPrimitive(grid[i]);
            double a = gas.soundSpeed(prim);
            max_speed = std::max(max_speed, std::abs(prim.u) + a);
        }

        double dt = CFL * grid.dx() / max_speed;
        if (time + dt > final_time) dt = final_time - time;

        solver.advance(dt);
        time += dt;
    }

    // Output final primitive states
    for (size_t i = 0; i < grid.num_cells(); ++i) {
        auto prim = gas.conservedToPrimitive(grid[i]);
        std::cout << i * grid.dx() << " " << prim.rho << " " << prim.u << " " << prim.p << "\n";
    }

    // VTKWriter1D writer(grid, gas);
    // writer.writePolyData("output.vtk");
    GnuplotWriter1D writerText(grid, gas);
    writerText.writeDataFile(output_path);

    return 0;
}
