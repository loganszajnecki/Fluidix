#include "VTKWriter1D.h"
#include <fstream>
#include <iomanip>

VTKWriter1D::VTKWriter1D(const Grid1D& grid, const GasModel& gas)
    : grid_(grid), gas_(gas) {}

void VTKWriter1D::writePolyData(const std::string& filename) const {
    std::ofstream out(filename);
    if (!out.is_open()) return;

    size_t N = grid_.num_cells();
    double dx = grid_.dx();

    out << "# vtk DataFile Version 3.0\n";
    out << "Hypersonic 1D CFD Output (POLYDATA)\n";
    out << "ASCII\n";
    out << "DATASET POLYDATA\n";

    out << "POINTS " << N << " float\n";
    for (size_t i = 0; i < N; ++i) {
        double x = i * dx;
        out << std::fixed << std::setprecision(6) << x << " 0 0\n";
    }

    out << "VERTICES " << N << " " << 2 * N << "\n";
    for (size_t i = 0; i < N; ++i) {
        out << "1 " << i << "\n";
    }

    out << "POINT_DATA " << N << "\n";

    auto write_scalar = [&](const std::string& name, auto accessor) {
        out << "SCALARS " << name << " float 1\n";
        out << "LOOKUP_TABLE default\n";
        for (size_t i = 0; i < N; ++i) {
            auto prim = gas_.conservedToPrimitive(grid_[i]);
            out << accessor(prim) << "\n";
        }
    };

    write_scalar("density", [](const auto& p) { return p.rho; });
    write_scalar("velocity", [](const auto& p) { return p.u; });
    write_scalar("pressure", [](const auto& p) { return p.p; });

    out.close();
}
