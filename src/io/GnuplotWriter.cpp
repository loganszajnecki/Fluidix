#include "GnuplotWriter1D.h"
#include <fstream>
#include <iomanip>

GnuplotWriter1D::GnuplotWriter1D(const Grid1D& grid, const GasModel& gas)
    : grid_(grid), gas_(gas) {}

void GnuplotWriter1D::writeDataFile(const std::string& filename) const {
    std::ofstream out(filename);
    if (!out.is_open()) return;

    size_t N = grid_.num_cells();
    double dx = grid_.dx();

    out << "# x rho u p\n";
    for (size_t i = 0; i < N; ++i) {
        double x = i * dx;
        auto prim = gas_.conservedToPrimitive(grid_[i]);
        out << std::fixed << std::setprecision(8)
            << x << " "
            << prim.rho << " "
            << prim.u   << " "
            << prim.p   << "\n";
    }

    out.close();
}
