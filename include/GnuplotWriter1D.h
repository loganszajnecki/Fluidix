#pragma once
#include <string>
#include "Grid1D.h"
#include "GasModel.h"

class GnuplotWriter1D {
public:
    GnuplotWriter1D(const Grid1D& grid, const GasModel& gas);

    // Writes a .dat file with columns: x, rho, u, p
    void writeDataFile(const std::string& filename) const;

private:
    const Grid1D& grid_;
    const GasModel& gas_;
};
