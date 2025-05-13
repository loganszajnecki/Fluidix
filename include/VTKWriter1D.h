#pragma once
#include <string>
#include "Grid1D.h"
#include "GasModel.h"

// Writes VTK POLYDATA format for 1D CFD results (points along x-axis)
class VTKWriter1D {
public:
    VTKWriter1D(const Grid1D& grid, const GasModel& gas);

    // Writes to a .vtk file compatible with ParaView
    void writePolyData(const std::string& filename) const;

private:
    const Grid1D& grid_;
    const GasModel& gas_;
};
