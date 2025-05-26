#pragma once
#include <vector>
#include "ConservedVariables.h"

class MUSCLReconstructor1D {
public:
    // Outputs reconstructed interface states:
    // UL[i] = left state at face i
    // UR[i] = right state at face i
    void reconstruct(const std::vector<ConservedVariables>& U,
                     std::vector<ConservedVariables>& UL,
                     std::vector<ConservedVariables>& UR) const;
};
