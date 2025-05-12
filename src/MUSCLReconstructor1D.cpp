#include "MUSCLReconstructor1D.h"
#include "SlopeLimiter.h"

void MUSCLReconstructor1D::reconstruct(const std::vector<ConservedVariables>& U,
                                       std::vector<ConservedVariables>& UL,
                                       std::vector<ConservedVariables>& UR) const {
    size_t N = U.size();
    UL.resize(N - 1);
    UR.resize(N - 1);

    for (size_t i = 1; i < N - 1; ++i) {
        // Slopes
        auto duL = U[i]   - U[i - 1];
        auto duR = U[i+1] - U[i];

        // Apply minmod limiter component-wise
        ConservedVariables slope(
            minmod(duL.rho,   duR.rho),
            minmod(duL.rho_u, duR.rho_u),
            minmod(duL.E,     duR.E)
        );

        // Reconstruct left/right states at face between i and i+1
        UR[i - 1] = U[i]   - slope * 0.5;
        UL[i]     = U[i]   + slope * 0.5;
    }

    // Edge extrapolation
    UL[0]       = U[0];
    UR[N - 2]   = U[N - 1];
}
