#include "HLLCFlux.h"
#include <algorithm>

ConservedVariables HLLCFlux::compute(const ConservedVariables& UL,
                                     const ConservedVariables& UR) const {
    auto VL = gas_.conservedToPrimitive(UL);
    auto VR = gas_.conservedToPrimitive(UR);

    double rhoL = VL.rho, uL = VL.u, pL = VL.p;
    double rhoR = VR.rho, uR = VR.u, pR = VR.p;

    double aL = gas_.soundSpeed(VL);
    double aR = gas_.soundSpeed(VR);

    // --- Estimate wave speeds S_L and S_R ---
    double SL = std::min(uL - aL, uR - aR);
    double SR = std::max(uL + aL, uR + aR);

    // --- Compute pressure in star region using Toro's average (optional)
    double pPVRS = 0.5 * (pL + pR) - 0.5 * (uR - uL) * 0.25 * (rhoL + rhoR) * (aL + aR);
    pPVRS = std::max(0.0, pPVRS);

    // --- Compute contact wave speed S_* ---
    double numerator = pR - pL + rhoL * uL * (SL - uL) - rhoR * uR * (SR - uR);
    double denominator = rhoL * (SL - uL) - rhoR * (SR - uR);
    double S_star = numerator / denominator;

    // --- Compute physical fluxes ---
    auto FL = physicalFlux(UL);
    auto FR = physicalFlux(UR);

    // --- Star region fluxes (F*L and F*R) ---
    double rho_star_L = rhoL * (SL - uL) / (SL - S_star);
    double E_star_L = UL.E / rhoL + (S_star - uL) * (S_star + pL / (rhoL * (SL - uL)));
    ConservedVariables U_star_L(
        rho_star_L,
        rho_star_L * S_star,
        rho_star_L * E_star_L
    );

    double rho_star_R = rhoR * (SR - uR) / (SR - S_star);
    double E_star_R = UR.E / rhoR + (S_star - uR) * (S_star + pR / (rhoR * (SR - uR)));
    ConservedVariables U_star_R(
        rho_star_R,
        rho_star_R * S_star,
        rho_star_R * E_star_R
    );

    if (0 <= SL) {
        return FL;
    } else if (SL <= 0 && 0 <= S_star) {
        return FL + (U_star_L - UL) * SL;
    } else if (S_star <= 0 && 0 <= SR) {
        return FR + (U_star_R - UR) * SR;
    } else {
        return FR;
    }
}

ConservedVariables HLLCFlux::physicalFlux(const ConservedVariables& U) const {
    auto V = gas_.conservedToPrimitive(U);
    double rho = V.rho;
    double u = V.u;
    double p = V.p;

    double rho_u = rho * u;
    double E = U.E;

    return {
        rho_u,
        rho_u * u + p,
        u * (E + p)
    };
}
