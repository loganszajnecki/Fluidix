from scipy.interpolate import interp1d
import numpy as np

# Load data
numerical = np.loadtxt("../build/output.dat", comments="#", skiprows=1)
exact = np.loadtxt("../validation/sod_exact.dat", comments="#", skiprows=1)

x_num = numerical[:, 0]
rho_num = numerical[:, 1]
u_num   = numerical[:, 2]
p_num   = numerical[:, 3]

x_exact = exact[:, 0]
rho_exact = exact[:, 1]
u_exact   = exact[:, 2]
p_exact   = exact[:, 3]

# Interpolate exact values onto numerical grid
interp_rho = interp1d(x_exact, rho_exact, kind='linear', fill_value="extrapolate")
interp_u   = interp1d(x_exact, u_exact,   kind='linear', fill_value="extrapolate")
interp_p   = interp1d(x_exact, p_exact,   kind='linear', fill_value="extrapolate")

rho_ref = interp_rho(x_num)
u_ref   = interp_u(x_num)
p_ref   = interp_p(x_num)

# Compute L2 norms
l2_rho = np.sqrt(np.mean((rho_num - rho_ref)**2))
l2_u   = np.sqrt(np.mean((u_num - u_ref)**2))
l2_p   = np.sqrt(np.mean((p_num - p_ref)**2))

# Print results
print(f"L2 error in density:   {l2_rho:.6e}")
print(f"L2 error in velocity:  {l2_u:.6e}")
print(f"L2 error in pressure:  {l2_p:.6e}")
