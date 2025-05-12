import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import os
import re

# === Load exact solution ===
exact_data = np.loadtxt("../validation/sod_exact.dat", skiprows=1)
x_exact = exact_data[:, 0]
rho_exact = exact_data[:, 1]
u_exact   = exact_data[:, 2]
p_exact   = exact_data[:, 3]

# === Discover available grid sizes dynamically ===
output_dir = "../validation/outputs"
pattern = re.compile(r"output_N(\d+)\.dat")

grid_files = []
for fname in os.listdir(output_dir):
    match = pattern.match(fname)
    if match:
        N = int(match.group(1))
        grid_files.append((N, os.path.join(output_dir, fname)))

# Sort by grid size
grid_files.sort(key=lambda x: x[0])
grid_sizes = [entry[0] for entry in grid_files]

# === Compute L2 error for each file ===
def compute_l2_error(filepath):
    numerical = np.loadtxt(filepath, skiprows=1)
    x_num = numerical[:, 0]
    rho_num = numerical[:, 1]
    u_num   = numerical[:, 2]
    p_num   = numerical[:, 3]

    interp_rho = interp1d(x_exact, rho_exact, kind='linear', fill_value='extrapolate')
    interp_u   = interp1d(x_exact, u_exact,   kind='linear', fill_value='extrapolate')
    interp_p   = interp1d(x_exact, p_exact,   kind='linear', fill_value='extrapolate')

    rho_ref = interp_rho(x_num)
    u_ref   = interp_u(x_num)
    p_ref   = interp_p(x_num)

    l2_rho = np.sqrt(np.mean((rho_num - rho_ref)**2))
    l2_u   = np.sqrt(np.mean((u_num - u_ref)**2))
    l2_p   = np.sqrt(np.mean((p_num - p_ref)**2))

    return l2_rho, l2_u, l2_p

# === Loop and compute errors ===
l2_rho_all, l2_u_all, l2_p_all = [], [], []

for N, path in grid_files:
    l2_rho, l2_u, l2_p = compute_l2_error(path)
    l2_rho_all.append(l2_rho)
    l2_u_all.append(l2_u)
    l2_p_all.append(l2_p)
    print(f"N = {N:4d} | L2 ρ = {l2_rho:.3e} | L2 u = {l2_u:.3e} | L2 p = {l2_p:.3e}")

# === Plotting ===
plt.figure(figsize=(6, 4))
plt.loglog(grid_sizes, l2_rho_all, 'o-', label='Density')
plt.loglog(grid_sizes, l2_u_all,   's--', label='Velocity')
plt.loglog(grid_sizes, l2_p_all,   'd-.', label='Pressure')

plt.grid(True, which='both', ls='--')
plt.xlabel("Grid size (N)")
plt.ylabel("L2 Error")
plt.title("Grid Convergence Study")
plt.legend()

# Show slope for density
if len(l2_rho_all) >= 2:
    order = np.log(l2_rho_all[-2] / l2_rho_all[-1]) / np.log(grid_sizes[-2] / grid_sizes[-1])
    plt.text(grid_sizes[-1], l2_rho_all[-1] * 1.5, f"ρ slope ≈ {order:.2f}", fontsize=9)

plt.tight_layout()
plt.savefig("../validation/convergence_plot.png", dpi=150)
plt.show()
