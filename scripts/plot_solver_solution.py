import numpy as np
import matplotlib.pyplot as plt
import os
import re

# === Paths ===
output_dir = "../validation/solution_plots"
input_dir = "../validation/outputs"
os.makedirs(output_dir, exist_ok=True)

# === File matching pattern ===
pattern = re.compile(r"output_N(\d+)\.dat")

# === Common plotting function ===
def plot_variable(x_exact, y_exact, x_num, y_num, color, label, title, fname):
    plt.figure(figsize=(8, 4))
    plt.plot(x_exact, y_exact, color=color, label=f"Exact {label}")
    plt.plot(x_num, y_num, linestyle='none', marker='o', markersize=4,
             color=color, label=f"Numerical {label}")
    plt.xlabel("x")
    plt.ylabel(label)
    plt.title(f"Sod Shock Tube – {title}")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(fname, dpi=150)
    plt.close()

# === Loop over matching output files ===
for fname in sorted(os.listdir(input_dir)):
    match = pattern.match(fname)
    if not match:
        continue

    N = match.group(1)
    numerical_path = os.path.join(input_dir, f"output_N{N}.dat")
    exact_path     = os.path.join(input_dir, f"sod_exact_N{N}.dat")

    if not os.path.exists(exact_path):
        print(f"[!] Skipping N={N}: missing exact solution file.")
        continue

    print(f"[+] Plotting for N = {N}")

    # Load data
    numerical = np.loadtxt(numerical_path, skiprows=1)
    exact = np.loadtxt(exact_path, skiprows=1)

    x_num = numerical[:, 0]
    rho_num = numerical[:, 1]
    u_num   = numerical[:, 2]
    p_num   = numerical[:, 3]

    x_exact = exact[:, 0]
    rho_exact = exact[:, 1]
    u_exact   = exact[:, 2]
    p_exact   = exact[:, 3]

    # Plot
    plot_variable(x_exact, rho_exact, x_num, rho_num, "blue", "Density",
                  "Density", os.path.join(output_dir, f"sod_density_N{N}.png"))
    plot_variable(x_exact, u_exact, x_num, u_num, "red", "Velocity",
                  "Velocity", os.path.join(output_dir, f"sod_velocity_N{N}.png"))
    plot_variable(x_exact, p_exact, x_num, p_num, "green", "Pressure",
                  "Pressure", os.path.join(output_dir, f"sod_pressure_N{N}.png"))
