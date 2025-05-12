import numpy as np
import matplotlib.pyplot as plt
import os

# === Paths ===
exact_path = "../validation/sod_exact.dat"
numerical_path = "../validation/outputs/output_N3200.dat"
output_dir = "../validation"
os.makedirs(output_dir, exist_ok=True)

# === Load data ===
exact = np.loadtxt(exact_path, skiprows=1)
numerical = np.loadtxt(numerical_path, skiprows=1)

x_exact = exact[:, 0]
rho_exact = exact[:, 1]
u_exact   = exact[:, 2]
p_exact   = exact[:, 3]

x_num = numerical[:, 0]
rho_num = numerical[:, 1]
u_num   = numerical[:, 2]
p_num   = numerical[:, 3]

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
    plt.savefig(os.path.join(output_dir, fname), dpi=150)
    plt.show()

# === Plot each variable ===
plot_variable(x_exact, rho_exact, x_num, rho_num, "blue", "Density", "Density", "sod_density.png")
plot_variable(x_exact, u_exact,   x_num, u_num,   "red",  "Velocity", "Velocity", "sod_velocity.png")
plot_variable(x_exact, p_exact,   x_num, p_num,   "green", "Pressure", "Pressure", "sod_pressure.png")
