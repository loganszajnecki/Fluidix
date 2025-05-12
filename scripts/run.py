import subprocess
import os
import numpy as np
# === Configuration ===
executable = "../build/Fluidix.exe"  # or "./Fluidix" on Linux
output_dir = "../validation/outputs/"
analytical_script = "sod_analytical.py"
grid_sizes = [int(n) for n in np.linspace(50, 3200, 5)]

# === Ensure output directory exists ===
os.makedirs(output_dir, exist_ok=True)

# === Loop over grid sizes ===
for N in grid_sizes:
    print(f"[+] Grid size N = {N}")

    # 1. Generate exact solution at matching resolution
    print(f"    [•] Generating sod_exact.dat with Nx = {N}")
    result = subprocess.run(["python", analytical_script, str(N)],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        print("    [!] Error in sod_analytical.py:")
        print(result.stderr)
        continue

    # 2. Run C++ CFD code
    output_file = os.path.join(output_dir, f"output_N{N}.dat")
    cmd = [executable, str(N), output_file]
    print(f"    [•] Running Fluidix → {output_file}")
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    if result.returncode != 0:
        print(f"    [!] Error running Fluidix N={N}:")
        print(result.stderr)
    else:
        print(f"    [✓] Finished N={N} → {output_file}")
