import subprocess
import os
import numpy as np
# === Configuration ===
executable = "../build/Fluidix.exe"  # or "./Fluidix" on Linux
output_dir = "../validation/outputs/"
analytical_script = "sod_analytical.py"
#grid_sizes = [int(n) for n in np.linspace(10, 500, 5)]
grid_sizes = [25, 50, 100, 200, 400, 800, 1600, 3200, 6400]

# === Ensure output directory exists ===
os.makedirs(output_dir, exist_ok=True)

# === Loop over grid sizes ===
for N in grid_sizes:
    print(f"[+] Grid size N = {N}")

    # 1. Generate exact solution at matching resolution
    exact_path = os.path.join(output_dir, f"sod_exact_N{N}.dat")
    print(f"    [•] Generating {exact_path} with Nx = {N}")
    result = subprocess.run(["py", analytical_script, str(N), exact_path],
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
