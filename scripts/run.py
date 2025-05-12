import subprocess
import os

# === Configuration ===
executable = "../build/Fluidix.exe"  # or "Fluidix.exe" on Windows
output_dir = "../validation/outputs/"
grid_sizes = [50, 100, 200, 300, 400, 600, 800, 1600, 3200]

# === Ensure output directory exists ===
os.makedirs(output_dir, exist_ok=True)

# === Loop over grid sizes ===
for N in grid_sizes:
    output_file = os.path.join(output_dir, f"output_N{N}.dat")
    cmd = [executable, str(N), output_file]

    print(f"[+] Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    if result.returncode != 0:
        print(f"    [!] Error running N={N}:")
        print(result.stderr)
    else:
        print(f"    [✓] Finished N={N} → {output_file}")
