# Set interactive terminal (auto-detect or set qt/wxt manually)
set terminal qt
set xlabel "x"
set ylabel "Value"
set key outside

# --- Plot interactively ---
plot "../validation/sod_exact.dat" using 1:2 with lines lc rgb "blue" title "Exact Density", \
     "../build/output.dat" using 1:2 with points pt 7 lc rgb "blue" title "Numerical Density", \
     "../validation/sod_exact.dat" using 1:3 with lines lc rgb "red" title "Exact u", \
     "../build/output.dat" using 1:3 with points pt 7 lc rgb "red" title "Numerical Velocity", \
     "../validation/sod_exact.dat" using 1:4 with lines lc rgb "green" title "Exact p", \
     "../build/output.dat" using 1:4 with points pt 7 lc rgb "green" title "Numerical Pressure"


pause -1 "Press Enter to continue to save plots..."

# --- Output files as PNG to validation/ folder ---
set terminal pngcairo size 1000,600 enhanced font 'Arial,10'

# --- Density ---
set output "../validation/density.png"
set title "Sod Shock Tube – Density"
plot "../validation/sod_exact.dat" using 1:2 with lines lc rgb "blue" title "Exact Density", \
     "../build/output.dat" using 1:2 with points pt 7 lc rgb "blue" title "Numerical Density"

# --- Velocity ---
set output "../validation/velocity.png"
set title "Sod Shock Tube – Velocity"
plot "../validation/sod_exact.dat" using 1:3 with lines lc rgb "red" title "Exact u", \
     "../build/output.dat" using 1:3 with points pt 7 lc rgb "red" title "Numerical Velocity"

# --- Pressure ---
set output "../validation/pressure.png"
set title "Sod Shock Tube – Pressure"
plot "../validation/sod_exact.dat" using 1:4 with lines lc rgb "green" title "Exact p", \
     "../build/output.dat" using 1:4 with points pt 7 lc rgb "green" title "Numerical Pressure"

# Reset output
set output
