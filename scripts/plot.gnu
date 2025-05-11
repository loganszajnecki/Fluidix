set xlabel "x"
set ylabel "Value"
plot "output.dat" using 1:2 with lines title "Density", \
     "output.dat" using 1:3 with lines title "Velocity", \
     "output.dat" using 1:4 with lines title "Pressure"
pause -1
