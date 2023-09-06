set title "Calculating Wave load using Extrapolation Method"
set nokey
set grid
set xlabel "Time"
set ylabel "wave Load"
plot 'forces_extrapolated.dat' using 1:2 with linespoints