set title "Calculating Wave load using Wheeler Stretching Method"
set nokey
set grid
set xlabel "Time"
set ylabel "wave Load"
plot 'forces_wheeler.dat' using 1:2 with linespoints