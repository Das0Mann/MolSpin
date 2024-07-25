# Gnuplot script to plot 4th (x) and 5th (y) columns of three .dat files

set terminal pngcairo size 800,600 enhanced font 'Verdana,10'
set output 'relaxation_plot.png'

set title "Relaxation Data"
set xlabel "Magnetic field [T]"
set ylabel "Singlet Yield"
set grid

# Define line styles
set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 ps 1.5   # Blue
set style line 2 lc rgb '#dd181f' lt 1 lw 2 pt 5 ps 1.5   # Red
set style line 3 lc rgb '#00aa00' lt 1 lw 2 pt 9 ps 1.5   # Green

plot 'staticss-timeevolution.dat' using 2:3 with linespoints linestyle 1 title 'Example Staticss', \
     'nz-relaxation-timeevolution.dat' using 2:3 with linespoints linestyle 2 title 'Nakajima-Zwanzig Relaxation', \
     'redfield-relaxation-timeevolution.dat' using 2:3 with linespoints linestyle 3 title 'Redfield Relaxation'

set output
