set terminal png\
    font ",36"\
    size 2048,1536
set output "plot_osc.png"
set grid mxtics mytics
set grid xtics ytics
set xrange [0:8]

set xlabel "Time"
set ylabel "Position"

plot "file_osc.dat" using 1:2 with lines