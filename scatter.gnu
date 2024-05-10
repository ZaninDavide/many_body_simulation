set terminal png \
    size 4096,2048 \
    font ",36" \
    linewidth 2
set output "scatter.png"
set grid mxtics mytics
set grid xtics ytics

set multiplot layout 1,2

# L = 7.469008  # M = 4, rho = 1.2
# L = 5.928156    # M = 2, rho = 1.2
# L = 4.705180  # M = 1, rho = 1.2
L = 8.549880
# L = 36.840317

set title sprintf("Particle grid at t = 0, L = %f", L)

set xlabel "x"
set ylabel "y"
set xrange [-L:L]
set yrange [-L:L]
set xtics L/5
set ytics L/5
unset key
plot "scatter.dat" using 1:2 title "Particles" pt 7 ps 3 

set xlabel "x"
set ylabel "z"
set xrange [-L:L]
set yrange [-L:L]
set xtics L/5
set ytics L/5
unset key
plot "scatter.dat" using 1:3 title "Particles" pt 7 ps 3 