set terminal png \
    size 2048,1536 \
    font ",36" \
    linewidth 2
set output "histogram.png"
set grid mxtics mytics
set grid xtics ytics
set xrange [0:1]
set xtics 0.1

set title "Particle distances distribution g(r)"
set xlabel "r/L"
set ylabel "g(r)"

# Lineplot
# plot "histogram.dat" using 1:2 with lines

# Histogram
# set offset graph 0.05,0.05,0.05,0.0
set boxwidth 0.00666666666667 # 1 / BINS
set style fill solid 0.20 #fillstyle
set tics out nomirror
plot "histogram.dat" using ($1 + 0.00333333333334):2 with boxes title "g(r) = (dn/dV)/ρ"