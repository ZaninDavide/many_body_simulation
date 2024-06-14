set terminal png \
    size 2048,1536 \
    font ",30" \
    linewidth 2
set output "varAvgB.png"
set grid mxtics mytics
set grid xtics ytics

set xlabel "Punti per blocco (B)"
set ylabel "Deviazione standard delle medie sui blocchi"

plot "varAvgB.dat" using 1:2 title "Std Avg Energy", \
     "varAvgB.dat" using 1:3 title "Std Avg Compressibility"