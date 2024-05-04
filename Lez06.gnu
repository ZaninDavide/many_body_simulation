set terminal png \
    size 2048,1536 \
    font ",30" \
    linewidth 2
set output "Lez06.png"
set grid mxtics mytics
set grid xtics ytics

# set logscale y
# plot "measures.dat" using 1:2 title "Temperatura"

plot "measures.dat" using 1:2 title "Temperatura", \
     "measures.dat" using 1:3 title "Energia cinetica media", \
     "measures.dat" using 1:4 title "Energia potenziale media", \
     "measures.dat" using 1:5 title "Energia per particella", \
     "measures.dat" using 1:6 title "Compressibilità"