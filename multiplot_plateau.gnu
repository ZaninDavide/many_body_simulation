set terminal png \
    size 2048*2,1536*3 \
    font ",45" \
    linewidth 2
set output "multiplot_plateau.png"
set grid mxtics mytics
set grid xtics ytics

# set multiplot layout 5,3
set multiplot layout 3,2

titles(i) = word("ρ=1.200 ρ=0.010 ρ=0.800 ρ=0.700 ρ=0.600 ρ=0.100 ρ=1.000 ρ=0.750 ρ=0.200 ρ=0.400 ρ=0.001 ρ=0.500 ρ=0.300 ρ=0.250",i+1)

do for [i=0:5] {
    set title titles(i)
    set xlabel "Punti per blocco (B)"
    set ylabel "Deviazione standard delle medie sui blocchi"
    set xrange [0:550]

    file = sprintf("varAvgB_%02d.dat", i)

    points_to_fit = 800
    
    f(x) = c*(1 - exp(-a*x))
    ff(x) = c - d*exp(-a*x)
    fit f(x) file using 1:2 every ::0::(points_to_fit - 1) via a, c
    fit ff(x) file using 1:2 every ::0::(points_to_fit - 1) via a, c, d
    title_f = sprintf('Fit Energia: %f - %f⋅exp(-B/%f))', c, d, 1/a)
    g(x) = cc*(1 - exp(-aa*x))
    gg(x) = cc - dd*exp(-aa*x)
    fit g(x) file using 1:3 every ::0::(points_to_fit - 1) via aa, cc
    fit gg(x) file using 1:3 every ::0::(points_to_fit - 1) via aa, cc, dd
    title_g = sprintf('Fit Compressibilità: %f - %f⋅exp(-B/%f))', cc, dd, 1/aa)

    set key horizontal top

    plot file using 1:2 title "Std Media Energia" lc "blue", \
         file using 1:3 title "Std Media Compressibilità" lc "red", \
         ff(x) with lines title sprintf("Energia: τ = %d, σ = %.4f", floor(1/a) + 1, c) lc "blue", \
         gg(x) with lines title sprintf("Compressibilità: τ = %d, σ = %.4f", floor(1/aa) + 1, cc) lc "red"
    
}

unset multiplot