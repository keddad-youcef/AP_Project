set title "Comparison between Optimization versions on different cache levels of nbody3D benchmark"

set auto x
set key left
set grid
set yrange [0:450]
set ylabel "Gflops"



set style data histogram
set style histogram cluster gap 1
set style fill solid border -1
set boxwidth 0.9
set xtic rotate by -45 scale 0
set multiplot layout 1, 1 rowsfirst


plot 'Performance/L1.dat' using 2:xtic(1) t "L1", 'Performance/L2.dat' u 2 t "L2", 'Performance/L3.dat' u 2 t "L3"


pause -1
