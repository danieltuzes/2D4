set xtics ("-0.5" -0.5, "-0.25" -0.25, "0" 0, "0.25" 0.25, "0.5" 0.5)
set ytics ("0" 0, "0.8" 0.8, "1" 1, "1.2" 1.2)

set term pngcairo size 400,300 font "DejaVu Sans,10"
set o "init_pattern_gen.png"

set key bottom
set grid

set title "patterned initial condition generation with\nA=0.2 and n=3" font ",12"
set xlabel "simulation area"
set ylabel "probability density function"

rn(x) = 1-A*sin(x*n*2*pi)
rp(x) = 1+A*sin(x*n*2*pi)
A=0.2
n=3

p [-0.5:0.5][0:1.3] rp(x) t "PDF(ρ_+) = 1 + A * sin(x * n * 2 π)" lc rgb "#c1000e" lw 2, rn(x) t "PDF(ρ_- ) = 1  -  A * sin(x * n * 2 π)" lc rgb "#005893" lw 2 dt "--"
set o
set term wxt
