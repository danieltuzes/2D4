# arra jó, hogy a kontúr vonalakat szebbnek lehessen rajzolni

reset
set encoding utf8

ifname = "1000_64.dconf"
i_dconf_name = ifname . "_wsts_r256"
ifname_rho_t = i_dconf_name . "_rho_t.txt"
ifname_kappa = i_dconf_name . "_kappa.txt"
ifname_rho_p = i_dconf_name . "_rho_p.txt"
ifname_rho_n = i_dconf_name . "_rho_n.txt"


# egy kis statisztika a felbontás, a diszlokációszám és a helyes κ térkép cbrange megállapításához
stats ifname_rho_t matrix

res = sqrt(STATS_records) # number of points
nofdisl = STATS_sum / res / res

stats ifname_kappa matrix
cb_kappa = STATS_max > -STATS_min ? STATS_max : -STATS_min

stats ifname_rho_p matrix
rho_p_max = STATS_max
stats ifname_rho_n matrix
rho_n_max = STATS_max

cb_rho_pn = rho_p_max > rho_n_max ? rho_p_max : rho_n_max

# formázási beállítások
set samples res, res
set size sq
set xtics -0.5,0.5,0.5
set ytics -0.5,0.5,0.5
unset key
set view map
set pm3d corners2color c1

# felszín és kontúr eltárolása
# felszín
set table "cells_a.tmp"
sp [-0.5:0.5][-0.5:0.5] ifname_rho_t u (($1+0.5)/res-0.5):(($2+0.5)/res-0.5):3 matrix
set table "cells_b.tmp"
sp [-0.5:0.5][-0.5:0.5] ifname_kappa u (($1+0.5)/res-0.5):(($2+0.5)/res-0.5):3 matrix
set table "cells_c.tmp"
sp [-0.5:0.5][-0.5:0.5] ifname_rho_p u (($1+0.5)/res-0.5):(($2+0.5)/res-0.5):3 matrix
set table "cells_d.tmp"
sp [-0.5:0.5][-0.5:0.5] ifname_rho_n u (($1+0.5)/res-0.5):(($2+0.5)/res-0.5):3 matrix


# kontúr
set contour base
unset surf
load "deb_levels_t.txt"
set table "contour_a.tmp"
sp [-0.5:0.5][-0.5:0.5] ifname_rho_t u (($1+0.5)/res-0.5):(($2+0.5)/res-0.5):3 matrix
load "deb_levels_s.txt"
set table "contour_b.tmp"
sp [-0.5:0.5][-0.5:0.5] ifname_kappa u (($1+0.5)/res-0.5):(($2+0.5)/res-0.5):3 matrix
load "deb_levels_p.txt"
set table "contour_c.tmp"
sp [-0.5:0.5][-0.5:0.5] ifname_rho_p u (($1+0.5)/res-0.5):(($2+0.5)/res-0.5):3 matrix
load "deb_levels_n.txt"
set table "contour_d.tmp"
sp [-0.5:0.5][-0.5:0.5] ifname_rho_n u (($1+0.5)/res-0.5):(($2+0.5)/res-0.5):3 matrix

unset contour
unset table


set term pdfcairo size 16,7 enhanced font "DejaVu Sans Disl,18" color
# set term pngcairo size 1650,900

set o i_dconf_name . "_ts.pdf"
set multiplot layout 1,2 title sprintf("Dislocation densities\n{/=10number of dislocations: %.0f; resolution: %.0fx%.0f}",nofdisl,res,res)

set title " total (ρ_t)"
set palette defined ( 0 0 0 0, 1 1 1 1 ) negative

p [-0.5:0.5][-0.5:0.5] "cells_a.tmp" w image, "contour_a.tmp" w l lt -1 lw 1, ifname u ($1):($3>0.0?$2:1/0) with points pt "┻" tc rgb '#c1000e', "" u ($1):($3<0.0?$2:1/0) with points pt "┳" tc rgb '#005893'


set title "signed (κ)"
set palette positive
load "diverging_bent_extra_map.pal"
set cbrange [-cb_kappa:cb_kappa]
p [-0.5:0.5][-0.5:0.5] "cells_b.tmp" w image, "contour_b.tmp" w l lt -1 lw 1, ifname u ($1):($3>0.0?$2:1/0) with points pt "┻" tc rgb '#c1000e', "" u ($1):($3<0.0?$2:1/0) with points pt "┳" tc rgb '#005893'
unset multiplot

print "ts done"
set o i_dconf_name . "_pn.pdf"
set multiplot layout 1,2 title sprintf("Dislocation densities\n{/=10number of dislocations: %.0f; resolution: %.0fx%.0f}",nofdisl,res,res)
set title "positive (ρ_+)"
set cbrange [-cb_rho_pn:cb_rho_pn]
p [-0.5:0.5][-0.5:0.5] "cells_c.tmp" w image, "contour_c.tmp" w l lt -1 lw 1, ifname u ($1):($3>0.0?$2:1/0) with points pt "┻" tc rgb '#c1000e'

set title "negative (ρ_-)"
set palette negative
p [-0.5:0.5][-0.5:0.5] "cells_d.tmp" w image, "contour_d.tmp" w l lt -1 lw 1, ifname u ($1):($3<0.0?$2:1/0) with points pt "┳" tc rgb '#005893'
unset multiplot

set o
set term wxt
