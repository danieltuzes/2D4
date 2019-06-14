reset
res = 256
i_dconf_name = "1000_64.dconf"
ifname_rho_p = i_dconf_name . "_rho_p.txt"
ifname_rho_n = i_dconf_name . "_rho_n.txt"

set view map

xmax = res-0.5
ymax = res-0.5
sp [-0.5:xmax][-0.5:ymax] ifname_rho_p matrix w image
cbrange_rho_p = GPVAL_CB_MAX > -GPVAL_CB_MIN ? GPVAL_CB_MAX : -GPVAL_CB_MIN

sp [-0.5:xmax][-0.5:ymax] ifname_rho_n matrix w image
cbrange_rho_n = GPVAL_CB_MAX > -GPVAL_CB_MIN ? GPVAL_CB_MAX : -GPVAL_CB_MIN

set term pngcairo size 1650,900
set o i_dconf_name . "_pn.png"



set multiplot layout 1,2 title "Dislocation densities"
set size sq

set title "positive (rho_p)"
load "diverging_bent_extra_map.pal"
set cbrange [-cbrange_rho_p:cbrange_rho_p]
sp [-0.5:xmax][-0.5:ymax] ifname_rho_p matrix w image


set title "negative (rho_n)"
load "diverging_bent_extra_map.pal"
set palette negative
set cbrange [-cbrange_rho_n:cbrange_rho_n]
sp [-0.5:xmax][-0.5:ymax] ifname_rho_n matrix w image
unset multiplot
set o
set term wxt
