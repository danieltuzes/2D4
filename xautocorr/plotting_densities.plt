reset
res = 256
i_dconf_name = "1000_64.dconf"
ifname_rho_t = i_dconf_name . "_rho_t.txt"
ifname_kappa = i_dconf_name . "_kappa.txt"

set view map

xmax = res-0.5
ymax = res-0.5
sp [-0.5:xmax][-0.5:ymax] ifname_kappa matrix w image
cbrange = GPVAL_CB_MAX > -GPVAL_CB_MIN ? GPVAL_CB_MAX : -GPVAL_CB_MIN

set term pngcairo size 1650,900
set o i_dconf_name . ".png"



set multiplot layout 1,2 title "Dislocation densities"
set size sq

set title " total (rho_t)"
set palette defined ( 0 0 0 0, 1 1 1 1 ) negative
sp [-0.5:xmax][-0.5:ymax] ifname_rho_t matrix w image


set title "signed (kappa)"
set palette positive
load "diverging_bent_extra_map.pal"
set cbrange [-cbrange:cbrange]
sp [-0.5:xmax][-0.5:ymax] ifname_kappa matrix w image
unset multiplot
set o
set term wxt
