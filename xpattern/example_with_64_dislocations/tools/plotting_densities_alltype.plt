# az összes módszerre ábrázolni a térképeket pdf-en
reset
reset session
set encoding utf8

ifn_pf = "1000_64.dconf"		# ifname prefix
ifn_rs = "_r8_"		# resolution string
print_contour = 0		# don't use if the resolution is low; high resolution comparison cannot be made with bc and gs
maptypes = "rho_t kappa rho_p rho_n"	# which maps would you like to draw
# methodes = "bc gs wsts wspn"		# which methodes of the followings: bc gs wsts wspn
methodes = "bc gs"		# which methodes of the followings: bc gs wsts wspn

ifname(maptype,methode) = sprintf("%s%s%s%s.txt",ifn_pf,methode,ifn_rs,maptype)
ilevelsfname(maptype,methode) = sprintf("%s%s%s%slevels.txt",ifn_pf,methode,ifn_rs,maptype)
odcfname(maptype) = sprintf("density_%s.tmp",maptype)
ocfname(maptype) = sprintf("contour_%s.tmp",maptype)
opfname(maptype,methode) = ifname(maptype,methode) . ".png"

normalize(x) = x > 0.5 ? normalize(x-1) : x < -0.5 ? normalize(x+1) : x

do for [methode in methodes] {
# egy kis statisztika a felbontás, a diszlokációszám és a helyes κ térkép cbrange megállapításához
print "\n\n". "\n\n". methode . "\n\n"
print "stats on " . ifname("rho_t",methode)
stats ifname("rho_t",methode) matrix

res = sqrt(STATS_records)		# fineness of the mesh
nofdisl = STATS_sum / res / res	# it will be printed out, good to check
print sprintf("res = %.0f; nofdisl = %.0f\n\n",res,nofdisl)

print "stats on " . ifname("kappa",methode)
stats ifname("kappa",methode) matrix
cb_kappa = STATS_max > -STATS_min ? STATS_max : -STATS_min

print "stats on " . ifname("rho_p",methode)
stats ifname("rho_p",methode) matrix
rho_p_max = STATS_max

print "stats on " . ifname("rho_n",methode)
stats ifname("rho_n",methode) matrix
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
do for [maptype in maptypes] {
unset contour
set surf
ofname = odcfname(maptype)
print ofname
set table ofname
sp [-0.5:0.5][-0.5:0.5] ifname(maptype,methode) u (($1+0.5)/res-0.5):(($2+0.5)/res-0.5):3 matrix
print ofname . " done"
}

# kontúr
if (print_contour && (methode eq "wspn" || methode eq "wsts")) {
set contour base
unset surf
do for [maptype in maptypes] {
if (methode eq "wspn" && (maptype eq "rho_t" || maptype eq "kappa")) {
}
else {
load ilevelsfname(maptype,methode)
ofname = ocfname(maptype)
print ofname
set table ofname
sp [-0.5:0.5][-0.5:0.5] ifname(maptype,methode) u (($1+0.5)/res-0.5):(($2+0.5)/res-0.5):3 matrix
print ofname . " done."
}
}
}
unset contour
unset table

# set term pdfcairo size 16,8 enhanced font "DejaVu Sans Disl,16"
set term pngcairo size 1650,800 enhanced font "DejaVu Sans Disl,10"

ofname = opfname("ts",methode)
print ofname
set o ofname
set multiplot layout 1,2 title sprintf("Dislocation densities, method: %s\n{/=10number of dislocations: %.0f; resolution: %.0fx%.0f}",methode,nofdisl,res,res)
print "total (ρₜ)"
set title " total (ρₜ)\nInput filename: " . ifname("rho_t",methode) noenhanced
set palette defined ( 0 0 0 0, 1 1 1 1 ) negative

if (print_contour && methode eq "wsts") {
print "contour"
p [-0.5:0.5][-0.5:0.5] odcfname("rho_t") w image, ocfname("rho_t") w l lt -1 lw 1, ifn_pf u (normalize($1)):($3>0.0?$2:1/0) with points pt "┻" tc rgb '#c1000e', "" u ($1):($3<0.0?$2:1/0) with points pt "┳" tc rgb '#005893'
}
else {
print "nocontour"
p [-0.5:0.5][-0.5:0.5] odcfname("rho_t") w image, ifn_pf u ($1):($3>0.0?$2:1/0) with points pt "┻" tc rgb '#c1000e', "" u (normalize($1)):($3<0.0?$2:1/0) with points pt "┳" tc rgb '#005893'
}
print "total (ρₜ) done"

set title "signed (κ)\nInput filename: " . ifname("kappa",methode) noenhanced
set palette positive
load "diverging_bent_extra_map.pal"
set cbrange [-cb_kappa:cb_kappa]
if (print_contour && methode eq "wsts") {
print "contour"
p [-0.5:0.5][-0.5:0.5] odcfname("kappa") w image, ocfname("kappa") w l lt -1 lw 1, ifn_pf u (normalize($1)):($3>0.0?$2:1/0) with points pt "┻" tc rgb '#c1000e', "" u (normalize($1)):($3<0.0?$2:1/0) with points pt "┳" tc rgb '#005893'
}
else {
print "nocontour"
p [-0.5:0.5][-0.5:0.5] odcfname("kappa") w image, ifn_pf u (normalize($1)):($3>0.0?$2:1/0) with points pt "┻" tc rgb '#c1000e', "" u (normalize($1)):($3<0.0?$2:1/0) with points pt "┳" tc rgb '#005893'
}
unset multiplot

print "ts done"
ofname = opfname("pn",methode)
print ofname
set o ofname
set multiplot layout 1,2 title sprintf("Dislocation densities, method: %s\n{/=10number of dislocations: %.0f; resolution: %.0fx%.0f}",methode,nofdisl,res,res)
print "positive (ρ₊)"
set title "positive (ρ₊)\nInput filename: " . ifname("rho_p",methode) noenhanced
set cbrange [-cb_rho_pn:cb_rho_pn]
if (print_contour && (methode eq "wsts" || methode eq "wspn")) {
print "contour"
p [-0.5:0.5][-0.5:0.5] odcfname("rho_p") w image, ocfname("rho_p") w l lt -1 lw 1, ifn_pf u (normalize($1)):($3>0.0?$2:1/0) with points pt "┻" tc rgb '#c1000e'
}
else {
print "nocontour"
p [-0.5:0.5][-0.5:0.5] odcfname("rho_p") w image, ifn_pf u (normalize($1)):($3>0.0?$2:1/0) with points pt "┻" tc rgb '#c1000e'
}
print "positive (ρ₊) done"

set title "negative (ρ₋)\nInput filename: " . ifname("rho_n",methode) noenhanced
set palette negative
if (print_contour && (methode eq "wsts" || methode eq "wspn")) {
print "contour"
p [-0.5:0.5][-0.5:0.5] odcfname("rho_n") w image, ocfname("rho_n") w l lt -1 lw 1, ifn_pf u (normalize($1)):($3<0.0?$2:1/0) with points pt "┳" tc rgb '#005893'
}
else {
print "nocontour"
p [-0.5:0.5][-0.5:0.5] odcfname("rho_n") w image, ifn_pf u (normalize($1)):($3<0.0?$2:1/0) with points pt "┳" tc rgb '#005893'
}
unset multiplot
}

set o
set term wxt
