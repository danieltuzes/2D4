reset
# ez a kód gyártja le a README.md-be a képet

set encoding utf8
set size 1,1

ifname = "1000_64.dconf"

set term pngcairo enhanced font "DejaVu Sans Disl,12" size 300,300
set out "dislocation_configuration.png"

set lmargin 3.2
set tmargin 2
set rmargin 0.6
set bmargin 1.7

set title "{/:Italic t}=0, {/:Italic N}=64, {/:Italic A}=1, {/:Italic n}=3,"
unset key
set xtics -0.5,0.5,0.5 out scale 0.5
set ytics -0.5,0.5,0.5 out scale 0.5
set size sq

p [-0.5:0.5][-0.5:0.5] ifname u ($1):($3>0.0?$2:1/0) with points pt "┻" tc rgb '#c1000e' ,\
                           '' u ($1):($3<0.0?$2:1/0) with points pt "┳" tc rgb '#005893' 

set o
set term wxt
