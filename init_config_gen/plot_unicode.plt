set encoding utf8
set size 1,1

ifname = "1000_64.dconf"

set term pdfcairo enhanced font "DejaVu Sans Disl,18" color
set out sprintf("%s_unicode.pdf",ifname)

set title "{/:Italic t}=0"
unset key
set xtics -0.5,0.5,0.5 out scale 0.5
set ytics -0.5,0.5,0.5 out scale 0.5
set size sq

p [-0.5:0.5][-0.5:0.5] ifname u ($1):($3>0.0?$2:1/0) with points pt "┻" tc rgb '#c1000e' ,\
                           '' u ($1):($3<0.0?$2:1/0) with points pt "┳" tc rgb '#005893' 

set o
set term wxt
