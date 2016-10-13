set terminal postscript landscape color enhanced "Times-Roman" 16
set output 'probabilityii+1.eps'
set autoscale 
set xtic auto 
#unset yrange #[0:0.5]                        
set ytic auto
set grid xtics
set grid ytics
set xlabel "{/Symbol e}_{T_i}" font "Times-Roman, 20"
set ylabel "P_{n_i --> n_f}" font "Times-Roman, 20" 
set key right bottom	#set the legend
# unset key

set style line 1 lt 1 lc rgb "#00eeee" lw 2 pt 7 ps 1
set style line 2 lt 1 lc rgb "#008040" lw 2 pt 11 ps 1
set style line 3 lt 1 lc rgb "#ffa040" lw 2 pt 13 ps 1
set style line 4 lt 0 lc rgb "#9400d3" lw 5 pt 2 ps 1
set style line 5 lt 0 lc rgb "#191970" lw 5 pt 1 ps 1
set style line 6 lt 0 lc rgb "#006400" lw 5 pt 64 ps 1
set style line 7 lt 0 lc rgb "#f055f0" lw 5 pt 65 ps 1
set style line 8 lt 0 lc rgb "#804014" lw 5 pt 73 ps 1

plot "prob01.txt" u 1:5 w p ls 1 t "semiclassical calculation",\
"prob01.txt" u 1:7 w p ls 2 t "1st order pert. th.",\
"prob01.txt" u 1:6 w p ls 3 t "2nd order pert. th.",\
"quantum.txt" u 1:2  smooth csplines \
w l ls 4 t "quantum prob 01",\
"prob12.txt" u 1:5 notitle w p ls 1,\
"prob12.txt" u 1:7 notitle w p ls 2,\
"prob12.txt" u 1:6 notitle w p ls 3,\
"quantum.txt" u 1:3  smooth csplines \
w l ls 5 t "quantum prob 12",\
"prob23.txt" u 1:5 notitle w p ls 1,\
"prob23.txt" u 1:7 notitle w p ls 2,\
"prob23.txt" u 1:6 notitle w p ls 3,\
"quantum.txt" u 1:4  smooth csplines \
w l ls 7 t "quantum prob 23",\
"prob34.txt" u 1:5 notitle w p ls 1,\
"prob34.txt" u 1:7 notitle w p ls 2,\
"prob34.txt" u 1:6 notitle w p ls 3,\
"quantum.txt" u 1:5  smooth csplines \
w l ls 6 t "quantum prob 34",\
"prob45.txt" u 1:5 notitle w p ls 1,\
"prob45.txt" u 1:7 notitle w p ls 2,\
"prob45.txt" u 1:6 notitle w p ls 3,\
"quantum.txt" u 1:6  smooth csplines \
w l ls 8 t "quantum prob 45",\
