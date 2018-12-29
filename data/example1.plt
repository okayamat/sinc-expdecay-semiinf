set term postscript eps enhanced "Times-Roman" 20
set output "example1.eps"
set size 0.85
set logscale y
set key spacing 3
#set key bottom left
#set xrange [0:100]
set yrange [1e-16:100]
set xlabel "{/Times-Italic=25 n}"
set ylabel "{/Times-Roman=25 Maximum Error}"
plot "Stenger1.dat" using 1:2 with lp title "Observed error (Stenger)" ls 1 lw 2 pt 2, "Stenger1.dat" using 1:3 with l title "Error bound (Theorem 2.1)" ls 3 lw 3, "Shintaku1.dat" using 1:2 with lp title "Observed error (New)" ls 1 lw 2 pt 4, "Shintaku1.dat" using 1:3 with l title "Error bound (Theorem 2.2)" ls 5 lw 3
