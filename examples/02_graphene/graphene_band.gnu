set style data dots
set nokey
set xrange [0: 4.01697]
set yrange [-10.42650 :  5.55540]
set arrow from  1.47031, -10.42650 to  1.47031,   5.55540 nohead
set arrow from  2.31920, -10.42650 to  2.31920,   5.55540 nohead
set xtics ("G"  0.00000,"M"  1.47031,"K"  2.31920,"G"  4.01697)
 plot "graphene_band.dat"
