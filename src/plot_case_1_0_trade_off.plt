 reset
#this script plots the monitor data from MOTS2
#title & labels
set title "Trade-off"
set xlabel "Pressure Drop"
set ylabel "Vorticity"
set grid
#set xrange [20:500]
#set yrange [0:14 ]
# Line styles
set border linewidth 1.5
set style line 1 linecolor rgb '#0000ff' linetype 1 linewidth 3  pointtype -1 pointsize default # blue
set style line 2 linecolor rgb '#ff0000' linetype 2 linewidth 3  pointtype 0 pointsize default # red
set style line 3 linecolor rgb '#00ff00' linetype 8 linewidth 3  pointtype -1 pointsize default # green
set style line 4 linecolor rgb '#ffff00' linetype 8 linewidth 2  pointtype 2 pointsize default # yellow
set style line 5 linecolor rgb '#ff0000' linetype 1 linewidth 2  pointtype 1 pointsize default # red
#Legend
set key under
#MOTS2-specific configuration settings
plot 'memories/HISTORY.txt' u 2:1 title 'History', 'memories/MTM.txt' u 2:1 title 'Pareto Front'
#selecting designs:
replot "<echo '-1 1'" title 'Datum design' pointtype 6 pointsize 3 linecolor rgb '#ff0000'
#replot "<echo '-1.7103 1.09332'" title 'Compromise design' pointtype 4 pointsize 3 linecolor rgb '#00ff00'
#replot "<echo '-3.88 2.09'" title 'Maximum Vorticity' pointtype 14 pointsize 3 linecolor rgb '#0000ff'		#extreme 1
#replot "<echo '-0.0687 0.083426'" title 'Minimum Pressure Drop' pointtype 12 pointsize 3  linecolor rgb '#ff00ff' #extreme 2
#EPS
set terminal postscript eps size 3.5,2.62 enhanced color \
   font 'Helvetica,14' linewidth 2
set output 'MOTS2_trade_off.eps'
replot
set term pop
set output 
