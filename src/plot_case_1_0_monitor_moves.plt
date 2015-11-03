reset
#this script plots the moves from MOTS2
#title & labels
set title "Monitor Moves" 
set xlabel "Optimisation step" 
set ylabel "Number of invocations" 
set grid
set xrange [0:771]
set yrange [0:39]
# Line styles
set border linewidth 1.5
set style line 1 linecolor rgb '#0000ff' linetype 1 linewidth 3  pointtype -1 pointsize default # blue 
set style line 2 linecolor rgb '#ff0000' linetype 2 linewidth 3  pointtype 0 pointsize default # red 
set style line 3 linecolor rgb '#00ff00' linetype 8 linewidth 3  pointtype -1 pointsize default # green 
set style line 4 linecolor rgb '#ffff00' linetype 8 linewidth 2  pointtype 2 pointsize default # yellow 
set style line 5 linecolor rgb '#ff0000' linetype 1 linewidth 2  pointtype 1 pointsize default # red 
#Legend
set key under
 plot 'monitor_data/intensify.out' with linespoints ls 5 title 'Intesification' , \
 'monitor_data/diversify.out' with linespoints ls 4 title 'Diversification', \
 'monitor_data/reduce.out' title 'Step size reduction' 
 #EPS
#set terminal postscript eps size 3.5,2.62 enhanced color \
#    font 'Helvetica,14' linewidth 2
#set output 'MOTS2_monitor_moves.eps' 
replot
set term pop
set output 
