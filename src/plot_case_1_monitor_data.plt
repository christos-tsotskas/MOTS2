reset
#this script plots the monitor data from MOTS2
#title & labels 
set title "Monitor Data" 
set xlabel "Optimisation step" 
set ylabel "Number of unsuccessful iterations" 
set grid
set xrange [0:2695]
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
#MOTS2-specific configuration settings
diversify_threshold(x)=20
intensify_threshold(x)=10
reduction_threshold(x)=35
plot 'monitor_data/case_1_i_local.out' with linespoints ls 5 title 'i-local' , \
	'monitor_data/case_1_im_size.out' with linespoints ls 4 title 'IM size' 
replot	diversify_threshold(x) with linespoints ls 1 title 'Diversification' 
replot	intensify_threshold(x) with linespoints ls 2 title 'Intensification' 
replot	reduction_threshold(x) with linespoints ls 3 title 'Step size reduction' 
#EPS 
set terminal postscript eps size 3.5,2.62 enhanced color \
    font 'Helvetica,14' linewidth 2 
set output 'MOTS2_monitor_data.eps'
replot 
set term pop
set output 
