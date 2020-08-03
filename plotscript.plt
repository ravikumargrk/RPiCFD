# set terminal pngcairo  transparent enhanced font "arial,10" fontscale 1.0 size 600, 400 
# set output 'heatmaps.2.png'
unset key
set style increment default
set view map scale 1
set style data lines
set xtics border in scale 0,0 mirror norotate  autojustify
set ytics border in scale 0,0 mirror norotate  autojustify
set ztics border in scale 0,0 nomirror norotate  autojustify
unset cbtics
set rtics axis in scale 0,0 nomirror norotate  autojustify
set title "Heat Map" 
set xrange [ 0 : 0.085 ] noreverse nowriteback
set x2range [ * : * ] noreverse writeback
set yrange [ 0 : 0.056 ] noreverse nowriteback
set y2range [ * : * ] noreverse writeback
set zrange [ * : * ] noreverse writeback
set cblabel "Temperature" 
#set cbrange [ 100 : 900] noreverse nowriteback
set rrange [ * : * ] noreverse writeback
set palette rgbformulae -21,-22,-23
## Last datafile plotted: "$map2"
plot 'result' using 2:1:3 with image
pause mouse any