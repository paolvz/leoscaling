reset
set autoscale xfixmin
set autoscale xfixmax
set autoscale yfixmin
set autoscale yfixmax
set size square
set pm3d map
set palette rgb 33,13,10
set terminal gif animate delay 50
set output 'test.gif'
do for [i=0:2000:40] { splot sprintf("./mat_files/matrix_%d.bin",i) bin array=62x62 format='%lf' rotate=90deg with image notitle; pause 0.1 }