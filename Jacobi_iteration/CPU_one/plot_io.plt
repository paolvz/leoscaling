reset
set autoscale xfixmin
set autoscale xfixmax
set autoscale yfixmin
set autoscale yfixmax
set size square
set pm3d map
set palette rgb 33,13,10
splot "matrix.bin" bin array=62x62 format='%lf' rotate=90deg notitle with image
