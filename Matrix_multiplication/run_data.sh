#!/bin/bash

#GENERATE DATA AND PLOT FILES


touch ./OUT_${METHOD}_${MATRIX_SIZE}/gen_data.sh
   
cat <<EOF > ./OUT_${METHOD}_${MATRIX_SIZE}/gen_data.sh
#!/bin/bash


for i in \$(seq 1 $NUM_AVG)
do
DATAFILE=OUT_${METHOD}_${MATRIX_SIZE}_\${i}.dat
touch \$DATAFILE

if [ $METHOD == "cublas" ]; then
    echo init_time comm_time comp_time cpu_gpu_time num_nodes >> \$DATAFILE
else
echo init_time comm_time comp_time num_nodes >> \$DATAFILE
fi

for ((file=1; file<=16; file*=2))
do
    tail -1 avg_\${file}_\${i}.out >> \$DATAFILE
done
done

module load python
source ~/.hpc/bin/activate
python3 average.py

rm -f *.out
rm -f OUT_${METHOD}_${MATRIX_SIZE}_*.dat

gnuplot plot_data.plt

EOF

touch ./OUT_${METHOD}_${MATRIX_SIZE}/average.py

cat <<EOF > ./OUT_${METHOD}_${MATRIX_SIZE}/average.py
import pandas 



dfs = []
for i in range(1, $NUM_AVG+1):
    dfs.append(pandas.read_csv(f"OUT_${METHOD}_${MATRIX_SIZE}_{i}.dat", sep=" "))
    

sum = 0
for i in range(len(dfs)):
    sum += dfs[i]

avg = sum / len(dfs)
avg["num_nodes"] = avg.num_nodes.astype(int)

avg.to_csv("OUT_${METHOD}_${MATRIX_SIZE}.dat", sep=" ", index=False, float_format='%f')

tot_time = avg.iloc[:, :3].sum(axis=1)
speedup = tot_time[0] / tot_time
speedup = pandas.Series(speedup, name="speedup")

n_nodes = pandas.Series([1, 2, 4, 8, 16], name="num_nodes")
speedup = pandas.concat([speedup, n_nodes], axis=1)

speedup.to_csv("speedup.dat", sep=" ", index=False, float_format='%f')
EOF


touch ./OUT_${METHOD}_${MATRIX_SIZE}/plot_data.plt

if [ $METHOD != "cublas" ]; then
cat <<EOF > ./OUT_${METHOD}_${MATRIX_SIZE}/plot_data.plt
set terminal pngcairo size 1600,768
set output 'plot.png'



set style data histograms
set style histogram rowstacked
set style fill solid border -1
set boxwidth 0.5 relative

set palette defined (0 "red", 1 "green", 2 "blue")
set key outside
set title "MAT MUL | Size=$MATRIX_SIZE | Method=$METHOD | Threads Per Task=$OMP_NUM_THREADS | Task Per Node=$NTASKS" font ",16" 
set xlabel "Number of Nodes" font ",16"
set ylabel "Average Time (s)" font ",16"

plot 'OUT_${METHOD}_${MATRIX_SIZE}.dat' using 1:xtic(4) title "Initialization", '' using 2:xtic(4) title "Commmunication", '' using 3:xtic(4) title "Computation"

set output 'plot_comp.png'
set title "MAT MUL | Size=$MATRIX_SIZE | Method=$METHOD | Threads Per Task=$OMP_NUM_THREADS | Task Per Node=$NTASKS" font ",16" 
set xlabel "Number of Nodes" font ",16"
set ylabel "Average Time (s)" font ",16"
plot 'OUT_${METHOD}_${MATRIX_SIZE}.dat' using 3:xtic(4) title "Computation"

set output 'plot_comm.png'
set title "MAT MUL | Size=$MATRIX_SIZE | Method=$METHOD | Threads Per Task=$OMP_NUM_THREADS | Task Per Node=$NTASKS" font ",16" 
set xlabel "Number of Nodes" font ",16"
set ylabel "Average Time (s)" font ",16"
plot 'OUT_${METHOD}_${MATRIX_SIZE}.dat' using 2:xtic(4) title "Communication"

set output 'plot_init.png'
set title "MAT MUL | Size=$MATRIX_SIZE | Method=$METHOD | Threads Per Task=$OMP_NUM_THREADS | Task Per Node=$NTASKS" font ",16"
set xlabel "Number of Nodes" font ",16"
set ylabel "Average Time (s)" font ",16"
plot 'OUT_${METHOD}_${MATRIX_SIZE}.dat' using 1:xtic(4) title "Initialization"  


set output 'plot_speed.png'
set title "MAT MUL SPEEDUP | Size=$MATRIX_SIZE | Method=$METHOD | Threads Per Task=$OMP_NUM_THREADS | Task Per Node=$NTASKS" font ",16" 
set xlabel "Number of Nodes" font ",16"
set ylabel "Speedup" font ",16"
set xrange [1:*]
set yrange [0:*]
set ytics 2
set xtics (1, 2, 4, 8, 16)

plot 'speedup.dat' using 2:1 with linespoints title "Speedup" lw 2, x title "Ideal Speedup" 






EOF

else
cat <<EOF > ./OUT_${METHOD}_${MATRIX_SIZE}/plot_data.plt
set terminal pngcairo size 1600,768
set output 'plot.png'



set style data histograms
set style histogram rowstacked
set style fill solid border -1
set boxwidth 0.5 relative

set palette defined (0 "red", 1 "green", 2 "blue")
set key outside
set title "MAT MUL | Size=$MATRIX_SIZE | Method=$METHOD | Threads Per Task=$OMP_NUM_THREADS | Task Per Node=$NTASKS" font ",16" 
set xlabel "Number of Nodes" font ",16"
set ylabel "Average Time (s)" font ",16"

plot 'OUT_${METHOD}_${MATRIX_SIZE}.dat' using 1:xtic(5) title "Initialization", '' using 2:xtic(5) title "Commmunication", '' using 3:xtic(5) title "Computation", '' using 4:xtic(5) title "CPU-GPU"

set output 'plot_comp.png'
set title "MAT MUL | Size=$MATRIX_SIZE | Method=$METHOD | Threads Per Task=$OMP_NUM_THREADS | Task Per Node=$NTASKS" font ",16" 
set xlabel "Number of Nodes" font ",16"
set ylabel "Average Time (s)" font ",16"
plot 'OUT_${METHOD}_${MATRIX_SIZE}.dat' using 3:xtic(5) title "Computation"

set output 'plot_comm.png'
set title "MAT MUL | Size=$MATRIX_SIZE | Method=$METHOD | Threads Per Task=$OMP_NUM_THREADS | Task Per Node=$NTASKS" font ",16" 
set xlabel "Number of Nodes" font ",16"
set ylabel "Average Time (s)" font ",16"
plot 'OUT_${METHOD}_${MATRIX_SIZE}.dat' using 2:xtic(5) title "Communication"

set output 'plot_init.png'
set title "MAT MUL | Size=$MATRIX_SIZE | Method=$METHOD | Threads Per Task=$OMP_NUM_THREADS | Task Per Node=$NTASKS" font ",16"
set xlabel "Number of Nodes" font ",16"
set ylabel "Average Time (s)" font ",16"
plot 'OUT_${METHOD}_${MATRIX_SIZE}.dat' using 1:xtic(5) title "Initialization"  

set output 'plot_cpu_gpu.png'
set title "MAT MUL | Size=$MATRIX_SIZE | Method=$METHOD | Threads Per Task=$OMP_NUM_THREADS | Task Per Node=$NTASKS" font ",16"
set xlabel "Number of Nodes" font ",16"
set ylabel "Average Time (s)" font ",16"
plot 'OUT_${METHOD}_${MATRIX_SIZE}.dat' using 4:xtic(5) title "CPU-GPU"

set output 'plot_speed.png'
set title "MAT MUL SPEEDUP | Size=$MATRIX_SIZE | Method=$METHOD | Threads Per Task=$OMP_NUM_THREADS | Task Per Node=$NTASKS" font ",16" 
set xlabel "Number of Nodes" font ",16"
set ylabel "Speedup" font ",16"
set xrange [1:*]
set yrange [0:*]
set ytics 2
set xtics (1, 2, 4, 8, 16)

plot 'speedup.dat' using 2:1 with linespoints title "Speedup" lw 2, x title "Ideal Speedup" 




EOF

fi

   


