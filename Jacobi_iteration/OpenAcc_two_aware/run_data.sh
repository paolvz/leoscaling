#!/bin/bash

#GENERATE DATA AND PLOT


touch ./OUT_${MATRIX_SIZE}_${ITERATIONS}/gen_data.sh
   
cat <<EOF > ./OUT_${MATRIX_SIZE}_$ITERATIONS/gen_data.sh
#!/bin/bash

set -e

for i in 1 2 3 4 5
do
DATAFILE=OUT_${MATRIX_SIZE}_${ITERATIONS}_\${i}.dat
touch \$DATAFILE
echo init_time comm_time comp_time num_nodes >> \$DATAFILE

for ((file=1; file<=16; file*=2))
do
    tail -1 avg_\${file}_\${i}.out >> \$DATAFILE
done
done

module load python
source ~/.hpc/bin/activate
python3 average.py



gnuplot plot_data.plt

rm -f avg_*.out
rm -f OUT_${MATRIX_SIZE}_${ITERATIONS}_*.dat
EOF

cat <<EOF > ./OUT_${MATRIX_SIZE}_$ITERATIONS/average.py
import pandas 



dfs = []
for i in range(1, 6):
    dfs.append(pandas.read_csv(f"OUT_${MATRIX_SIZE}_${ITERATIONS}_{i}.dat", sep=" "))
    

sum = 0
for i in range(len(dfs)):
    sum += dfs[i]

avg = sum / len(dfs)
avg["num_nodes"] = avg.num_nodes.astype(int)

avg.to_csv("OUT_${MATRIX_SIZE}_${ITERATIONS}.dat", sep=" ", index=False, float_format='%f')

tot_time = avg.iloc[:, :3].sum(axis=1)
speedup = tot_time[0] / tot_time
speedup = pandas.Series(speedup, name="speedup")

n_nodes = pandas.Series([1, 2, 4, 8, 16], name="num_nodes")
speedup = pandas.concat([speedup, n_nodes], axis=1)

speedup.to_csv("speedup.dat", sep=" ", index=False, float_format='%f')
EOF


touch ./OUT_${MATRIX_SIZE}_${ITERATIONS}/plot_data.plt

cat <<EOF > ./OUT_${MATRIX_SIZE}_${ITERATIONS}/plot_data.plt
set terminal pngcairo size 1600,768
set output 'plot.png'



set style data histograms
set style histogram rowstacked
set style fill solid border -1
set boxwidth 0.5 relative

set palette defined (0 "red", 1 "green", 2 "blue")
set key outside
set title "GPU JACOBI OpenACC | Size=$MATRIX_SIZE | Iterations=$ITERATIONS | Threads Per Task=$OMP_NUM_THREADS | Task Per Node=$NTASKS" font ",16" 
set xlabel "Number of Nodes" font ",16"
set ylabel "Average Time (s)" font ",16"

plot 'OUT_${MATRIX_SIZE}_${ITERATIONS}.dat' using 1:xtic(4) title "Initialization", '' using 2:xtic(4) title "Commmunication", '' using 3:xtic(4) title "Computation"



set output 'plot_comp.png'
set title "GPU JACOBI OpenACC | Size=$MATRIX_SIZE | Iterations=$ITERATIONS | Threads Per Task=$OMP_NUM_THREADS | Task Per Node=$NTASKS" font ",16" 
set xlabel "Number of Nodes" font ",16"
set ylabel "Average Time (s)" font ",16"
plot 'OUT_${MATRIX_SIZE}_${ITERATIONS}.dat' using 3:xtic(4) title "Computation"

set output 'plot_comm.png'
set title "GPU JACOBI OpenACC | Size=$MATRIX_SIZE | Iterations=$ITERATIONS | Threads Per Task=$OMP_NUM_THREADS | Task Per Node=$NTASKS" font ",16" 
set xlabel "Number of Nodes" font ",16"
set ylabel "Average Time (s)" font ",16"
plot 'OUT_${MATRIX_SIZE}_${ITERATIONS}.dat' using 2:xtic(4) title "Communication"

set output 'plot_init.png'
set title "GPU JACOBI OpenACC | Size=$MATRIX_SIZE | Iterations=$ITERATIONS | Threads Per Task=$OMP_NUM_THREADS | Task Per Node=$NTASKS" font ",16" 
set xlabel "Number of Nodes" font ",16"
set ylabel "Average Time (s)" font ",16"
plot 'OUT_${MATRIX_SIZE}_${ITERATIONS}.dat' using 1:xtic(4) title "Initialization"

set output 'speedup.png'
set title "GPU JACOBI OpenACC | Size=$MATRIX_SIZE | Iterations=$ITERATIONS | Threads Per Task=$OMP_NUM_THREADS | Task Per Node=$NTASKS" font ",16"
set xlabel "Number of Nodes" font ",16"
set ylabel "Speedup" font ",16"
set xrange [1:*]
set yrange [0:*]
set ytics 2
set xtics (1, 2, 4, 8, 16)

plot 'speedup.dat' using 2:1 with linespoints title "Speedup" lw 2, x title "Ideal Speedup" 


EOF


   


