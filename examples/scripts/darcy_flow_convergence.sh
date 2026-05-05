#!/bin/bash
# darcy_flow_convergence.sh - Automates testing and plotting for Mixed Darcy solver

EXEC="./build/examples/darcy_flow"
DATA_FILE="darcy_convergence.dat"
PLOT_FILE="darcy_convergence.png"

N_VALUES="20 30 40 60 80 120 160"

# Clean up previous runs
rm -f $DATA_FILE

echo "# N    h              Err_P          Err_U          Err_V" > $DATA_FILE
echo "Running Mixed Darcy Convergence Tests..."

for N in $N_VALUES; do
    h=$(awk -v N=$N 'BEGIN {printf "%.6e", 1.0/N}')
    echo -n "  Testing N=$N (h=$h)... "

    OUTPUT=$($EXEC $N tmp_out.dat)

    # Extract numeric values from the end of the line
    ERR_P=$(echo "$OUTPUT" | grep "(P)" | awk '{print $NF}')
    ERR_U=$(echo "$OUTPUT" | grep "(U)" | awk '{print $NF}')
    ERR_V=$(echo "$OUTPUT" | grep "(V)" | awk '{print $NF}')

    if [ -z "$ERR_P" ]; then
        echo "Failed to parse output for N=$N"
        continue
    fi

    echo "Done."
    echo "$N  $h  $ERR_P  $ERR_U  $ERR_V" >> $DATA_FILE
done

rm -f tmp_out.dat

gnuplot <<- EOF
    set terminal pngcairo size 900,700 enhanced font 'Arial,12'
    set output '$PLOT_FILE'

    # Simplified Title
    set title 'Mixed Darcy Flow Convergence' font ',14'

    set xlabel 'Grid spacing h' font ',12'
    set ylabel 'Max Pointwise Error' font ',12'

    set logscale xy
    set grid mytics ytics xtics
    set format x "10^{%L}"
    set format y "10^{%L}"
    set key bottom right spacing 1.3 box opaque

    set fit quiet; set fit logfile '/dev/null'

    # Initial guesses
    pp = 2.0; cp = 1.0
    pu = 1.0; cu = 1.0

    # Fit data
    f_p(x) = cp * (x**pp); fit f_p(x) '$DATA_FILE' using 2:3 via cp, pp
    f_u(x) = cu * (x**pu); fit f_u(x) '$DATA_FILE' using 2:4 via cu, pu

    plot '$DATA_FILE' using 2:3 with points pt 7 ps 1.5 lc rgb '#1f77b4' title 'Pressure (L_{inf})', \
         f_p(x) with lines lw 2 lc rgb '#1f77b4' title sprintf('P Fit: O(h^{%.2f})', pp), \
         '$DATA_FILE' using 2:4 with points pt 5 ps 1.5 lc rgb '#d62728' title 'U-Velocity (L_{inf})', \
         f_u(x) with lines lw 2 lc rgb '#d62728' title sprintf('U Fit: O(h^{%.2f})', pu)
EOF

echo "Success! Plot saved to $PLOT_FILE."