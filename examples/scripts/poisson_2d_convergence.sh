#!/bin/bash
# poisson_2d_convergence.sh - Automates testing and plotting for the 9-pt Poisson solver

EXEC="./build/examples/poisson_9pt"

# --- Default Parameters ---
ERROR_TYPE="max"
N_VALUES=""

# --- Command Line Parsing ---
while [[ $# -gt 0 ]]; do
    case $1 in
        --error)
            ERROR_TYPE="$2"
            shift 2
            ;;
        -h|--help)
            echo "Usage: $0 [--error max|rms] \"N1 N2 N3 ...\""
            echo "Example: $0 --error rms \"\$(seq 10 10 100)\""
            exit 0
            ;;
        *)
            if [ -z "$N_VALUES" ]; then
                N_VALUES="$1"
                shift
            else
                echo "Error: Unexpected argument '$1'"
                exit 1
            fi
            ;;
    esac
done

# If no N values provided, default to a sensible logarithmic-ish list
if [ -z "$N_VALUES" ]; then
    N_VALUES="5 10 20 40 80 160 200"
fi

# --- Set Variables Based on Error Type ---
if [ "$ERROR_TYPE" = "rms" ]; then
    GREP_STR="RMS pointwise error"
    Y_LABEL="RMS Pointwise Error (L_2)"
    PREFIX="poisson_rms"
else
    GREP_STR="Max pointwise error"
    Y_LABEL="Max Pointwise Error (L_{inf})"
    PREFIX="poisson_max"
fi

# --- File Names ---
DATA_FILE="${PREFIX}_convergence.dat"
PLOT_FILE="${PREFIX}_convergence.png"
FIT_LOG="${PREFIX}_fit.log"

# Initialize the data file
echo "# N       h             Error" > $DATA_FILE

echo "Running convergence tests for Poisson 9-point ($Y_LABEL)..."
echo "Grid sizes (N): $N_VALUES"

for N in $N_VALUES; do
    # Calculate h = 1 / (N - 1)
    h=$(awk -v N=$N 'BEGIN {printf "%.6e", 1.0/(N-1)}')
    echo -n "  Testing N=$N (h=$h)... "
    
    # Run the solver. 
    # $N is the first positional argument, /dev/null is the second (output_file).
    ERR=$($EXEC $N /dev/null | grep "$GREP_STR" | awk '{print $5}')
    
    echo "Done."
    echo "$N  $h  $ERR" >> $DATA_FILE
done

echo "Testing complete. Generating plot and fitting convergence rate..."
echo "Fit details will be saved to: $FIT_LOG"

# Clear out the old log file if it exists
rm -f "$FIT_LOG"

# Run gnuplot with a heredoc
gnuplot <<- EOF
    # Set up the terminal and output
    set terminal pngcairo size 800,600 enhanced font 'Arial,12'
    set output '$PLOT_FILE'
    
    set title 'Grid Convergence of Poisson Solver (9-point)' font ',14'
    set xlabel 'Grid spacing h = 1/(N-1)' font ',12'
    set ylabel '$Y_LABEL' font ',12'
    
    # Log-Log scale
    set logscale xy
    set grid mytics ytics xtics
    set format x "10^{%L}"
    set format y "10^{%L}"
    set key top left box opaque
    
    # Redirect fit output to the local log file
    set fit quiet
    set fit logfile '$FIT_LOG'

    # Define the power-law fit function
    f(x) = c * (x**p)
    
    # Provide initial guesses (p=2 is a safe starting point)
    c = 1.0; p = 2.0
    
    # Perform the fit using column 2 (h) and column 3 (Error)
    fit f(x) '$DATA_FILE' using 2:3 via c, p
    
    # Plot the data and the fitted line
    plot '$DATA_FILE' using 2:3 with points pt 7 ps 1.5 lc rgb '#d62728' title '9-point Data', \
         f(x) with lines lw 2 dt 2 lc rgb '#d62728' title sprintf('Power Fit: O(h^{%.2f})', p)
EOF

echo "Success! Plot saved to $PLOT_FILE."