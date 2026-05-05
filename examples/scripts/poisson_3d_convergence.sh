#!/bin/bash
# poisson3d_convergence.sh - Automates testing and plotting for the 3D Poisson solver

EXEC="./build/examples/poisson_3d"

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
            echo "Example: $0 --error rms \"5 10 15 20 25\""
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

# If no N values provided, default to a sensible list for 3D
# (Keeping values small due to rapid memory fill-in and solve times)
if [ -z "$N_VALUES" ]; then
    N_VALUES="8 12 16 20 24"
fi

# --- Set Variables Based on Error Type ---
if [ "$ERROR_TYPE" = "rms" ]; then
    GREP_STR="RMS pointwise error"
    Y_LABEL="RMS Pointwise Error (L_2)"
    PREFIX="poisson3d_rms"
else
    GREP_STR="Max pointwise error"
    Y_LABEL="Max Pointwise Error (L_{inf})"
    PREFIX="poisson3d_max"
fi

# --- File Names ---
DATA_FILE="${PREFIX}_convergence.dat"
PLOT_FILE="${PREFIX}_convergence.png"
FIT_LOG="${PREFIX}_fit.log"

# Initialize the data file
echo "# N       h             Error" > $DATA_FILE

echo "Running convergence tests for 3D Poisson solver ($Y_LABEL)..."
echo "Grid sizes (N): $N_VALUES"

for N in $N_VALUES; do
    # Calculate h = 1 / (N - 1)
    h=$(awk -v N=$N 'BEGIN {printf "%.6e", 1.0/(N-1)}')
    echo -n "  Testing N=$N (h=$h)... "
    
    # Run the solver. 
    # kx, ky, kz will default to 1 1 1 as per the solver's help menu.
    ERR=$($EXEC $N | grep "$GREP_STR" | awk '{print $5}')
    
    # Check if ERR is empty (e.g., if the solver crashed due to OOM)
    if [ -z "$ERR" ]; then
        echo "Failed! (Solver likely ran out of memory or crashed)"
        continue
    fi

    echo "Done (Error: $ERR)."
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
    
    set title 'Grid Convergence of 3D Poisson Solver' font ',14'
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
    
    # Provide initial guesses (p=2 is a safe starting point for Poisson)
    c = 1.0; p = 2.0
    
    # Perform the fit using column 2 (h) and column 3 (Error)
    fit f(x) '$DATA_FILE' using 2:3 via c, p
    
    # Plot the data and the fitted line
    plot '$DATA_FILE' using 2:3 with points pt 7 ps 1.5 lc rgb '#1f77b4' title '3D Solver Data', \
         f(x) with lines lw 2 dt 2 lc rgb '#1f77b4' title sprintf('Power Fit: O(h^{%.2f})', p)
EOF

echo "Success! Plot saved to $PLOT_FILE."
