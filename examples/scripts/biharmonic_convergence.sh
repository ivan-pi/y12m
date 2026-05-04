#!/bin/bash
# plot_convergence.sh - Automates testing and plotting for the Biharmonic solver

EXEC="./build/examples/biharmonic_13pt"

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

# If no N values provided, default to a sensible list
if [ -z "$N_VALUES" ]; then
    N_VALUES="10 20 30 40 60 80 100"
fi

# --- Set Variables Based on Error Type ---
if [ "$ERROR_TYPE" = "rms" ]; then
    GREP_STR="RMS pointwise error"
    Y_LABEL="RMS Pointwise Error (L_2)"
    PREFIX="rms"
else
    GREP_STR="Max pointwise error"
    Y_LABEL="Max Pointwise Error (L_{inf})"
    PREFIX="max"
fi

# --- File Names ---
DATA_FILE="biharmonic_${PREFIX}_convergence.dat"
PLOT_FILE="biharmonic_${PREFIX}_convergence.png"
FIT_LOG="biharmonic_${PREFIX}_fit.log"

# Initialize the data file
echo "# N       h             ${PREFIX}_Error_13pt   ${PREFIX}_Error_21pt" > $DATA_FILE

echo "Running convergence tests for $Y_LABEL..."
echo "Grid sizes (N): $N_VALUES"

for N in $N_VALUES; do
    h=$(awk -v N=$N 'BEGIN {printf "%.6e", 1.0/(N-1)}')
    echo -n "  Testing N=$N (h=$h)... "

    # Run 13-point and extract the specified error
    ERR_13=$($EXEC --stencil 13 $N /dev/null | grep "$GREP_STR" | awk '{print $5}')

    # Run 21-point and extract the specified error
    ERR_21=$($EXEC --stencil 21 $N /dev/null | grep "$GREP_STR" | awk '{print $5}')

    echo "Done."
    echo "$N  $h  $ERR_13  $ERR_21" >> $DATA_FILE
done

echo "Testing complete. Generating plot and fitting convergence rates..."
echo "Fit details will be saved to: $FIT_LOG"

# Clear out the old log file if it exists so we only see the latest fit
rm -f "$FIT_LOG"

# Run gnuplot with a heredoc
gnuplot <<- EOF
    # Set up the terminal and output
    set terminal pngcairo size 800,600 enhanced font 'Arial,12'
    set output '$PLOT_FILE'

    set title 'Grid Convergence of Biharmonic Solver ($PREFIX error)' font ',14'
    set xlabel 'Grid spacing h = 1/(N-1)' font ',12'
    set ylabel '$Y_LABEL' font ',12'

    # Log-Log scale
    set logscale xy
    set grid mytics ytics xtics
    set format x "10^{%L}"
    set format y "10^{%L}"
    set key top left box opaque

    # Redirect fit output to the local log file instead of terminal
    set fit quiet
    set fit logfile '$FIT_LOG'

    # Define the power-law fit functions
    f13(x) = c13 * (x**p13)
    f21(x) = c21 * (x**p21)

    # Provide initial guesses
    c13 = 1.0; p13 = 2.0
    c21 = 1.0; p21 = 2.0

    # Perform the fits
    fit f13(x) '$DATA_FILE' using 2:3 via c13, p13
    fit f21(x) '$DATA_FILE' using 2:4 via c21, p21

    # Plot
    plot '$DATA_FILE' using 2:3 with points pt 7 ps 1.5 lc rgb '#1f77b4' title '13-point Data', \
         f13(x) with lines lw 2 dt 2 lc rgb '#1f77b4' title sprintf('13-pt Fit: O(h^{%.2f})', p13), \
         '$DATA_FILE' using 2:4 with points pt 5 ps 1.5 lc rgb '#ff7f0e' title '21-point Data', \
         f21(x) with lines lw 2 dt 2 lc rgb '#ff7f0e' title sprintf('21-pt Fit: O(h^{%.2f})', p21)
EOF

echo "Success! Plot saved to $PLOT_FILE."