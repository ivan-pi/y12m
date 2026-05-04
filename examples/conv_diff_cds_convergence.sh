#!/bin/bash
# plot_convdiff_convergence.sh - Automates testing and plotting for the Convection-Diffusion solver

EXEC="./build/examples/conv_diff_cds"

# --- Default Parameters ---
ERROR_TYPE="max"
N_VALUES=""
ALPHA="0.1"
TEST_CASE="2" # Default to the mixed Dirichlet/Neumann case

# --- Command Line Parsing ---
while [[ $# -gt 0 ]]; do
    case $1 in
        --error)
            ERROR_TYPE="$2"
            shift 2
            ;;
        --alpha)
            ALPHA="$2"
            shift 2
            ;;
        --case)
            TEST_CASE="$2"
            shift 2
            ;;
        -h|--help)
            echo "Usage: $0 [--error max|rms] [--alpha value] [--case 1|2] \"N1 N2 N3 ...\""
            echo "Example: $0 --error rms --alpha 0.05 --case 2 \"\$(seq 10 10 100)\""
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
    N_VALUES="10 20 40 60 80 100"
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

# --- File Names (Now include the test case number) ---
DATA_FILE="convdiff_case${TEST_CASE}_${PREFIX}_convergence.dat"
PLOT_FILE="convdiff_case${TEST_CASE}_${PREFIX}_convergence.png"
FIT_LOG="convdiff_case${TEST_CASE}_${PREFIX}_fit.log"

# Initialize the data file
echo "# N        h              ${PREFIX}_Error_CDS" > $DATA_FILE

echo "Running convergence tests for Case $TEST_CASE, $Y_LABEL (Alpha = $ALPHA)..."
echo "Grid sizes (N): $N_VALUES"

for N in $N_VALUES; do
    h=$(awk -v N=$N 'BEGIN {printf "%.6e", 1.0/(N-1)}')
    echo -n "  Testing N=$N (h=$h)... "

    # Run CDS driver and extract the specified error
    # Note: We pass $N, $ALPHA, and /dev/null (to suppress the .dat file generation)
    ERR_CDS=$($EXEC --case $TEST_CASE $N $ALPHA /dev/null | grep "$GREP_STR" | awk '{print $5}')

    echo "Done."
    echo "$N  $h  $ERR_CDS" >> $DATA_FILE
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

    set title 'Grid Convergence of Convection-Diffusion Solver (Case $TEST_CASE, $PREFIX error, {/Symbol a}=$ALPHA)' font ',14'
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

    # Define the power-law fit function
    f_cds(x) = c_cds * (x**p_cds)

    # Provide initial guesses
    c_cds = 1.0; p_cds = 2.0

    # Perform the fit
    fit f_cds(x) '$DATA_FILE' using 2:3 via c_cds, p_cds

    # Plot
    plot '$DATA_FILE' using 2:3 with points pt 7 ps 1.5 lc rgb '#1f77b4' title 'Central Diff (CDS)', \
         f_cds(x) with lines lw 2 dt 2 lc rgb '#1f77b4' title sprintf('Fit: O(h^{%.2f})', p_cds)
EOF

echo "Success! Plot saved to $PLOT_FILE."