#!/bin/bash
# plot_fem_convergence.sh - Automates testing and plotting for the FEM Anisotropic solver

EXEC="./build/examples/fem_anisotropic"

# --- Default Parameters ---
ERROR_TYPE="max"
# We define an array of strings to hold the (NX, NY) pairs
NX_NY_PAIRS=("10 20" "20 40" "40 80" "80 160" "100 200")

# --- Command Line Parsing ---
while [[ $# -gt 0 ]]; do
    case $1 in
        --error)
            ERROR_TYPE="$2"
            shift 2
            ;;
        -h|--help)
            echo "Usage: $0 [--error max|rms]"
            exit 0
            ;;
        *)
            echo "Error: Unexpected argument '$1'"
            exit 1
            ;;
    esac
done

# --- Set Variables Based on Error Type ---
if [ "$ERROR_TYPE" = "rms" ]; then
    GREP_STR="RMS pointwise error"
    Y_LABEL="RMS Pointwise Error (L_2)"
    PREFIX="fem_rms"
else
    GREP_STR="Max pointwise error"
    Y_LABEL="Max Pointwise Error (L_{inf})"
    PREFIX="fem_max"
fi

# --- File Names ---
DATA_FILE="${PREFIX}_convergence.dat"
PLOT_FILE="${PREFIX}_convergence.png"
FIT_LOG="${PREFIX}_fit.log"

# Initialize the data file
echo "# NX      NY      Error" > "$DATA_FILE"

echo "Running convergence tests for FEM Anisotropic Diffusion ($Y_LABEL)..."

for pair in "${NX_NY_PAIRS[@]}"; do
    # Extract NX and NY from the string pair
    read -r NX NY <<< "$pair"
    
    echo -n "  Testing NX=$NX, NY=$NY ... "
    
    # Run the solver. /dev/null prevents dumping the massive solution arrays to disk.
    ERR=$($EXEC $NX $NY /dev/null | grep "$GREP_STR" | awk '{print $5}')
    
    echo "Done."
    echo "$NX  $NY  $ERR" >> "$DATA_FILE"
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
    
    set title 'Grid Convergence of FEM Solver (Rectangular Elements)' font ',14'
    set xlabel 'Grid size in Y (NY = 2*NX)' font ',12'
    set ylabel '$Y_LABEL' font ',12'
    
    # Log-Log scale
    set logscale xy
    set grid mytics ytics xtics
    set format x "10^{%L}"
    set format y "10^{%L}"
    set key top right box opaque
    
    # Redirect fit output to the local log file
    set fit quiet
    set fit logfile '$FIT_LOG'

    # Define the power-law fit function.
    # Since we plot against N instead of h, we expect a negative exponent (p ≈ -2)
    f(x) = c * (x**p)
    
    # Provide initial guesses
    c = 1.0; p = -2.0
    
    # Perform the fit using column 2 (NY) as x, and column 3 (Error) as y
    fit f(x) '$DATA_FILE' using 2:3 via c, p
    
    # Plot the data and the fitted line
    plot '$DATA_FILE' using 2:3 with points pt 7 ps 1.5 lc rgb '#2ca02c' title 'FEM Q1 Data', \
         f(x) with lines lw 2 dt 2 lc rgb '#2ca02c' title sprintf('Power Fit: O(NY^{%.2f})', p)
EOF

echo "Success! Plot saved to $PLOT_FILE."
