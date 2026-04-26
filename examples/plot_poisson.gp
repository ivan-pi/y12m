# plot_poisson.gp — Visualise the output of the poisson_9pt example.
#
# Produces a two-panel PNG with a colour-map (heat-map) of the numerical
# solution and the exact solution side by side.
#
# Usage:
#   gnuplot plot_poisson.gp
#
# Requires poisson_9pt.dat (produced by the poisson_9pt program).

set terminal pngcairo size 1200,520 font "Sans,12"
set output "poisson_9pt.png"

# --- palette: blue (cold) → white → red (hot) ---
set palette defined (0 "#2166ac", 0.25 "#92c5de", 0.5 "white", \
                     0.75 "#f4a582", 1.0 "#d6604d")
set cbrange [0:1]
set cbtics 0.2

set xlabel "x" offset 0,0.5
set ylabel "y" offset 1,0
set xrange [0:1]
set yrange [0:1]
set xtics 0.25
set ytics 0.25
set size ratio 1

set pm3d map interpolate 4,4

set multiplot layout 1,2 \
    title "Steady-state heat diffusion on unit square\n9-point isotropic stencil, BC: u = 4x(1-x) on top, u = 0 elsewhere" \
    font "Sans,13"

# --- Left panel: numerical solution ---
set title "Numerical solution  (y12ma solver)"
splot "poisson_9pt.dat" using 1:2:3 notitle

# --- Right panel: exact solution ---
set title "Exact solution  u = {/Symbol S}_{n odd} [32/(n{/Symbol p})^3] sin(n{/Symbol p}x) sinh(n{/Symbol p}y)/sinh(n{/Symbol p})"
splot "poisson_9pt.dat" using 1:2:4 notitle

unset multiplot
