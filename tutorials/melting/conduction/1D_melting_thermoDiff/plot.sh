#!/bin/bash

foamLog log >/dev/null

gnuplot -persist > /dev/null 2>&1 <<EOF
	set title "residual vs. Iteration"
	set xlabel "Iteration"
	set ylabel "Residual"
	plot "logs/TFinalRes_1" with lines
EOF
