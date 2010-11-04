#!/bin/sh

gnuplot -persist <<EOF
plot "wave/fdtd_point.csv" using 1:2 w l
EOF
