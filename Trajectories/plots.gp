set terminal pdf
set output "trajectory_rho_eff.pdf"

set xlabel "t"
set ylabel "r"
p 'data.csv' u 1:3 t "r(t)" w l lt 1



















