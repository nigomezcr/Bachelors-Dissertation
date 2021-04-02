set terminal pdf
#set output "trajectory_rho_0.pdf"
#set output "force_rho_0.pdf"
set output "force_rho_eff.pdf"


set xlabel "t"
set ylabel "r"
p 'data.csv' u 1:4 t "f(t)" w l lt 3 lw 2

set output "trajectory_rho_eff.pdf"
set xlabel "t (s)"
set ylabel "r (m)"
p 'data.csv' u 1:2 t "r(t)" w l lt 6 lw 2


















