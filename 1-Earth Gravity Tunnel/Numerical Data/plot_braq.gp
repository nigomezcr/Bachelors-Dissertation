set polar
set size square
unset border
unset tics
unset raxis
set key center left

set term pdf
set output "brachistochrones.pdf"
p 1 t "" lw 2 black, 'braq13.csv' u 1:2 t "Density 1" w l lt rgb "pink", 'braq13.csv' u (-1)*$1:2 t "" w l lt rgb "pink", 'braq23.csv' u 1:2 t "Density 2" w l lw 1.5 lt rgb "purple", 'braq23.csv' u ((-1)*$1):2 t "" w l lw 1.5 lt rgb "purple", 'braqprem3.csv' u 1:2 t "PREM" w l lt rgb "blue" , 'braqprem3.csv' u ((-1)*$1):2 t "" w l lt rgb "blue", 'braq16.csv' u 1:2 t "" w l lt rgb "pink", 'braq16.csv' u (-1)*$1:2 t "" w l lt rgb "pink", 'braq26.csv' u 1:2 t "" w l lw 1.5 lt rgb "purple", 'braq26.csv' u ((-1)*$1):2 t "" w l lw 1.5 lt rgb "purple", 'braqprem6.csv' u 1:2 t "" w l lt rgb "blue",'braqprem6.csv'  u (-1*$1):2 t "" w l lt rgb "blue" 

