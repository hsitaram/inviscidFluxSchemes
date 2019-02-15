plot 'final.dat' u 1:3 w p ps 1 title "computed",'../exact_soln/exact_soln' u ($1+0.5):3 w l title "exact"
