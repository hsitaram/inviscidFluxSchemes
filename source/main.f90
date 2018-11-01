program shock_tube
	use inputoutput
        use vars
	use calc	
	implicit none


	call readinpfile("stube.inp")
	call initializevars()	
	call printoutput("initial.dat")
	call timestepping()
	call printoutput("final.dat")

end program shock_tube
