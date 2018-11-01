module inputoutput
	use vars
        use inputdata
  implicit none
contains
  subroutine readinpfile(filename)

    character(LEN=*) :: filename
    
    namelist /allinputs/ npts_grid,pL,rL,vL,pR,rR,vR,gama,final_time,cfl,&
    secondorderflux,secondordertime,scheme_num,limiter_num
    open(unit=2,file=filename,form='formatted',status='old')

    read (2,allinputs)

  end subroutine readinpfile
  !=====================================================================
  subroutine printoutput(filename)

    character(LEN=*) :: filename
    integer :: i
    real*8 :: x,dx

    open(unit=3,file=filename)

    dx=1.d0/npts_grid !hard coded domain from 0 to 1

    do i=1,npts_grid
       x=(i-0.5)*dx
       write(3,'(f12.4 f12.4 f12.4 f12.4 f12.4 f12.4 f12.4)') x,pres(i),dens(i),vel(i), &
       consvar(i,RHO_INDX),consvar(i,RHOU_INDX),consvar(i,RHOE_INDX)
    enddo

    close(3)

  end subroutine printoutput
  !=====================================================================

end module inputoutput
