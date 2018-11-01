module vars
  use constants
  use inputdata
  implicit none
  real*8,allocatable, save :: pres(:),dens(:),vel(:),consvar(:,:),invfluxresidual(:,:)
  real*8,allocatable, save :: consvar_old(:,:) !for second order RK
  real*8, save :: dx_grid
contains
  !===========================================================================
  subroutine initializevars()

    integer :: i

    dx_grid = 1.d0/npts_grid !domain is from 0 to 1 (hard coded)

    allocate(pres(npts_grid))
    allocate(dens(npts_grid))
    allocate(vel(npts_grid))
    allocate(consvar(npts_grid,NCVRS))
    allocate(consvar_old(npts_grid,NCVRS))
    allocate(invfluxresidual(npts_grid,NCVRS))

    do i=1,npts_grid

       if(i<=npts_grid/2) then
          pres(i)=pL
          dens(i)=rL
          vel(i)=vL
       else 
          pres(i)=pR
          dens(i)=rR
          vel(i)=vR
       endif

    enddo

   call updateconsvars()

  end subroutine initializevars
  !===========================================================================
  subroutine updateconsvars()

    integer :: i

    do i=1,npts_grid
       consvar(i,RHO_INDX)  = dens(i)
       consvar(i,RHOU_INDX) = dens(i)*vel(i)
       consvar(i,RHOE_INDX) = pres(i)/(gama-1) + 0.5*dens(i)*vel(i)*vel(i)
    enddo

    consvar_old=consvar

  end subroutine updateconsvars
  !===========================================================================
  subroutine updateprimvars(flag)

    integer :: i
    integer,intent(out) :: flag
    integer :: presnegflag=0

    do i=1,npts_grid

       dens(i) = consvar(i,RHO_INDX)
       vel(i) = consvar(i,RHOU_INDX)/dens(i)
       pres(i) = (gama-1)*(consvar(i,RHOE_INDX) - 0.5*dens(i)*vel(i)*vel(i))

       if(pres(i)<0) then
	       presnegflag=1
	       exit
       endif

    enddo

    flag = presnegflag

  end subroutine updateprimvars
  !===========================================================================
    
end module vars
