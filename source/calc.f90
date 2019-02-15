module calc
  use vars
  use fluxschemes
  use ausm_schemes
  use cusp_schemes
  use hll_schemes
  use flux_utils
  use constants
  use inputdata
  implicit none
contains
  !======================================================================
  subroutine findtimestep(tstep)

    real*8, intent(inout) :: tstep
    integer :: i
    real*8 :: cs,delt

    cs = sqrt(gama*pres(1)/dens(1))
    tstep = dx_grid/cs 

    do i=1,npts_grid
       cs = sqrt(gama*pres(i)/dens(i))
       delt = dx_grid/cs

       if(tstep<delt) then
          tstep=delt
       endif

    enddo

  end subroutine findtimestep
  !===========================================================================
  subroutine updateresidual()

    integer :: i
    real*8  :: cvarL(NCVRS),cvarR(NCVRS)
    real*8  :: flux(NCVRS)

    real*8  :: cim2(NCVRS),cim1(NCVRS),ci(NCVRS),cip1(NCVRS)

    !interior flux
    do i=2,npts_grid

       if(secondorderflux .eq. 1) then

            cim1=consvar(i-1,:)
            ci  =consvar(i,  :)

            if(i .gt. 2) then
                cim2=consvar(i-2,:)
            else
                cim2=cim1
            endif

            if(i .lt. npts_grid) then
                cip1=consvar(i+1,:)
            else
                cip1=ci
            endif

            call get_higherorder_states(gama,cim2,cim1,ci,cip1,cvarL,cvarR,limiter_num)

        else
         cvarL(:)=consvar(i-1,:)
         cvarR(:)=consvar(i,  :)
        endif

       if(scheme_num .eq. LAXF_SCHEME) then
           call laxfrflux(cvarL,cvarR,gama,flux)

        else if(scheme_num .eq. AUSM_SCHEME) then
           call ausm(cvarL,cvarR,gama,flux)

        else if(scheme_num .eq. AUSM_PLUS_SCHEME) then
           call ausm_plus(cvarL,cvarR,gama,flux)

        else if(scheme_num .eq. AUSM_P_UP_SCHEME) then
           call ausm_plus_up(cvarL,cvarR,gama,flux)

        else if(scheme_num .eq. ECUSP_SCHEME) then
           call ecusp(cvarL,cvarR,gama,flux)

        else if(scheme_num .eq. HCUSP_SCHEME) then
           call hcusp(cvarL,cvarR,gama,flux)
        
       else if(scheme_num .eq. HLL_SCHEME) then
           call hll(cvarL,cvarR,gama,flux)

        else
           print *,"scheme not implemented yet"
           stop
        endif

        invfluxresidual(i-1,:) = invfluxresidual(i-1, :) + flux(:)
        invfluxresidual(i,  :) = invfluxresidual(i,   :) - flux(:)
    enddo


    !boundaries

    !left side
    invfluxresidual(1,  RHO_INDX) = invfluxresidual(1,  RHO_INDX) - dens(1)*vel(1)
    invfluxresidual(1, RHOU_INDX) = invfluxresidual(1, RHOU_INDX) - (dens(1)*vel(1)**2 + pres(1))
    invfluxresidual(1, RHOE_INDX) = invfluxresidual(1, RHOE_INDX) - (gama*pres(1)/(gama-1) + 0.5*dens(1)*vel(1)**2)*vel(1)

    !right side
    invfluxresidual(npts_grid,  RHO_INDX) = invfluxresidual(npts_grid,  RHO_INDX) + dens(npts_grid)*vel(npts_grid)
    
    invfluxresidual(npts_grid, RHOU_INDX) = invfluxresidual(npts_grid, RHOU_INDX) + &
        (dens(npts_grid)*vel(npts_grid)**2 + pres(npts_grid))
    
    invfluxresidual(npts_grid, RHOE_INDX) = invfluxresidual(npts_grid, RHOE_INDX) + &
        (gama*pres(npts_grid)/(gama-1) + 0.5*dens(npts_grid)*vel(npts_grid)**2)*vel(npts_grid)


  end subroutine updateresidual
  !===========================================================================
  subroutine timestepping()

      real*8 :: total_time
      real*8 :: delt
      integer :: flag
      integer :: rkstage

      flag = 0
      total_time = 0.d0
      invfluxresidual = 0.d0

      do while(total_time < final_time)

      invfluxresidual = 0.d0
      consvar_old = consvar
      call findtimestep(delt)

      !first order

      if(secondordertime .eq. 1) then

          do rkstage=1,RK23_STAGES
            invfluxresidual = 0.d0
            call updateresidual()
            consvar(:,:) = consvar_old(:,:) - RK23_COEFFS(rkstage)*(cfl*delt/dx_grid)*invfluxresidual(:,:)
          enddo

      else

          call updateresidual()
          consvar(:,:) = consvar_old(:,:) - (cfl*delt/dx_grid)*invfluxresidual(:,:)

      endif

      call updateprimvars(flag)

      if(flag==1) then
          exit
      endif


      total_time = total_time + cfl*delt

      print *,"time step and total time is ",cfl*delt,total_time

      enddo


      if(flag==1) then
          print *,"ERROR::pressure became negative"
      endif


  end subroutine timestepping
  !===========================================================================
end module calc
