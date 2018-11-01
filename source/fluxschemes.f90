module fluxschemes  
    use constants
    use flux_utils
    implicit none
contains
    !======================================================================
    subroutine laxfrflux(ul,ur,g,flux)

        real*8, intent(in) :: ul(NCVRS),ur(NCVRS)
        real*8, intent(inout):: flux(NCVRS)
        real*8, intent(in) :: g

        real*8 :: presL,presR
        real*8 :: fluxL(NCVRS),fluxR(NCVRS)
        real*8 :: cs
        integer :: i
        real*8 :: vmax,vL,vR

        vL    = ul(RHOU_INDX)/ul(RHO_INDX)
        presL = (g-1)*(ul(RHOE_INDX) - 0.5*ul(RHO_INDX)*vL*vL)

        vR    = ur(RHOU_INDX)/ur(RHO_INDX)
        presR = (g-1)*(ur(RHOE_INDX) - 0.5*ur(RHO_INDX)*vR*vR)

        call getflux(ul,g,fluxL)
        call getflux(ur,g,fluxR)

        cs   = sqrt(g*max(presL/ul(RHO_INDX),presR/ur(RHO_INDX)))
        vmax = max(abs(vL),abs(vR)) 
        !print *,cs,vmax

        do i=1,3
        flux(i)=0.5*(fluxL(i)+fluxR(i)) - 0.5*(cs+vmax)*(ur(i) - ul(i))
        enddo

    end subroutine laxfrflux
    !======================================================================

end module fluxschemes
