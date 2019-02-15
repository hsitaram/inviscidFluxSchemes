module hll_schemes
      use constants
      use flux_utils
      implicit none
      contains
    !=====================================================
    subroutine hll(ul,ur,g,flux)

        real*8 :: g
        real*8 :: ul(NCVRS),ur(NCVRS)
        real*8 :: pl(NCVRS),pr(NCVRS)
        real*8 :: flux(NCVRS),fluxL(NCVRS),fluxR(NCVRS)

        real*8 :: ahalf,vhalf,SL,SR
        real*8 :: Uhll(NCVRS),Fhll(NCVRS)

        call cons_to_primitive(ul,g,pl)
        call cons_to_primitive(ur,g,pr)

        !averages 
        ahalf = 0.5*(sqrt(g*pl(PRES_INDX)/pl(DENS_INDX)) + &
                sqrt(g*pr(PRES_INDX)/pr(DENS_INDX)) )
        vhalf = 0.5*(pl(VEL_INDX) + pr(VEL_INDX))

        SL = vhalf-ahalf
        SR = vhalf+ahalf
        
        call getflux(ul,g,fluxL)
        call getflux(ur,g,fluxR)

        if(SL .gt. 0) then
            flux=fluxL
        else
            if(SR .gt. 0) then

                Uhll=(SR*ur-SL*ul+fluxL-fluxR)/(SR-SL)
                Fhll=fluxR+SR*(Uhll-UR)
                flux=Fhll
            else
                flux=fluxR
            endif
        endif

    end subroutine hll
    !=====================================================
end module hll_schemes
