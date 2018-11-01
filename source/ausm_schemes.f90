module ausm_schemes
      use constants
      use flux_utils
      implicit none
      contains
    !=====================================================
    subroutine ausm(ul,ur,g,flux)

        real*8 :: g
        real*8 :: ul(NCVRS),ur(NCVRS)
        real*8 :: pl(NCVRS),pr(NCVRS)
        real*8 :: flux(NCVRS),fluxL(NCVRS),fluxR(NCVRS)

        real*8 :: aL,aR,vL,vR,ML,MR
        real*8 :: Mhalf,MLplus,MRminus
        real*8 :: PLplus,PRminus,mod_Mhalf

        real*8,parameter :: delta  = 0.3d0
        real*8,parameter :: beta   = 0.125d0 !1/8
        real*8,parameter :: alpha  = 0.1875d0 !3/16

        call cons_to_primitive(ul,g,pl)
        call cons_to_primitive(ur,g,pr)

        aL = sqrt(g*pl(PRES_INDX)/pl(DENS_INDX))
        aR = sqrt(g*pr(PRES_INDX)/pr(DENS_INDX))

        vL = pl(VEL_INDX)
        vR = pr(VEL_INDX)

        ML = vL/aL
        MR = vR/aR

        call split_Mach(ML, beta,  1.d0 ,MLplus);
        call split_Mach(MR, beta, -1.d0,MRminus);

        Mhalf = MLplus + MRminus

        fluxL(RHO_INDX)  = ul(RHO_INDX)*aL	
        fluxL(RHOU_INDX) = ul(RHOU_INDX)*aL	
        fluxL(RHOE_INDX) = (ul(RHOE_INDX)+pl(PRES_INDX))*aL

        fluxR(RHO_INDX)  = ur(RHO_INDX)*aR	
        fluxR(RHOU_INDX) = ur(RHOU_INDX)*aR	
        fluxR(RHOE_INDX) = (ur(RHOE_INDX)+pr(PRES_INDX))*aR

        mod_Mhalf=abs(Mhalf)

        if(mod_Mhalf .le. delta) then
            mod_Mhalf=(Mhalf**2+delta**2)/(2.d0*delta)
        endif

        flux = 0.5*Mhalf*(fluxL+fluxR)-0.5*mod_Mhalf*(fluxR-fluxL)

        !pressure term
        call split_pres(ML, alpha,  1.d0,PLplus)
        call split_pres(MR, alpha, -1.d0,PRminus)

        flux(RHOU_INDX) = flux(RHOU_INDX) +  (PLplus*pl(PRES_INDX)+PRminus*pr(PRES_INDX))


    end subroutine ausm
    !=====================================================
    subroutine ausm_plus(ul,ur,g,flux)

        real*8 :: g
        real*8 :: ul(NCVRS),ur(NCVRS)
        real*8 :: pl(NCVRS),pr(NCVRS)
        real*8 :: flux(NCVRS),fluxL(NCVRS),fluxR(NCVRS)

        real*8 :: aL,aR,vL,vR,ML,MR,HL,HR,astar2,a_half
        real*8 :: ahatL,ahatR
        real*8 :: Mhalf,MLplus,MRminus
        real*8 :: PLplus,PRminus,mod_Mhalf

        real*8,parameter :: delta  = 0.3d0
        real*8,parameter :: beta   = 0.125d0 !1/8
        real*8,parameter :: alpha  = 0.1875d0 !3/16

        call cons_to_primitive(ul,g,pl)
        call cons_to_primitive(ur,g,pr)

        aL = sqrt(g*pl(PRES_INDX)/pl(DENS_INDX))
        aR = sqrt(g*pr(PRES_INDX)/pr(DENS_INDX))

        HL=g/(g-1)*(pl(PRES_INDX)/pl(DENS_INDX)) + 0.5*pl(VEL_INDX)**2
        HR=g/(g-1)*(pr(PRES_INDX)/pr(DENS_INDX)) + 0.5*pr(VEL_INDX)**2

        vL = pl(VEL_INDX)
        vR = pr(VEL_INDX)

        !left side
        astar2 = 2.d0*(g-1.d0)/(g+1.d0)*HL
        ahatL  = astar2/max(sqrt(astar2),abs(vL))

        !right side
        astar2 = 2.d0*(g-1.d0)/(g+1.d0)*HR
        ahatR  = astar2/max(sqrt(astar2),abs(vR)) 

        a_half = min(ahatL,ahatR)

        ML = vL/a_half
        MR = vR/a_half

        call split_Mach(ML, beta,  1.d0 ,MLplus)
        call split_Mach(MR, beta, -1.d0,MRminus)

        Mhalf = MLplus + MRminus

        fluxL(RHO_INDX)  = ul(RHO_INDX)
        fluxL(RHOU_INDX) = ul(RHOU_INDX)
        fluxL(RHOE_INDX) = (ul(RHOE_INDX)+pl(PRES_INDX))

        fluxR(RHO_INDX)  = ur(RHO_INDX)
        fluxR(RHOU_INDX) = ur(RHOU_INDX)
        fluxR(RHOE_INDX) = (ur(RHOE_INDX)+pr(PRES_INDX))

        mod_Mhalf=abs(Mhalf)

        if(mod_Mhalf .le. delta) then
            mod_Mhalf=(Mhalf**2+delta**2)/(2.d0*delta)
        endif

        flux = (0.5*Mhalf*(fluxL+fluxR)-0.5*mod_Mhalf*(fluxR-fluxL))*a_half

        !pressure term
        call split_pres(ML, alpha,  1.d0,PLplus)
        call split_pres(MR, alpha, -1.d0,PRminus)

        flux(RHOU_INDX) = flux(RHOU_INDX) +  (PLplus*pl(PRES_INDX)+PRminus*pr(PRES_INDX))

    end subroutine ausm_plus
    !=====================================================
    subroutine ausm_plus_up(ul,ur,g,flux,Minf)

        real*8, intent(in) :: g
        real*8, intent(in) :: ul(NCVRS),ur(NCVRS)
        real*8, intent(in), optional   :: Minf

        real*8 :: pl(NCVRS),pr(NCVRS)
        real*8 :: flux(NCVRS),fluxL(NCVRS),fluxR(NCVRS)
        real*8 :: aL,aR,vL,vR,ML,MR,HL,HR,astar2,a_half
        real*8 :: ahatL,ahatR
        real*8 :: Mhalf,MLplus,MRminus,Minf2
        real*8 :: PLplus,PRminus,mod_Mhalf
        real*8 :: alpha,rho_half,phalf

        real*8 :: M_bar2,Mo2,Mo,fa_Mo

        real*8,parameter :: delta   = 0.3d0
        real*8,parameter :: beta    = 0.125d0 !1/8
        real*8,parameter :: sigma   = 1.d0

        !for testing if ausm+up will become ausm+
        !real*8, parameter :: Kp=0.d0
        !real*8, parameter :: Ku=0.d0

        real*8, parameter :: Kp=0.25d0 !0.75
        real*8, parameter :: Ku=0.75d0 !0.75

        if(present(Minf)) then
            Minf2=Minf**2
        else
            Minf2=0.09 !0.3**2
        endif

        call cons_to_primitive(ul,g,pl)
        call cons_to_primitive(ur,g,pr)

        aL = sqrt(g*pl(PRES_INDX)/pl(DENS_INDX))
        aR = sqrt(g*pr(PRES_INDX)/pr(DENS_INDX))

        HL=g/(g-1)*(pl(PRES_INDX)/pl(DENS_INDX)) + 0.5*pl(VEL_INDX)**2
        HR=g/(g-1)*(pr(PRES_INDX)/pr(DENS_INDX)) + 0.5*pr(VEL_INDX)**2

        vL = pl(VEL_INDX)
        vR = pr(VEL_INDX)

        rho_half=0.5d0*(pl(DENS_INDX)+pr(DENS_INDX))

        !left side
        astar2 = 2.d0*(g-1.d0)/(g+1.d0)*HL
        ahatL  = astar2/max(sqrt(astar2),abs(vL))

        !right side
        astar2 = 2.d0*(g-1.d0)/(g+1.d0)*HR
        ahatR  = astar2/max(sqrt(astar2),abs(vR)) 

        a_half = min(ahatL,ahatR)

        M_bar2 = 0.5d0*(vL**2+vR**2)/a_half**2
        Mo2 = min(1.d0,max(M_bar2,Minf2))
        Mo = sqrt(Mo2)
        fa_Mo = Mo*(2.d0-Mo)

        ML = vL/a_half
        MR = vR/a_half

        call split_Mach(ML, beta,  1.d0 ,MLplus)
        call split_Mach(MR, beta, -1.d0,MRminus)

        Mhalf = MLplus + MRminus
        Mhalf = Mhalf - (Kp/fa_Mo)*max(1-sigma*M_bar2,0.d0)* &
            (pr(PRES_INDX)-pl(PRES_INDX))/(rho_half*a_half**2)

        fluxL(RHO_INDX)  = ul(RHO_INDX)
        fluxL(RHOU_INDX) = ul(RHOU_INDX)
        fluxL(RHOE_INDX) = (ul(RHOE_INDX)+pl(PRES_INDX))

        fluxR(RHO_INDX)  = ur(RHO_INDX)
        fluxR(RHOU_INDX) = ur(RHOU_INDX)
        fluxR(RHOE_INDX) = (ur(RHOE_INDX)+pr(PRES_INDX))

        mod_Mhalf=abs(Mhalf)

        if(mod_Mhalf .le. delta) then
            mod_Mhalf=(Mhalf**2+delta**2)/(2.d0*delta)
        endif

        flux = (0.5*Mhalf*(fluxL+fluxR)-0.5*mod_Mhalf*(fluxR-fluxL))*a_half

        !pressure term
        !alpha=3.d0/16.d0
        alpha = (3.d0/16.d0)*(-4.d0 + 5.d0*fa_Mo**2)
        call split_pres(ML, alpha,  1.d0,PLplus)
        call split_pres(MR, alpha, -1.d0,PRminus)
        phalf = (PLplus*pl(PRES_INDX)+PRminus*pr(PRES_INDX))
        phalf = phalf - Ku*PLplus*PRminus*(pl(DENS_INDX)+pr(DENS_INDX))&
            *(fa_Mo*a_half)*(vR-vL)

        flux(RHOU_INDX) = flux(RHOU_INDX) + phalf 


    end subroutine ausm_plus_up 
    !=====================================================
    subroutine split_Mach(M,beta,plus_or_minus,Msplit)

        real*8,intent(in)  :: M,beta,plus_or_minus
        real*8,intent(out) :: Msplit

        !Msplit(+-) = (+-)0.25*(M (+-) 1)**2   !subsonic
        !Msplit(+-) = 0.5 * (M (+-) |M|)       !supersonic

        if(plus_or_minus .lt. 0) then

            if(abs(M) .le. 1) then
                Msplit = -0.25*(M-1)**2-beta*(M**2-1)**2
            else
                Msplit = 0.5*(M-abs(M))
            endif
        else
            if(abs(M) .le. 1) then
                Msplit = 0.25*(M+1)**2+beta*(M**2-1)**2
            else
                Msplit = 0.5*(M+abs(M))
            endif
        endif

    end subroutine split_Mach
    !=====================================================
    subroutine split_pres(M,alpha,plus_or_minus,Psplit)

        real*8, intent(in)  :: M,alpha,plus_or_minus
        real*8, intent(out) :: Psplit

        !Psplit(+-) = 0.25 (M (+-) 1)**2 (2 (-+) M)  !subsonic
        !Psplit(+-) = 0.5 * (M (+-) |M|)/M           !supersonic

        if(plus_or_minus .lt. 0) then

            if(abs(M) .le. 1) then
                Psplit = 0.25 * (M-1)**2 * (2+M) - alpha*M*(M**2-1)**2
            else
                Psplit =  0.5 * (M-abs(M))/M
            endif
        else
            if(abs(M) .le. 1) then
                Psplit = 0.25 * (M+1)**2 * (2-M) + alpha*M*(M**2-1)**2
            else
                Psplit = 0.5 * (M+abs(M))/M
            endif
        endif

    end subroutine split_pres
    !=====================================================
end module ausm_schemes
