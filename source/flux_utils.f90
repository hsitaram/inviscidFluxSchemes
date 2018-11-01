module flux_utils
      use constants
      use limiters
      implicit none
  contains
    !=========================================================================
    subroutine cons_to_primitive(u,g,primvars)

        real*8 :: u(NCVRS),primvars(NCVRS)
        real*8 :: g

        primvars(DENS_INDX) = u(RHO_INDX)
        primvars(VEL_INDX)  = u(RHOU_INDX)/u(RHO_INDX)
        primvars(PRES_INDX) = (g-1)*(u(RHOE_INDX) - 0.5*u(RHO_INDX)*primvars(VEL_INDX)**2)

    end subroutine cons_to_primitive
    !=========================================================================
    subroutine prim_to_conservative(pvars,g,u)

        real*8 :: u(NCVRS),pvars(NCVRS)
        real*8 :: g

        u(RHO_INDX) = pvars(DENS_INDX)
        u(RHOU_INDX) = pvars(DENS_INDX)*pvars(VEL_INDX)	
        u(RHOE_INDX) = pvars(PRES_INDX)/(g-1) + 0.5*pvars(DENS_INDX)*pvars(VEL_INDX)**2

    end subroutine prim_to_conservative
    !=========================================================================
    subroutine getflux(u,g,flux)

        real*8 :: g
        real*8 :: u(NCVRS)
        real*8 :: flux(NCVRS)
        real*8 :: p,vel,rho,rho_e
        real*8 :: pvars(NCVRS)

        call cons_to_primitive(u,g,pvars)

        rho  = pvars(DENS_INDX)
        vel  = pvars(VEL_INDX)
        p    = pvars(PRES_INDX)
        rho_e = u(RHOE_INDX)

        flux(RHO_INDX)  = rho*vel
        flux(RHOU_INDX) = (rho*vel**2+p)
        flux(RHOE_INDX) = (rho_e+p)*vel

    end subroutine getflux
    !=====================================================
    subroutine roeavg(leftval,rightval,rhol,rhor,avgval)

        real*8 :: leftval,rightval
        real*8 :: rhol,rhor
        real*8 :: avgval

        avgval = (leftval*sqrt(rhol)+rightval*sqrt(rhor))/(sqrt(rhol)+sqrt(rhor))

    end subroutine roeavg
    !=====================================================
    subroutine get_higherorder_states(g,u_im2,u_im1,u_i,u_ip1,uL,uR,limiter_num)

        !we are looking at face between i and i-1

        !-----------|----x-----|-----x------|-----x------|-----x------|
        !-----------|---(i-2)--|----(i-1)---|----i-------|----(i+1)---|
        !--------------------------------(my face)--------------------|

        !use reconstruction of primitive variables with limiter
        real*8, intent(in) :: u_im2(NCVRS),u_im1(NCVRS),u_i(NCVRS),u_ip1(NCVRS)
        real*8, intent(in) :: g
        integer, intent(in) :: limiter_num
        real*8, intent(out) :: ul(NCVRS),ur(NCVRS)

        real*8 :: p_im2(NCVRS),p_im1(NCVRS),p_i(NCVRS),p_ip1(NCVRS)
        real*8 :: psiL(NCVRS),psiR(NCVRS) !limiters
        real*8 :: pl(NCVRS),pr(NCVRS)
        integer :: v

        call cons_to_primitive(u_im2,g, p_im2)
        call cons_to_primitive(u_im1,g, p_im1)
        call cons_to_primitive(u_i  ,g, p_i)
        call cons_to_primitive(u_ip1,g, p_ip1)


        do v=1,NCVRS
            if(limiter_num .eq. MINMOD_LIMITER) then 
                call minmod(p_i(v)-p_im1(v), p_im1(v)-p_im2(v),psiL(v))
                call minmod(p_i(v)-p_im1(v), p_ip1(v)-p_i(v)  ,psiR(v))
            else if(limiter_num .eq. VAN_ALBADA_LIMITER) then
                call van_albada(p_i(v)-p_im1(v), p_im1(v)-p_im2(v),psiL(v))
                call van_albada(p_i(v)-p_im1(v), p_ip1(v)-p_i(v)  ,psiR(v))
            else
                print *,"limiter not implemented yet"
                stop
            endif
        enddo


        pl(:) = p_im1(:) +  0.5*psiL(:)*(p_im1-p_im2)
        pr(:) = p_i(:)   -  0.5*psiR(:)*(p_ip1-p_i)

        !ul(:) = u_im1(:) +  0.5*psiL(:)*(u_im1-u_im2)
        !ur(:) = u_i(:)   -  0.5*psiR(:)*(u_ip1-u_i)

        call prim_to_conservative(pl,g,ul)
        call prim_to_conservative(pr,g,ur)
        

    end subroutine get_higherorder_states
    !=====================================================
  end module flux_utils
