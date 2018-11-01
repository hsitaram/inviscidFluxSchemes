module cusp_schemes
      use constants
      use flux_utils
      implicit none
      contains
    !=====================================================
    subroutine ecusp(ul,ur,g,flux)

        !refer to page 103,chap 4 in Blazek's textbook
        !about CUSP scheme by Jameson


        !@inproceedings{tatsumi1995new,
        !  title={A new high resolution scheme for compressible viscous flows with shocks},
        !    author={Tatsumi, S and Martinelli, L and Jameson, A},
        !    booktitle={33rd Aerospace Sciences Meeting and Exhibit},
        !    pages={466},
        !    year={1995}
        !    }

        real*8, intent(in)    :: ul(NCVRS),ur(NCVRS)
        real*8, intent(in)    :: g
        real*8, intent(inout) :: flux(NCVRS)

        real*8 :: pl(NCVRS),pr(NCVRS)
        real*8 :: fl(NCVRS),fr(NCVRS)
        real*8 :: Vn,Mn,c
        real*8 :: diss(NCVRS),beta,alpha
        real*8 :: lambda_plus,lambda_minus

        call cons_to_primitive(ul,g,pl)
        call cons_to_primitive(ur,g,pr)

        call getflux(ul,g,fl) 
        call getflux(ur,g,fr)

        flux(:) = 0.5*(fl(:)+fr(:))

        Vn = 0.5*(pl(VEL_INDX)+pr(VEL_INDX))
        c  = 0.5*(sqrt(g*pl(PRES_INDX)/pl(DENS_INDX))+sqrt(g*pr(PRES_INDX)/pr(DENS_INDX)))
        Mn = Vn/c

        lambda_plus  = Vn+c
        lambda_minus = Vn-c

        alpha = abs(Mn)

        if(abs(Mn) .ge. 1.d0) then
            beta = Mn/abs(Mn) !sign of Mn
        else if( (Mn .ge. 0) .and. (Mn .lt. 1)) then
            beta = max(0.d0,(Vn+lambda_minus)/(Vn-lambda_minus))
        else
            beta = -max(0.d0,(Vn+lambda_plus)/(Vn-lambda_plus))
        endif

        diss = 0.5*alpha*c*(ur-ul) + 0.5*beta*(fr-fl)

        flux=flux-diss

    end subroutine ecusp
    !=====================================================
    subroutine hcusp(ul,ur,g,flux)

        !refer to page 103,chap 4 in Blazek's textbook
        !about CUSP scheme by Jameson


        !@inproceedings{tatsumi1995new,
        !  title={A new high resolution scheme for compressible viscous flows with shocks},
        !    author={Tatsumi, S and Martinelli, L and Jameson, A},
        !    booktitle={33rd Aerospace Sciences Meeting and Exhibit},
        !    pages={466},
        !    year={1995}
        !    }

        real*8, intent(in)    :: ul(NCVRS),ur(NCVRS)
        real*8, intent(in)    :: g
        real*8, intent(inout) :: flux(NCVRS)

        real*8 :: pl(NCVRS),pr(NCVRS)
        real*8 :: fl(NCVRS),fr(NCVRS)
        real*8 :: vn,Mn,c
        real*8 :: diss(NCVRS),beta,alpha
        real*8 :: lambda_plus,lambda_minus
        real*8 :: u_half,H_half
        real*8 :: Hl,Hr

        call cons_to_primitive(ul,g,pl)
        call cons_to_primitive(ur,g,pr)

        call getflux(ul,g,fl) 
        call getflux(ur,g,fr)

        Hl = (ul(RHOE_INDX)+pl(PRES_INDX))/pl(DENS_INDX)
        Hr = (ur(RHOE_INDX)+pr(PRES_INDX))/pr(DENS_INDX) 

        flux(:) = 0.5*(fl(:)+fr(:))

        call roeavg(pl(VEL_INDX),pr(VEL_INDX),pl(DENS_INDX),pr(DENS_INDX),u_half)
        call roeavg(Hl,Hr,pl(DENS_INDX),pr(DENS_INDX),H_half)

        Vn = u_half
        c  = sqrt( (g-1)*(H_half-0.5*u_half**2) )
        Mn = Vn/c

        lambda_plus  = (g+1)/(2.d0*g)*Vn + sqrt( ((g-1)/(2.d0*g)*Vn)**2 + c**2/g )
        lambda_minus = (g+1)/(2.d0*g)*Vn - sqrt( ((g-1)/(2.d0*g)*Vn)**2 + c**2/g )

        alpha = abs(Mn)

        if(abs(Mn) .ge. 1.d0) then
            beta = Mn/abs(Mn) !sign of Mn
        else if( (Mn .ge. 0) .and. (Mn .lt. 1)) then
            beta = max(0.d0,(Vn+lambda_minus)/(Vn-lambda_minus))
        else
            beta = -max(0.d0,(Vn+lambda_plus)/(Vn-lambda_plus))
        endif

        diss = 0.5*alpha*c*(ur-ul) + 0.5*beta*(fr-fl)

        !in H cusp the (ur-ul) term needs the 
        !enthalpy and not internal energy
        diss(RHOE_INDX)=diss(RHOE_INDX) + 0.5*alpha*c*(pr(PRES_INDX)-pl(PRES_INDX))

        flux=flux-diss

    end subroutine hcusp
    !=====================================================

end module cusp_schemes
