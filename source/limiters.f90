module limiters 
    use constants
    implicit none
contains
    !======================================================================
    subroutine minmod(ntr,dtr,phi)

        real*8,intent(in) :: ntr,dtr
        real*8,intent(out) :: phi

        real*8 :: r

        if(dtr .ne. 0) then
            r   = ntr/dtr
            phi = max(0.d0,min(1.d0,r))
        else
            phi = 1.d0
        endif

    end subroutine minmod
    !======================================================================
    subroutine van_albada(ntr,dtr,phi)

        real*8,intent(in) :: ntr,dtr
        real*8,intent(out) :: phi

        real*8 :: r

        if(dtr .ne. 0) then
            r   = ntr/dtr
            phi =(r**2+r)/(r**2+1)
        else
            phi = 1.d0
        endif

    end subroutine van_albada
    !======================================================================
end module limiters
