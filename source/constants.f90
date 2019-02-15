module constants
    implicit none

    integer, parameter :: NCVRS = 3
    integer, parameter :: RHO_INDX  = 1
    integer, parameter :: RHOU_INDX = 2 
    integer, parameter :: RHOE_INDX = 3

    integer, parameter :: DENS_INDX = 1
    integer, parameter :: VEL_INDX  = 2
    integer, parameter :: PRES_INDX = 3

    integer, parameter :: RK23_STAGES = 2
    real*8, parameter  :: RK23_COEFFS(2) = (/0.5,1.0/)

    integer, parameter :: LAXF_SCHEME = 1
    integer, parameter :: AUSM_SCHEME = 2
    integer, parameter :: AUSM_PLUS_SCHEME = 3
    integer, parameter :: AUSM_P_UP_SCHEME = 4
    integer, parameter :: ECUSP_SCHEME = 5
    integer, parameter :: HCUSP_SCHEME = 6
    integer, parameter :: HLL_SCHEME = 7

    integer, parameter :: MINMOD_LIMITER = 1
    integer, parameter :: VAN_ALBADA_LIMITER = 2
end module constants
