module inputdata
      integer, save  :: secondorderflux !0 or 1
      integer, save  :: secondordertime ! 0 or 1
      integer, save  :: scheme_num
      integer, save  :: limiter_num
      real*8,  save  :: pL,rL,vL
      real*8,  save  :: pR,rR,vR
      integer, save  :: npts_grid
      real*8,  save  :: gama
      real*8,  save  :: final_time,cfl
end module inputdata
