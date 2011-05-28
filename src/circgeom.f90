!--------------------------------------------------------------------
!> this routine suplies the parameters that describe the circu-
!> lar equilibrium. in principle it is the user responcibility
!> to give the correct values of all the quantities. it is here
!> used in a simple way. the desired parameters are saved and 
!> the calls from geom than read out the quantities
!>   
!> input  :  it       : control parameter (integer) if 
!>                      1 = quantities are saved
!>                      2 = read out of saved quantities
!>                      the routine should be changed such
!>                      that all quantities are calculated 
!>                      as a function of rho
!>           rho      : flux surface label.
!> output :  rn       : the major radius of the magnetic axis
!>           e        : the inverse aspect ratio of the surface
!>           q        : the safety factor of the surface
!>           bn       : the magnetic field strength in chi = 
!>                      pi / 2. with chi being the poloidal 
!>                      angle
!--------------------------------------------------------------------
subroutine circgeom(it,rho,rn,e,q,bn)

  implicit none

  real e, rn, q, bn, rho, ek, rnk, qk, bnk
  integer it

  save ek, rnk, qk, bnk

  if (it .eq.1) then
    ek  = e
    rnk = rn
    qk  = q
    bnk = bn
  else
    e  = ek
    rn = rnk
    q  = qk
    bn = bnk
  endif

return
end subroutine circgeom 
