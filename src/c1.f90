!> function defined by hirshman and sigmar eq. 6.72
function c1(alpha,ginv)

  implicit none

  real c1, alpha, ginv

  c1 = 1. - 0.52*alpha/(0.59+alpha+1.34*ginv**2)

return 
end function c1 
