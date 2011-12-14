!> function defined by hirshman and sigmar eq. 6.74
function c3(alpha,ginv)

  implicit none

  real c3, alpha, ginv

  c3 = 1.41d0 + 3.25d0*alpha - (0.41d0 + 1.66d0*alpha)**2 &
     & /(0.59d0+alpha+1.34d0*ginv**2)

return 
end function c3
