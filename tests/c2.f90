!> function defined by hirshman and sigmar eq. 6.73
function c2(alpha,ginv)

  implicit none

  real c2, alpha, ginv

  c2 = 1.5 - (0.29+1.20*alpha)/(0.59+alpha+1.34*ginv**2)

return 
end function c2
