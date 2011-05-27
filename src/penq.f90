!--------------------------------------------------------------------
!> subroutine that calculates the coupling coefficients p^ab and 
!> q^ab 
!>
!> input    ta   : the temperature of species a
!>          tb   : the temperature of species b
!>          ma   : the mass of species a
!>          mb   : the mass of species b
!> output   pp   : array(2,2) containing the coupling coef. p
!>          qq   : array(2,2) containing the coupling coef. q
!--------------------------------------------------------------------
subroutine penq(ta,tb,ma,mb,pp,qq)

  implicit none

  real ta, tb, ma, mb, pp, qq, xab, tab, mab
  dimension pp(2,2), qq(2,2)

  tab = ta / tb
  mab = ma / mb
  xab = mab / tab

  pp(1,1) = - 3.*mab / (1+xab)**1.5
  pp(1,2) = -4.5 * xab / (1+xab)**2.5
  pp(2,1) = pp(1,2)
  pp(2,2) = - 3 * xab*(13./4.+ 2*xab + 2.5*xab**2)/(1+xab)**3.5

  qq(1,1) = - pp(1,1) / tab
  qq(1,2) = - xab*pp(1,2)
  qq(2,1) = - pp(1,2)
  qq(2,2) = 45. * xab**2 / (4.*(1+xab)**3.5)

return
end subroutine penq
