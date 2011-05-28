!     SUBROUTINE THAT CALCULATES THE TRANSPORT COEFFICIENTS 
!     FOR THE PARTICLE FLUX OF A TRACE IMPURITY IN AN HYDROGEN
!     PLASMA WITH AN IMPURITY. (PFIRSCH SCHLUETER REGIME)
subroutine psflux(t,nsm,ncm,mas,den,zsp,tau,l1,l2,l3,l11tt,l11th,l12th)
 
  implicit none

  real d0,d1,d2,d3,t,mm,nn,mmh,nnh,mas,mt,mi,zhi
  real l1,l2,l3,l11tt,l11th,l12th,c1,c2,den,tau,zsp
  integer i,j, nsm, ncm

  dimension  mm(3,3), nn(3,3), mmh(3,3), nnh(3,3), mas(nsm)
  dimension  tau(nsm,nsm),zsp(nsm,ncm),den(nsm,ncm)

  mt = mas(3)
  mi = mas(2)
  call MENN(T,T,MT,MI,MM,NN)

  do i = 1, 3 
    do j = 1, 3
      mmh(i,j) = mm(i,j)-mm(i,3)*mm(3,j)/mm(3,3)
      nnh(i,j) = nn(i,j)-mm(i,3)*nn(3,j)/mm(3,3)
    end do 
  end do 
  
  d0 = - sqrt(mt/mi)*(mmh(1,1)- mmh(1,2)*mmh(2,1)/mmh(2,2))
  d1 = 0.88*sqrt(mt/mi)*(nnh(1,2)-mmh(1,2)*nnh(2,2)/mmh(2,2) &
     &     -4./15. * (nnh(1,3)-mmh(1,2)*nnh(2,3)/mmh(2,2)) )
  d2 = mmh(2,1)/mmh(2,2)
  d3 = ((zsp(2,1)/zsp(3,1))**2*d2 - d1)/(0.293 + 0.592)
  zhi = den(2,1)*zsp(2,1)**2/den(1,1)

  l1 = - tau(2,3)*(d0+ 0.236*d3)- tau(1,3)*(1-c1(zhi,0))
  l2 = tau(2,3)*0.236*d3
  l3 = tau(2,3)*d3
  l11tt = tau(2,3)*d0 + tau(1,3)
  l11th = tau(1,3)*c1(zhi,0)
  l12th = tau(1,3)*c2(zhi,0)

  return 
end subroutine psflux 
