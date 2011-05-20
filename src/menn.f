
      SUBROUTINE MENN(TA,TB,MA,MB,MM,NN)
C--------------------------------------------------------------------
C     SUBROUTINE THAT CALCULATES THE FRICTION COEFFICIENTS M^AB AND 
C     N^AB 
C
C     INPUT    TA   : THE TEMPERATURE OF SPECIES A
C              TB   : THE TEMPERATURE OF SPECIES B
C              MA   : THE MASS OF SPECIES A
C              MB   : THE MASS OF SPECIES B
C     OUTPUT   MM   : ARRAY(3,3) CONTAINING THE TEST PARTICLE FRICTION
C                     COEFFICIENTS
C              NN   : ARRAY(3,3) CONTAINING THE FIELD PARTICLE FRICTION
C                     COEFFICIENTS
C--------------------------------------------------------------------

      IMPLICIT NONE

      REAL*8 TA, TB, MA, MB, MM, NN, XAB, TAB, MAB
      DIMENSION MM(3,3), NN(3,3)

      TAB = TA / TB
      MAB = MA / MB
      XAB = MAB / TAB

      MM(1,1) = - (1. + MAB) / (1 + XAB)**1.5
      MM(1,2) =   1.5*(1 + MAB)/(1 + XAB)**2.5
      MM(2,1) =   MM(1,2)
      MM(2,2) = - (3.25 + 4*XAB + 7.5*XAB**2) / (1+ XAB)**2.5
      MM(1,3) = - 15. *(1.+MAB)/(8.*(1+XAB)**3.5)
      MM(3,1) =   MM(1,3)
      MM(2,3) =   (69./16.+6.*XAB + 15.75*XAB**2)/(1+XAB)**3.5
      MM(3,2) =   MM(2,3)
      MM(3,3) = - (433./64.+ 17.*XAB + 459.*XAB**2/8. + 28.*XAB**3 +
     +            175.*XAB**4/8.)/(1 + XAB)**4.5

      NN(1,1) = - MM(1,1)
      NN(2,1) = - MM(2,1)
      NN(3,1) = - MM(3,1)
      NN(1,2) = - XAB*MM(1,2)
      NN(1,3) = - XAB**2*MM(1,3)
      NN(2,2) = 6.75*SQRT(MAB*XAB)/(1+XAB)**2.5
      NN(2,3) = -225.*TAB*XAB**2/(16.*(1+XAB)**3.5)
      NN(3,2) = NN(2,3)/(XAB*TAB)
      NN(3,3) = 2625.*SQRT(TAB)*XAB**2/(64.*(1+XAB)**4.5)

      RETURN 
      END
