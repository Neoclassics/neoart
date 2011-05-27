

      SUBROUTINE PENQ(TA,TB,MA,MB,PP,QQ)
C--------------------------------------------------------------------
C     SUBROUTINE THAT CALCULATES THE COUPLING COEFFICIENTS P^AB AND 
C     Q^AB 
C
C     INPUT    TA   : THE TEMPERATURE OF SPECIES A
C              TB   : THE TEMPERATURE OF SPECIES B
C              MA   : THE MASS OF SPECIES A
C              MB   : THE MASS OF SPECIES B
C     OUTPUT   PP   : ARRAY(2,2) CONTAINING THE COUPLING COEF. P
C              QQ   : ARRAY(2,2) CONTAINING THE COUPLING COEF. Q
C--------------------------------------------------------------------

      IMPLICIT NONE

      REAL TA, TB, MA, MB, PP, QQ, XAB, TAB, MAB
      DIMENSION PP(2,2), QQ(2,2)

      TAB = TA / TB
      MAB = MA / MB
      XAB = MAB / TAB

      PP(1,1) = - 3.*MAB / (1+XAB)**1.5
      PP(1,2) = -4.5 * XAB / (1+XAB)**2.5
      PP(2,1) = PP(1,2)
      PP(2,2) = - 3 * XAB*(13./4.+ 2*XAB + 2.5*XAB**2)/
     +            (1+XAB)**3.5

      QQ(1,1) = - PP(1,1) / TAB
      QQ(1,2) = - XAB*PP(1,2)
      QQ(2,1) = - PP(1,2)
      QQ(2,2) = 45. * XAB**2 / (4.*(1+XAB)**3.5)

      RETURN
      END
