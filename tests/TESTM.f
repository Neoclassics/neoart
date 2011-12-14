      PROGRAM TESTM

      IMPLICIT NONE

      INTEGER I,J
      REAL*8 MM, NN, TA, TB, MA, MB, V 
      DIMENSION MM(3,3), NN(3,3)

      TA = 1.
      TB = 0.5
      MA = 1.6727E-27
      MB = 12*MA

      CALL MENN(TA,TB,MA,MB,MM,NN)

      WRITE(*,*)'THE MM_(AB) COEFFICIENTS'
      DO 200 I = 1, 3
        WRITE(*,100)(MM(I,J), J = 1, 3)
 200  CONTINUE

      WRITE(*,*)'THE NN_(AB) COEFFICIENTS'
      DO 300 I = 1, 3
        WRITE(*,100)(NN(I,J), J = 1, 3)
 300  CONTINUE

      CALL MENN(TB,TA,MB,MA,MM,NN)

      WRITE(*,*)'THE NN_(BA) COEFFICIENTS'
      WRITE(*,*)'PRINTED IN MIRROR ODER AND NORMALIZED'
      V = SQRT(MB/MA)*EXP(1.5*LOG(TA/TB))
      DO 500 I = 1, 3
        WRITE(*,100)(V*NN(J,I), J = 1, 3)
 500  CONTINUE

 100  FORMAT(3(2X,1PE13.5))

      STOP
      END
