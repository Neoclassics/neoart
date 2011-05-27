

      SUBROUTINE TEST17
C-----------------------------------------------------------------
C     TEST THE HAMADA COORDINATE VERSION OF THE GEOMETRY DEPEN-
C     DENT QUANTITIES
C-----------------------------------------------------------------
      IMPLICIT NONE

      INTEGER NMG
      PARAMETER(NMG = 1000)

      INTEGER ISEL, MMX, I, ISHOT, J
      REAL RHO, E, Q, BN, RN, EPS, BAV1, B2AV1, BI2A1, RBT1, 
     +       BGRADP1, DPSIDR1, RNQ1, FC1, FM1, BAV2, B2AV2, 
     +       BI2A2, RBT2, BGRADP2, DPSIDR2, RNQ2, FC2, FM2,
     +       R2I
      DIMENSION FM1(NMG), FM2(NMG)

      ishot = 1001

      do 200 j = 1, 10

C     THE GEOMETRY IS DETERMINED BY THE FOLLOWING QUANTITIES
C     RHO IS THE INVERSE ASPECT RATIO 
      E = 3e-03 + (j-1)*0.5 / 16.5
C     Q IS THE SAFETY FACTOR
      Q = 1.4
C     BN IS THE TOTAL MAGNETIC FIELD AT THE LOCATION CHI = PI/2 
      BN = 2.5 * sqrt((1 + e**2 * (1/q**2 -1))/(1-e**2))
C     RN IS THE MAJOR RADIUS OF THE MAGNETIC AXIS
      RN = 1.65
C     COPY THEM INTO THE VALUES USED BY THE CODE
      CALL CIRCGEOM(1,RHO,RN,E,Q,BN)

      rho = e * 1.65 

      ISEL = 3
      EPS = 1E-6

      CALL GEOM(NMG,ISEL,ISHOT,RHO,EPS,BAV1,B2AV1,BI2A1,RBT1,BGRADP1,
     +          DPSIDR1,RNQ1,FC1,FM1,MMX,R2I)
      
      ISEL = 2
      EPS = 1E-6

      CALL GEOM(NMG,ISEL,ISHOT,RHO,EPS,BAV2,B2AV2,BI2A2,RBT2,BGRADP2,
     +          DPSIDR2,RNQ2,FC2,FM2,MMX,R2I)
      
      write(*,*)'Epsilon ',e
      WRITE(*,1000)'BAV    ',BAV1, BAV2
      WRITE(*,1000)'B2AV   ',B2AV1, B2AV2
      WRITE(*,1000)'BI2A   ',BI2A1, BI2A2
      WRITE(*,1000)'RBT    ',RBT1, RBT2
      WRITE(*,1000)'BGRADP ',BGRADP1, BGRADP2
      WRITE(*,1000)'DPSIDR ',DPSIDR1, DPSIDR2
      WRITE(*,1000)'RNQ    ',RNQ1, RNQ2
      WRITE(*,1000)'FC     ',FC1, FC2
      WRITE(*,1000)'FM     '
      DO 100 I = 1, 3 
        WRITE(*,1100)I,FM1(I),FM2(I)
 100  CONTINUE
      write(*,*)

200   continue
 1000 FORMAT(A7,2(1X,1PE13.5))
 1100 FORMAT(I3,2(1X,1PE13.5))

      RETURN
      END
