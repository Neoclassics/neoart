
      SUBROUTINE RALF
C-----------------------------------------------------------------
C     THE ASDEX UPGRADE IMPURITY TRANSPORT EXPERIMENTS. 
C-----------------------------------------------------------------
      IMPLICIT NONE

      INTEGER NS,NC,NAR,ISEL,NZM,I,J,L,ISHOT
      REAL*8 M,T,DEN,DS,EPS, DD, VV, ZSP
      REAL*8 RHO, RN, E, Q, BN

      PARAMETER(NAR = 5)
      PARAMETER(NZM = 30)

      DIMENSION NC(NAR),ZSP(NAR,NZM),M(NAR),T(NAR),DEN(NAR,NZM),
     +          DS(NAR,NZM,2),DD(NAR,NZM),VV(NAR,NZM)


C     SET THE PARAMETERS FOR THE EQUILIBRIUM
      ISHOT = 10502
      E = 0.08
      Q = 1
      RN = 1.65
      BN = 2.5
      RHO = E*RN

C     COPY THEM INTO THE VALUES USED BY THE CODE
      CALL CIRCGEOM(1,RHO,RN,E,Q,BN)
C     READ THE EQUILIBRIUM VALUES FROM FILE
      ISEL = 3
C     SET THE ACCURACY
      EPS = 1E-4

C     THE NUMBER OF SPECIES IS 3
      NS = 3
C     IONS AND ELECTRONS HAVE ONLY ONE CHARGE 
      NC(1) = 1
      NC(2) = 1
      NC(3) = 1
C     THE MASS OF THE ELECTRON AND DEUTERIUM
      M(1) = 9.1096E-31
      M(2) = 2*1.6727E-27
C     THE CHARGE OF THE ELECTRON AND DEUTERIUM
      ZSP(1,1) = -1
      ZSP(2,1) = 1

C     FOUR CASES
      DO 100 L = 1, 4
  
      IF (L.EQ.1) THEN
        M(3) = 20*1.6727E-27
        ZSP(3,1) = 10
        DEN(1,1) = 8.92
        DEN(3,1) = 0.03
        DEN(2,1) = DEN(1,1) - ZSP(3,1)*DEN(3,1)
        T(1) = 3.2
        T(2) = 3.2
        T(3) = 3.2
      ENDIF   

      IF (L.EQ.2) THEN
        M(3) = 40*1.6727E-27
        ZSP(3,1) = 17
        DEN(1,1) = 8.44
        DEN(3,1) = 0.0043
        DEN(2,1) = DEN(1,1) - ZSP(3,1)*DEN(3,1)
        T(1) = 2.8
        T(2) = 2.8
        T(3) = 2.8
      ENDIF   

      IF (L.EQ.3) THEN
        M(3) = 84*1.6727E-27
        ZSP(3,1) = 30
        DEN(1,1) = 8.29
        DEN(3,1) = 0.00072
        DEN(2,1) = DEN(1,1) - ZSP(3,1)*DEN(3,1)
        T(1) = 2.7
        T(2) = 2.7
        T(3) = 2.7
      ENDIF   

      IF (L.EQ.4) THEN
        M(3) = 131*1.6727E-27
        ZSP(3,1) = 37
        DEN(1,1) = 8.24
        DEN(3,1) = 0.000099
        DEN(2,1) = DEN(1,1) - ZSP(3,1)*DEN(3,1)
        T(1) = 2.4
        T(2) = 2.4
        T(3) = 2.4
      ENDIF   

C     THE THERMODYNAMIC FORCES
      DO 201 I = 1, NS
        DO 201 J = 1, NC(I)
          DS(I,J,1) = 1.+ 2.5
          DS(I,J,2) = 2.5
 201  CONTINUE

      CALL DANDV(NAR,NZM,NS,NC,ZSP,ISEL,ISHOT,
     +           M,T,DEN,DS,RHO,EPS,DD,VV)

      WRITE(*,1000)ZSP(3,1),DD(3,1),VV(3,1),VV(3,1)/DD(3,1)

 100  CONTINUE

 1000 FORMAT(4(1X,1PE13.5))

      RETURN 
      END
