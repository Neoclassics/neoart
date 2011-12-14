

      SUBROUTINE TEST16
C-----------------------------------------------------------------
C     TEST THE DANDV ROUTINE. A PURE PLASMA IN THE BANANA 
C     PLATEAU REGIME 
C-----------------------------------------------------------------
      IMPLICIT NONE

      INTEGER NS,NC,NAR,ISEL,NZM,I,J,NMAXGR,IC,ID,ISHOT
      REAL M,T,DEN,DS,XI,TAU, EPS, NORM, DD, VV, ZSP
      REAL RHO, RN, E, Q, BN, VT, EPARR

      include 'elem_config.inc'
      
      parameter(NAR = NELMAX+2)
      parameter(NZM = NIONMAX)
      PARAMETER(NMAXGR = 1000)

      DIMENSION NC(NAR),ZSP(NAR,NZM),M(NAR),T(NAR),DEN(NAR,NZM),
     +          DS(NAR,NZM,2),XI(NAR,NZM),TAU(NAR,NAR),DD(NAR,NZM),
     +          VV(NAR,NZM),VT(NAR,NZM)


C     SET THE PARAMETERS FOR THE CIRCULAR GEOMETRY
      RHO = 0.05
      E = 1E-4
      Q = 2
      RN = 1.65
      BN = 2.5
C     COPY THEM INTO THE VALUES USED BY THE CODE
      CALL CIRCGEOM(1,RHO,RN,E,Q,BN)
C     USE THE CIRCULAR GEOMETRY APPROXIMATION
      ISEL = 2
C     SET THE ACCURACY
      EPS = 1E-5

C     THE NUMBER OF SPECIES IS 2
      NS = 2
C     IONS AND ELECTRONS HAVE ONLY ONE CHARGE 
      NC(1) = 1
      NC(2) = 1
C     THE MASS OF THE ELECTRON AND PROTON
      M(1) = 9.1096E-31
      M(2) = 1.6727E-27
C     THE CHARGE OF THE ELECTRON AND ION
      ZSP(1,1) = -1
      ZSP(2,1) = 1
C     THE DENSITY OF THE SPECIES IN 10^19 M^-3
      DEN(1,1) = 5.
      DEN(2,1) = 5. 
C     THE TEMPERATURE IN KEV
      T(1) = 10.
      T(2) = 10.

C     THE THERMODYNAMIC FORCES
      DO 201 I = 1, NS
        DO 201 J = 1, NC(I)
          DS(I,J,1) = 0.
          DS(I,J,2) = 0.
 201  CONTINUE

C     CALCULATE THE BP CONTRIBUTION ONLY
      IC = 1
C     NO VT CONTRIBUTION
      ID = 0

c     SET THE LOOP VOLTAGE
      EPARR = RN

      CALL DANDV(NAR,NZM,NS,NC,ZSP,ISEL,ISHOT,IC,ID,M,T,DEN,
     +           DS,RHO,EPARR,EPS,DD,VV,VT)

C     NOW CALCULATE THE NORMALIZED COLLISION FREQUENCIES
      CALL COLXI(NAR,NZM,NS,NC,ZSP,DEN,T,M,TAU,XI)


c     NORM = 1.6E-3*BN**2*DEN(1,1)*SQRT(E**3)/(2*(Q)**2
c    +       *T(1)*TAU(1,2))

      NORM = sqrt(E)*BN/Q
      DO 300 I = 1, NS
        WRITE(*,*)'D NORM ',DD(I,1)*NORM,' VV NORM ',VV(I,1)*NORM
 300  CONTINUE

      RETURN 
      END











