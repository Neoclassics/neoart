
      SUBROUTINE TEST13
C-----------------------------------------------------------------
C     COMPARE THE 2 AND 3 LAGUERRE POLYNOMAL EXPANSIONS OF THE 
C     ION HEAT FLUX AT FINITE INVERSE ASPECT RATIO.
C-----------------------------------------------------------------
      IMPLICIT NONE

      INTEGER NS,NC,NAR,ISEL,NREG,NLEG,NENERGY,NCOF,
     +        IC,NZM,I,J,NMAXGR,L,NMG, mmx, ISHOT 
      REAL M,T,DEN,DS,CFF4,XI,TAU, EPS, SIGMA, NORM, RESUL
      REAL RHO, RN, E, Q, BN, BAV, B2AV, BI2A, RBT, BGRADP,
     +       DPSIDR, RNQ, FC, FM, DUM, ZSP, EPARR, R2I, GCLASS
      LOGICAL NEOGEO, NEOFRC

      include 'elem_config.inc'
      
      parameter(NAR = NELMAX+2)
      parameter(NZM = NIONMAX)
      PARAMETER(NMAXGR = 1000)

      DIMENSION NC(NAR),ZSP(NAR,NZM),M(NAR),T(NAR),DEN(NAR,NZM),
     +          DS(NAR,NZM,2),CFF4(NAR,NZM,4),XI(NAR,NZM),
     +          TAU(NAR,NAR),SIGMA(4),RESUL(22,6),FM(NMAXGR)


C     USE THE CIRCULAR GEOMETRY APPROXIMATION
      ISEL = 2
C     SET THE ACCURACY
      EPS = 1E-5
c     FORCE THE BANANA REGIME IN THE CALCULATION OF THE TRANSPORT
C     FLUXES (IN THE CALCULATION OF THE VISCOSITY)
      NREG = 1 
c     THE NUMBER OF LEGENDRE HARMONICS
      NLEG = 2
C     USE ENERGY SCATTERING IN THE CALC. OF VISCOSITY
      NENERGY = 1
C     SWITCH OFF ION-ELECTRON COLLISIONS
      NCOF = 0
C     IN THE FIRST CALL THE MATRICES HAVE TO BE CALCULATED
      NEOFRC = .FALSE.
C     RECALCULATE THE GEOMETRY DEPENDENT PARAMETERS EVERY TIME
      NEOGEO = .TRUE.
C     NO SHOT NUMBER 
      ISHOT = 0 
C     ZERO PARALLEL ELECTRIC FIELD 
      EPARR = 0. 

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
      DEN(1,1) = 1.
      DEN(2,1) = 1. 
C     THE TEMPERATURE IN KEV
      T(1) = 10.
      T(2) = 10.
C     SET THE COUPLING
      SIGMA(1) = 0
      SIGMA(2) = 0
      SIGMA(3) = 0 
      SIGMA(4) = 0

C     LOOP OVER INVERSE ASPECT RATIO
      DO 10000 L = 1, 22

C     SET THE PARAMETERS FOR THE CIRCULAR GEOMETRY
      RHO = 0.05
      E = 1E-4 + 0.4999*(L-2.)/20.
      IF (L.EQ.1) E = 1E-4
      IF (L.EQ.2) E = 1E-2
      Q = 2
      RN = 1.65
      BN = 2.5
C     COPY THEM INTO THE VALUES USED BY THE CODE
      CALL CIRCGEOM(1,RHO,RN,E,Q,BN)


C     CALCULATE THE TOTAL FLUX WITH 2 LEGENDRE HARM.
      IC = 3
      NLEG = 2

      DO 304 I = 1, NS
        DO 304 J = 1, NC(I)
          DS(I,J,1) = 0.
          DS(I,J,2) = 0.
 304  CONTINUE
      DS(2,1,2) = 1.


      CALL NEOART(NS,NC,NAR,NZM,ZSP,M,T,DEN,DS,RHO,EPS,
     +            ISEL,ISHOT,NREG,SIGMA,NLEG,NENERGY,NCOF,
     +            NEOGEO,NEOFRC,IC,EPARR,CFF4)


      CALL COLXI(NAR,NZM,NS,NC,ZSP,DEN,T,M,TAU,XI)

      NORM = 1.6E-22*BN**2*SQRT(E**3)/(2*(Q)**2*(1.-E**2)
     +       *T(2)*TAU(2,2))

C     USE INVERSE ASPECT RATIO AS X AXIS
      RESUL(L,1) = E
      RESUL(L,2) = -SQRT(2.)*CFF4(2,1,2)*NORM

      NMG = NMAXGR
      CALL GEOM(NMG,ISEL,ISHOT,RHO,EPS,BAV,B2AV,BI2A,RBT,BGRADP,
     +          DPSIDR,RNQ,FC,GCLASS,FM,MMX,R2I)

      DUM = (1.-FC)/FC

C     TWO TERM LAGUERRE
      RESUL(L,4) = (SQRT(1.-E**2)*0.46*DUM/(1.+ 0.46*DUM) +
     +              1 + 1.5*E**2 - SQRT(1.-E**2))/SQRT(E)
C     THREE TERM LAGUERRE
      RESUL(L,5) = (SQRT(1.-E**2)*(0.4615*DUM+0.1029*DUM**2)/
     +             (1.+0.7305*DUM+0.1029*DUM**2) +
     +              1 + 1.5*E**2 - SQRT(1.-E**2))/SQRT(E)
      RESUL(L,6) = 0.675

C     CALCULATE THE TOTAL CONTRIBUTION WITH NLEG = 3
      IC = 3
      NLEG = 3


C     THE THERMODYNAMIC FORCES
      DO 204 I = 1, NS
        DO 204 J = 1, NC(I)
          DS(I,J,1) = 0.
          DS(I,J,2) = 0.
 204  CONTINUE
      DS(2,1,2) = 1.


      CALL NEOART(NS,NC,NAR,NZM,ZSP,M,T,DEN,DS,RHO,EPS,
     +            ISEL,ISHOT,NREG,SIGMA,NLEG,NENERGY,NCOF,
     +            NEOGEO,NEOFRC,IC,EPARR,CFF4)

      RESUL(L,3) = -SQRT(2.)*CFF4(2,1,2)*NORM


10000 CONTINUE

C     OPEN(11, FILE = 'TEST13.DAT')
C     THE COLLUMS HAVE THE FOLLOWING MEANING 
C     COLLUM 1. : THE INVERSE ASPECT RATIO
C     COLLUM 2. : THE CODE RESULT FOR 2 LAGUERRE HARM.
C     COLLUM 3. : THE CODE RESULT FOR 3 LAGUERRE HARM.
C     COLLUM 4. : THE ANALYTIC RESULT FOR 2 LAGUERRE HARM.
C     COLLUM 5. : THE ANALYTIC RESULT FOR 3 LAGUERRE HARM.
      DO 30000 L = 1, 22
        WRITE(*,30001)(RESUL(L,I),I = 1,6)
30000 CONTINUE
30001 FORMAT(9(1X,1PE13.5))

      RETURN
      END
