
      SUBROUTINE TEST9
C-----------------------------------------------------------------
C     CALCULATE THE HEAT TRANSPORT OF A PLASMA CONTAINING
C     OF HYDROGEN AND OXYGEN, AND COMPARE THEM WITH THE ANALYTIC 
C     EXPRESSIONS OF HIRSHMANN AND SIGMAR.
C     ALSO ENERGY COUPLING IS CONSIDERED HERE.
C-----------------------------------------------------------------
      IMPLICIT NONE

      INTEGER NS,NC,NAR,ISEL,NREG,NLEG,NENERGY,NCOF,
     +        IC,NZM,I,J,NMAXGR,K, L, ISHOT
      REAL M,T,DEN,DS,CFF4,XI,TAU,K22,C3,ZSP, EPARR,
     +       EPS, SIGMA, NORM, ALPHA, OMTH, OMTI, RESUL
      REAL RHO, RN, E, Q, BN
      LOGICAL NEOGEO, NEOFRC

      include 'elem_config.inc'
      
      parameter(NAR = NELMAX+2)
      parameter(NZM = NIONMAX)
      PARAMETER(NMAXGR = 1000)

      DIMENSION NC(NAR),ZSP(NAR,NZM),M(NAR),T(NAR),DEN(NAR,NZM),
     +          DS(NAR,NZM,2),CFF4(NAR,NZM,4),XI(NAR,NZM),
     +          TAU(NAR,NAR),SIGMA(4),RESUL(40,10)


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
C     ION-ELECTRON COLLISIONS
      NCOF = 1
C     IN THE FIRST CALL THE MATRICES HAVE TO BE CALCULATED
      NEOFRC = .FALSE.
C     RECALCULATE THE GEOMETRY 
      NEOGEO = .TRUE. 
C     NO SHOT NUMBER 
      ISHOT = 0 
C     ZERO PARALLEL ELECTRIC FIELD 
      EPARR = 0. 
C     CALCULATE THE PS CONTRIBUTION
      IC = 2
C     SET THE COUPLING 
      SIGMA(1) = 0 
      SIGMA(2) = 0
      SIGMA(3) = 1
      SIGMA(4) = 1

C     THE NUMBER OF SPECIES IS 2
      NS = 2
C     IONS HAVE ONLY ONE CHARGE STATE 
      NC(1) = 1
      NC(2) = 1
C     THE MASS OF THE HYDROGEN AND OXYGEN
      M(1) = 1.6727E-27
      M(2) = 16*1.6727E-27
C     THE CHARGE OF THE HYDROGEN AND OXYGEN
      ZSP(1,1) = 1
      ZSP(2,1) = 4
C     THE DENSITY OF THE SPECIES IN 10^19 M^-3
      DEN(1,1) = 1.
      DEN(2,1) = DEN(1,1)/ZSP(2,1)**2

C     DO WITH AND WITHOUT ENERGY COUPLING
      DO 10000 K = 1, 2
      
      IF (K.EQ.2) SIGMA(3) = 0.
      
C     LOOP OVER OMTI
      DO 20000 L = 1, 40
 
      T(1) = 0.25 / SQRT(REAL(L))
      T(2) = T(1)

      IF (L.EQ.1) THEN
        T(1) = 1.
        T(2) = 1.
      ENDIF
      IF (L.EQ.2) THEN
        T(1) = 0.25
        T(2) = 0.25
      ENDIF
      IF (L.EQ.3) THEN
        T(1) = 0.18
        T(2) = 0.18
      ENDIF

      NEOFRC = .FALSE.
C     THE THERMODYNAMIC FORCES
      DO 204 I = 1, NS
        DO 204 J = 1, NC(I)
          DS(I,J,1) = 0.
          DS(I,J,2) = 0.
 204  CONTINUE
      DS(2,1,2) = 1.  
      DS(1,1,2) = 1.


      CALL NEOART(NS,NC,NAR,NZM,ZSP,M,T,DEN,DS,RHO,EPS,
     +            ISEL,ISHOT,NREG,SIGMA,NLEG,NENERGY,NCOF,
     +            NEOGEO,NEOFRC,IC,EPARR,CFF4)


      CALL COLXI(NAR,NZM,NS,NC,ZSP,DEN,T,M,TAU,XI)


      NORM = 1.6E-22*BN**2/(2.*Q**2*T(1))
      ALPHA  = DEN(2,1)*ZSP(2,1)**2/(DEN(1,1)*ZSP(1,1)**2)


      OMTH = SQRT(2*1.6E-16*T(1)*M(1))*DEN(1,1)*1E19/TAU(1,1)
     +       /(Q*RN)
      OMTI = SQRT(2*1.6E-16*T(2)*M(2))*DEN(2,1)*1E19/TAU(2,2)
     +       /(Q*RN)

C     USE 1/OMTI AS X COORDINATE
      RESUL(L,1) = 1./OMTI

C     CODE RESULT
      RESUL(L,1+K) = -ZSP(2,1)**2*CFF4(2,1,2)*NORM/(TAU(2,2)*OMTI)

C     RESULT OF HIRSHMANN SIGMAR
      RESUL(L,4) = K22(ALPHA,OMTH,OMTI,M(1),M(2),ZSP(1,1),
     +                 ZSP(2,1))/OMTI

      RESUL(L,5) = C3(0.D0,OMTI)/OMTI

20000 CONTINUE 
10000 CONTINUE

C     OPEN(11,FILE = 'TEST9.DAT')
C     THE COLLUMS HAVE THE FOLLOWING MEANING
C     COLLUM  1. : 1/ OMEGA_I TAU_II
C     COLLUM  2. : CODE RESULT INCLUDING ENERGY COUPLING
C     COLLUM  3. : CODE RESULT WITHOUT ENERGY COUPLING
C     COLLUM  4. : ANALYTIC RESULT WITH ENERGY COUPLING
C     COLLUM  5. : ANALYTIC RESULT WITHOUT ENERGY COUPLING
      DO 30000 I = 1, 40
        WRITE(*,1000)(RESUL(I,L),L=1,5)
30000 CONTINUE

 1000 FORMAT(14(1X,1PE13.5))

      RETURN 
      END
