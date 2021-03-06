
      SUBROUTINE TEST4
C-----------------------------------------------------------------
C     CALCULATE THE TRANSPORT COEFFICIENTS OF AN PLASMA CONTAINING
C     OF HYDROGEN AND OXYGEN, AND COMPARE THEM WITH THE ANALYTIC 
C     EXPRESSIONS OF HIRSHMANN AND SIGMAR.
C-----------------------------------------------------------------
      IMPLICIT NONE

      INTEGER NS,NC,NAR,ISEL,NREG,NLEG,NENERGY,NCOF,
     +        IC,NZM,I,J,NMAXGR,K,L,ISHOT
      REAL M,T,DEN,DS,CFF1,CFF2,CFF3,CFF4,XI,TAU,C1,C2,C3,
     +       EPS, SIGMA, NORM, ALPHA, OMTH, OMTI, RESUL,
     +       RESUA,ZSP,EPARR
      REAL RHO, RN, E, Q, BN
      LOGICAL NEOGEO, NEOFRC

      include 'elem_config.inc'
      
      parameter(NAR = NELMAX+2)
      parameter(NZM = NIONMAX)
      PARAMETER(NMAXGR = 1000)

      DIMENSION NC(NAR),ZSP(NAR,NZM),M(NAR),T(NAR),DEN(NAR,NZM),
     +          DS(NAR,NZM,2),CFF1(NAR,NZM,4),CFF2(NAR,NZM,4),
     +          CFF3(NAR,NZM,4),CFF4(NAR,NZM,4),XI(NAR,NZM),
     +          TAU(NAR,NAR),SIGMA(4),RESUL(20,10),RESUA(20,4)


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
C     RECALCULATE THE GEOMETRY COEFFICIENTS 
      NEOGEO = .TRUE. 
C     NO SHOT NUMBER 
      ISHOT = 0
C     ZERO PARALLEL ELECTRIC FIELD 
      EPARR = 0 
C     CALCULATE THE PS CONTRIBUTION
      IC = 2
C     SET THE COUPLING
      SIGMA(1) = 0 
      SIGMA(2) = 0
      SIGMA(3) = 0
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
      DEN(2,1) = 2.*DEN(1,1)/ZSP(2,1)**2 

C     DO TWICE WITH DIFFERENT MASS RATIO
      DO 10000 K = 1, 2

      IF (K.EQ.2) M(2) = 1E6*M(1)

C     LOOP OVER TEMPERATURE
      DO 20000 L = 1, 20
        
      T(1) = 1E-2*EXP(2*(L-1)*LOG(10.)/19)
      T(2) = T(1)

      NEOFRC = .FALSE.
C     THE THERMODYNAMIC FORCES
      DO 201 I = 1, NS
        DO 201 J = 1, NC(I)
          DS(I,J,1) = 0.
          DS(I,J,2) = 0.
 201  CONTINUE
      DS(1,1,1) = 1.

      CALL NEOART(NS,NC,NAR,NZM,ZSP,M,T,DEN,DS,RHO,EPS,
     +            ISEL,ISHOT,NREG,SIGMA,NLEG,NENERGY,NCOF,
     +            NEOGEO,NEOFRC,IC,EPARR,CFF1)

C     AFTER THE FIRST CALL IT IS NOT NECESSARY TO RECALCULTE
C     THE MATRICES
      NEOFRC = .TRUE.

      DO 202 I = 1, NS
        DO 202 J = 1, NC(I)
          DS(I,J,1) = 0.
          DS(I,J,2) = 0.
 202  CONTINUE
      DS(2,1,1) = 1.

      CALL NEOART(NS,NC,NAR,NZM,ZSP,M,T,DEN,DS,RHO,EPS,
     +            ISEL,ISHOT,NREG,SIGMA,NLEG,NENERGY,NCOF,
     +            NEOGEO,NEOFRC,IC,EPARR,CFF2)

      DO 203 I = 1, NS
        DO 203 J = 1, NC(I)
          DS(I,J,1) = 0.
          DS(I,J,2) = 0.
 203  CONTINUE
      DS(1,1,2) = 1.

      CALL NEOART(NS,NC,NAR,NZM,ZSP,M,T,DEN,DS,RHO,EPS,
     +            ISEL,ISHOT,NREG,SIGMA,NLEG,NENERGY,NCOF,
     +            NEOGEO,NEOFRC,IC,EPARR,CFF3)

      DO 204 I = 1, NS
        DO 204 J = 1, NC(I)
          DS(I,J,1) = 0.
          DS(I,J,2) = 0.
 204  CONTINUE
      DS(2,1,2) = 1.


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


      RESUL(L,1) = 1./OMTH
      RESUL(L,2) = T(1)
C     SHOULD BE EQUAL TO C1(ALPHA,OMTH)      
      RESUL(L,4*K-1) = -ZSP(1,1)**2*CFF1(1,1,1)*NORM/TAU(1,2)
C     SHOULD BE EQUAL TO C2(ALPHA,OMTH)
      RESUL(L,4*K) = ZSP(1,1)**2*CFF3(1,1,1)*NORM/TAU(1,2)
C     SHOULD BE EQUAL TO C3(ALPHA,OMTH)
      RESUL(L,4*K+1) = -ZSP(1,1)**2*CFF3(1,1,2)*NORM/TAU(1,1)
C     SHOULD BE EQUAL TO C3(0,OMTI)
      RESUL(L,4*K+2) = -ZSP(2,1)**2*CFF4(2,1,2)*NORM/TAU(2,2)

      RESUA(L,1) = C1(ALPHA,OMTH)
      RESUA(L,2) = C2(ALPHA,OMTH)
      RESUA(L,3) = C3(ALPHA,OMTH)
      RESUA(L,4) = C3(0.D0,OMTI)

20000 CONTINUE
10000 CONTINUE

C     OPEN(11,FILE = 'TEST4.DAT')
C     THE COLLUMS OF THE FILE CONTAIN THE FOLLOWING RESULTS
C     COLLUM  1. : 1/OMEGA_H TAU_HH
C     COLLUM  2. : TEMPERATRUE (KEV)
C     COLLUM  3. : CODE RESULT OXYGEN EQUIV OF C1(ALPHA,OMTH)
C     COLLUM  4. : CODE RESULT OXYGEN EQUIV OF C2(ALPHA,OMTH)
C     COLLUM  5. : CODE RESULT OXYGEN EQUIV OF C3(ALPHA,OMTH)
C     COLLUM  6. : CODE RESULT OXYGEN EQUIV OF C3(0,OMTI)
C     COLLUM  7. : CODE RESULT HEAVY IMP EQUIV OF COLLUM 3
C     COLLUM  8. : CODE RESULT HEAVY IMP EQUIV OF COLLUM 4
C     COLLUM  9. : CODE RESULT HEAVY IMP EQUIV OF COLLUM 5
C     COLLUM 10. : CODE RESULT HEAVY IMP EQUIV OF COLLUM 6
C     COLLUM 11. : ANALYTIC FUNCTION C1(ALPHA,OMTH)
C     COLLUM 12. : ANALYTIC FUNCTION C2(ALPHA,OMTH)
C     COLLUM 13. : ANALYTIC FUNCTION C3(ALPHA,OMTH)
C     COLLUM 14. : ANALYTIC FUNCTION C3(0,OMTI)
      DO 30000 I = 1, 20
        WRITE(*,1000)(RESUL(I,L),L=1,10),(RESUA(I,L),L=1,4)
30000 CONTINUE

 1000 FORMAT(14(1X,1PE13.5))


      RETURN 
      END 
