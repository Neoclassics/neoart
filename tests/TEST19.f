      SUBROUTINE TEST19
C-----------------------------------------------------------------
C     Based on TEST1
C     For a circular equilibrium, compares the results obtained with
C     the metric read from CHEASE to the analytical formula
C-----------------------------------------------------------------
      IMPLICIT NONE

      INTEGER NS,NC,NAR,ISEL,NREG,NLEG,NENERGY,NCOF,
     +        IC,NZM,I,J,NMAXGR,ISHOT
      REAL M,T,DEN,DS,CFF1,CFF2,CFF3,CFF4,XI,TAU,COEFF,
     +       EPS, SIGMA, NORM, ZSP
      REAL RHO, RN, E, Q, BN, EPARR
      LOGICAL NEOGEO, NEOFRC

      include 'elem_config.inc'
      
      parameter(NAR = NELMAX+2)
      parameter(NZM = NIONMAX)
      PARAMETER(NMAXGR = 1000)

      DIMENSION NC(NAR),ZSP(NAR,NZM),M(NAR),T(NAR),DEN(NAR,NZM),
     +          DS(NAR,NZM,2),CFF1(NAR,NZM,4),CFF2(NAR,NZM,4),
     +          CFF3(NAR,NZM,4),CFF4(NAR,NZM,4),XI(NAR,NZM),
     +          TAU(NAR,NAR),SIGMA(4)


C     PARAMETERS FOR THE CIRCULAR GEOMETRY
C     USED HERE ONLY TO NORMALISE THE OUTPUT
C     SET THE PARAMETERS FOR THE CIRCULAR GEOMETRY
      RHO = 0.0485 
      E = 0.0485
      Q = 1.1328
      RN = 1.3379
      BN = 3.1074
C     READS THE GEOMETRY FROM THE CHEASE OUTPUT FILE
      RHO = 0.0485 ! RHO =  amin/RN
      ISEL = 4
C     set the electric field to zero
      EPARR = 0
C     SET THE ACCURACY
      EPS = 1E-5
c     FORCE THE BANANA REGIME IN THE CALCULATION OF THE TRANSPORT
C     FLUXES (IN THE CALCULATION OF THE VISCOSITY)
      NREG = 1
c     THE NUMBER OF LEGENDRE HARMONICS
      NLEG = 3
C     USE ENERGY SCATTERING IN THE CALC. OF VISCOSITY
      NENERGY = 1
C     SWITCH OFF ION-ELECTRON COLLISIONS
      NCOF = 0
C     IN THE FIRST CALL THE MATRICES HAVE TO BE CALCULATED
      NEOFRC = .FALSE.
C     CALCULATE THE BANANA PLATEAU CONTRIBUTION
      IC = 1

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
      DS(1,1,1) = 1.

C     SET EPARALLEL TO 0 
      EPARR = 0.

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

      NORM = 1.6E-22*BN**2*SQRT(E**3)/(2*(Q)**2
     +       *T(2)*TAU(2,2))

      WRITE(*,1000)'THE ION HEAT FLUX                 '
      WRITE(*,1000)'ION TEMPERATURE GRADIENT          ',
     +    SQRT(2.)*CFF4(2,1,2)*NORM
      WRITE(*,1000)'ELECTRON TEMPERATURE GRADIENT     ',
     +    SQRT(2.)*CFF3(2,1,2)*NORM
      WRITE(*,1000)'ION PRESSURE GRADIENT             ',
     +    SQRT(2.)*CFF2(2,1,2)*NORM
      WRITE(*,1000)'ELECTRON PRESSURE GRADIENT        ',
     +    SQRT(2.)*CFF1(2,1,2)*NORM

      NORM = 1.6E-22*BN**2*SQRT(E**3)/(2*(Q)**2
     +       *T(1)*TAU(1,2))

      WRITE(*,*)
      WRITE(*,1000)'THE ELECTRON HEAT FLUX            '
      WRITE(*,1000)'ION TEMPERATURE GRADIENT          ',
     +   CFF4(1,1,2)*NORM
      WRITE(*,1000)'ELECTRON TEMPERATURE GRADIENT     ',
     +   CFF3(1,1,2)*NORM
      WRITE(*,1000)'ION PRESSURE GRADIENT             ',
     +   CFF2(1,1,2)*NORM
      WRITE(*,1000)'ELECTRON PRESSURE GRADIENT        ',
     +   CFF1(1,1,2)*NORM

      WRITE(*,*)
      WRITE(*,1000)'THE PARTICLE FLUX                 '
      WRITE(*,1000)'ION TEMPERATURE GRADIENT          ',
     +    CFF4(1,1,1)*NORM
      WRITE(*,1000)'ELECTRON TEMPERATURE GRADIENT     ',
     +    CFF3(1,1,1)*NORM
      WRITE(*,1000)'ION PRESSURE GRADIENT             ',
     +    CFF2(1,1,1)*NORM
      WRITE(*,1000)'ELECTRON PRESSURE GRADIENT        ',
     +    CFF1(1,1,1)*NORM

 1000 FORMAT(A34,1X,1PE13.5,1X,A19)

      NORM = SQRT(E)*BN/(1600.*Q*DEN(1,1)*T(1))
      WRITE(*,*)
      WRITE(*,1000)'THE BOOTSTRAP CURRENT             '

      COEFF = 0.
      DO 211 I = 1, NS
        DO 211 J = 1, NC(I)
          COEFF = COEFF + CFF4(I,J,3)
 211  CONTINUE
      WRITE(*,1000)'ION TEMPERATURE GRADIENT          ',
     +    COEFF*NORM

      COEFF = 0.
      DO 212 I = 1, NS
        DO 212 J = 1, NC(I)
          COEFF = COEFF + CFF3(I,J,3)
 212  CONTINUE
      WRITE(*,1000)'ELECTRON TEMPERATURE GRADIENT     ',
     +    COEFF*NORM

      COEFF = 0.
      DO 213 I = 1, NS
        DO 213 J = 1, NC(I)
          COEFF = COEFF + CFF2(I,J,3)
 213  CONTINUE
      WRITE(*,1000)'ION PRESSURE GRADIENT             ',
     +    COEFF*NORM

      COEFF = 0.
      DO 214 I = 1, NS
        DO 214 J = 1, NC(I)
          COEFF = COEFF + CFF1(I,J,3)
 214  CONTINUE
      WRITE(*,1000)'ELECTRON PRESSURE GRADIENT        ',
     +    COEFF*NORM

      NORM = 1.E-3*BN**2*E/(Q*T(2))
      WRITE(*,*)
      WRITE(*,1000)'THE ION ROTATION                  '
      WRITE(*,1000)'ION TEMPERATURE GRADIENT          ',
     +    CFF4(2,1,4)*NORM
      WRITE(*,1000)'ELECTRON TEMPERATURE GRADIENT     ',
     +    CFF3(2,1,4)*NORM
      WRITE(*,1000)'ION PRESSURE GRADIENT              ',
     +    CFF2(2,1,4)*NORM
      WRITE(*,1000)'ELECTRON PRESSURE GRADIENT         ',
     +    CFF1(2,1,4)*NORM


      RETURN
      END      
