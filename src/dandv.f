C--------------------------------------------------------------------
C     START OF THE ACTUAL PACKAGE
C--------------------------------------------------------------------

      SUBROUTINE DANDV(NAR,NZM,NS,NC,ZSP,ISEL,ISHOT, IC, ID, 
     +                 M,T,DEN,DS,RHO,EPARR,EPS,DD,VV,VT)
C--------------------------------------------------------------------
C     SUBROUTINE THAT CALCULATES THE DIFFUSION AND PINCH COEFFICIENTS
C
C     INPUT    : NAR   THE MAXIMUM NUMBER OF SPECIES (USED FOR 
C                      CONSISTENCY CHECK.
C                NZM   THE MAXIMUM NUMBER OF CHARGE STATES (USED
C                      FOR CONSISTENCY CHECK)
C                NS    THE ACTUAL NUMBER OF SPECIES
C                NC    ARRAY(NS) THE NUMBER OF CHARGES STATES PER
C                      SPECIES
C                ZSP   ARRAY(NAR,NZM) REAL.  THE CHARGE OF EVERY
C                      COMPONENTS (NORMALIZED TO E)
C                ISEL  INTEGER THAT SELECTS THE REPRESENTATION OF 
C                      THE GEOMETRY 1 = HAMADA COORDINATES 2 = 
C                      CIRCULAR APPROXIMATION, 3 = READ FROM FILE
C                ISHOT SHOT NUMBER, ONLY USED WHEN ISEL = 3
C                IC    THE CONTRIBUTION FOR WHICH THE COEFFICIENTS
C                      ARE CALCULATED. 
C                      1 BANANA PLATEAU CONTRIBUTION
C                      2 PFIRSCH SCHLUETER CONTRIBUTION
C                      3 BOTH BANANA-PLATEAU AND P.S.
C                ID    THE DRIFT CONTRIBUTION
C                      0 DO NOT GIVE TEMPERATURE DRIFT CONTRIBUTIONS VT
C                      1 GIVE TEMPERATURE DRIFT CONTRIBUTION VT
C                M     ARRAY(NS) THE MASS OF THE SPECIES IN KG
C                T     ARRAY(NS) THE TEMPERATURE OF THE SPECIES IN
C                      KEV
C                DEN   ARRAY(NAR,NZM) THE DENSITY OF EVERY COMPONENT
C                      IN UNITS OF 1.10^-19 M^-3
C                DS    ARRAY(NAR,NZM,2) THE THERMODYNAMIC FORCES 
C                      SEE THE COMMENTS OF THE ROUTINE NEOART TO 
C                      FIND THE EXACT DEFINITION
C                RHO   THE FLUX SURFACE LABEL
C                EPARR LOOP VOLTAGE 
C                EPS   THE DESIRED ACCURACY
C     OUTPUT     DD    ARRAY(NAR,NZM) THE DIFFUSION COEFFICIENTS 
C                      OF EVERY COMPONENTS IN UNITS OF RHO^2/SEC
C                VV    ARRAY(NAR,NZM) THE TOTAL PINCH COEFFICIENTS OF 
C                      EVERY COMPONENT IN UNITS OF RHO/SEC
C                VT    ARRAY(NAR,NZM) TEMPERATURE DEPENDENT PINCH      
C                      EVERY COMPONENT IN UNITS OF RHO/SEC
C                      (ONLY CALCULATED FOR ID=1)      
C
C     THE ROUTINE CALLS THE FOLLOWING SUBROUTINES
C     PERR    : FOR ERROR HANDLING
C     NEOART  : TO CALCULATE THE FLUXES
C--------------------------------------------------------------------
      IMPLICIT NONE

      include 'elem_config.inc'
      
      integer NSM, NCM 	
      parameter(NSM = NELMAX+2)
      parameter(NCM = NIONMAX)

      INTEGER NS, NC, ISEL, NREG, NLEG, NENERGY, NCOF, IC, ID,
     +        I, J, K, L, MT, NAR, NZM, ISHOT
      REAL M, T, DEN, DS, DSL, RHO, EPS, SIGMA, COEFF, DD, VV,
     +       ZSP, VT, EPARR, EPARN
      LOGICAL NEOFRC, NEOGEO

      DIMENSION M(NSM), T(NSM), DEN(NSM,NCM), DS(NSM,NCM,2), 
     +          DSL(NSM,NCM,2), SIGMA(4), COEFF(NSM,NCM,4),
     +          DD(NSM,NCM), VV(NSM,NCM), NC(NSM), ZSP(NSM,NCM),
     +          VT(NSM,NCM) 

C     CHECK CONSITENCY OF SETTINGS
      IF ((NAR.NE.NSM).OR.(NZM.NE.NCM)) CALL PERR(1)

C     ALWAYS USE ALL THE FULL EXPRESSION OF THE VISCOSITY
      NREG = 0
C     USE 3 LEGENDRE POLYNOMALS
      NLEG = 3
C     USE ENERGY SCATTERING IN THE COLLISION OPERATOR
      NENERGY = 1
C     USE ION - ELECTRON COLLISIONS
      NCOF = 1
C     USE ALL COUPLINGS IN THE PFIRSCH SCHLUETER REGIME
      SIGMA(1) = 1
      SIGMA(2) = 1
      SIGMA(3) = 1
      SIGMA(4) = 1
C     ON THE FIRST CALL RECACLULATE THE FRICTION/VISCOSITY AND
C     GEOMETRY PARAMETERS
      NEOFRC = .FALSE.
      NEOGEO = .TRUE. 

C     FIRST INITIALIZE D AND V TO ZERO
      DO 1 I = 1, NS 
        DO 1 J = 1, NC(I)
          DD(I,J) = 0.
          VV(I,J) = 0.
 1    CONTINUE

      IF (ID.EQ.1) THEN
         DO I=1,NS
            DO J=1,NC(I)
               VT(I,J) = 0.
            ENDDO
         ENDDO
      ENDIF

C     SET THE EPAR TO ZERO FOR THE TIME BEING
      EPARN = 0.
C     LOOP OVER THE SPECIES 
      DO 100 I = 1, NS
        DO 100 J = 1, NC(I)

C       DETERMINE D
        DO 200 K = 1, NS
          DO 200 L = 1, NC(K)
            DO 200 MT = 1, 2
              DSL(K,L,MT) = 0.
 200    CONTINUE
        DSL(I,J,1) = 1.

        CALL NEOART(NS,NC,NAR,NZM,ZSP,M,T,DEN,DSL,RHO,EPS,
     +              ISEL,ISHOT,NREG,SIGMA,NLEG,NENERGY,NCOF,
     +              NEOGEO,NEOFRC,IC,EPARN,COEFF)

C       THE FLUX OF COMPONENT I,J APPEARS IN THE DIFFUSION
        DD(I,J) = -COEFF(I,J,1)/(DEN(I,J)*1E19)

C       THE CONTRIBUTION TO ALL OTHER FLUXES APPEAR AS PINCH
           DO 300 K = 1, NS
              DO 300 L = 1, NC(K)
                 IF ((K.EQ.I).AND.(L.EQ.J)) GOTO 300
                 VV(K,L) = VV(K,L) + COEFF(K,L,1)*
     +                (DS(I,J,1)-DS(I,J,2))/(DEN(K,L)*1E19)
 300       CONTINUE

C       DO NOT CALCULATE A NEW FRICTION AND VISCOSITY
        NEOGEO = .FALSE. 
        NEOFRC = .TRUE.

C       NOW CALCULATE THE PINCH CONTRIBUTION DUE TO THE TEMPERATURE
C       GRADIENT.
        DO 400 K = 1, NS
          DO 400 L = 1, NC(K)
            DO 400 MT = 1, 2
              DSL(K,L,MT) = 0.
 400    CONTINUE
        DSL(I,J,1) = DS(I,J,2)
        DSL(I,J,2) = DS(I,J,2)

        CALL NEOART(NS,NC,NAR,NZM,ZSP,M,T,DEN,DSL,RHO,EPS,
     +              ISEL,ISHOT,NREG,SIGMA,NLEG,NENERGY,NCOF,
     +              NEOGEO,NEOFRC,IC,EPARN,COEFF)

C       ALL THESE CONTRIBUTIONS APPEAR IN THE PINCH
        DO 500 K = 1, NS
          DO 500 L = 1, NC(K)
            VV(K,L) = VV(K,L) + COEFF(K,L,1)/(DEN(K,L)*1E19)
 500    CONTINUE

C       THE TEMPERATURE DEPENDENT PART IS STORED SEPERATELY IN VT
        IF (ID.EQ.1) THEN
           DO  K = 1, NS
              DO  L = 1, NC(K)
                 VT(K,L) = VT(K,L) + COEFF(K,L,1)/(DEN(K,L)*1E19)
              ENDDO
           ENDDO
        ENDIF

 100  CONTINUE

C     ADD THE WARE PINCH
      DO 700 K = 1, NS
        DO 700 L = 1, NC(K)
          DO 700 MT = 1, 2
            DSL(K,L,MT) = 0.
 700  CONTINUE


      CALL NEOART(NS,NC,NAR,NZM,ZSP,M,T,DEN,DSL,RHO,EPS,
     +            ISEL,ISHOT,NREG,SIGMA,NLEG,NENERGY,NCOF,
     +            NEOGEO,NEOFRC,IC,EPARR,COEFF)

C     ADD TO THE PINCH
      DO 600 K = 1, NS
        DO 600 L = 1, NC(K)
          VV(K,L) = VV(K,L) + COEFF(K,L,1)/(DEN(K,L)*1E19)
 600  CONTINUE


      write(*,*)eparr,coeff(1,1,1)
      RETURN 
      END







