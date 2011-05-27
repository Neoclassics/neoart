      
      SUBROUTINE CHIS(NAR,NZM,NS,NC,ZSP,ISEL,FILE,IC, 
     +                M,T,DEN,DS,RHO,EPS,CHI)
C--------------------------------------------------------------------
C     SUBROUTINE THAT CALCULATES THE HEAT DIFFUSION COEFFICIENT
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
C                FILE  GEOEMETRY FILE, ONLY USED WHEN ISEL = 3
C                IC    THE CONTRIBUTION FOR WHICH THE COEFFICIENTS
C                      ARE CALCULATED.
C                      0 THE CLASSICAL PARTICLE FLUX            
C                      1 BANANA PLATEAU CONTRIBUTION
C                      2 PFIRSCH SCHLUETER CONTRIBUTION
C                      3 BANANA-PLATEAU, CLASSICAL AND P.S.
C                M     ARRAY(NS) THE MASS OF THE SPECIES IN KG
C                T     ARRAY(NS) THE TEMPERATURE OF THE SPECIES IN
C                      KEV
C                DEN   ARRAY(NAR,NZM) THE DENSITY OF EVERY COMPONENT
C                      IN UNITS OF 1.10^-19 M^-3
C                DS    ARRAY(NAR,NZM,2) THE THERMODYNAMIC FORCES 
C                      SEE THE COMMENTS OF THE ROUTINE NEOART TO 
C                      FIND THE EXACT DEFINITION
C                RHO   THE FLUX SURFACE LABEL     
C                EPS   THE DESIRED ACCURACY
C     OUTPUT     CHI   ARRAY(NAR,NZM) OF THE HEAT DIFFUSION COEFFICIENTS 
C                      OF EVERY COMPONENTS IN UNITS OF RHO^2/SEC    
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

      INTEGER NS, NC, ISEL, NREG, NLEG, NENERGY, NCOF, IC,
     +        I, J, K, L, MT, NAR, NZM
      REAL M, T, DEN, DS, DSL, RHO, EPS, SIGMA, COEFF, CHI, 
     +       ZSP, EPARN
      CHARACTER*(*) FILE
      LOGICAL NEOFRC, NEOGEO

      DIMENSION M(NSM), T(NSM), DEN(NSM,NCM), DS(NSM,NCM,2), 
     +          DSL(NSM,NCM,2), SIGMA(4), COEFF(NSM,NCM,4),
     +          CHI(NSM,NCM), NC(NSM), ZSP(NSM,NCM)
         
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

C     FIRST INITIALIZE CHI TO ZERO
      DO 1 I = 1, NS 
        DO 1 J = 1, NC(I)
          CHI(I,J) = 0.
 1    CONTINUE
      
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
        DSL(I,J,1) = 1.  !only temperature in pressure term
        DSL(I,J,2) = 1.  !temperature term

        CALL NEOART(NS,NC,NAR,NZM,ZSP,M,T,DEN,DSL,RHO,EPS,
     +              ISEL,FILE,NREG,SIGMA,NLEG,NENERGY,NCOF,
     +              NEOGEO,NEOFRC,IC,EPARN,COEFF)

C       THE FLUX OF COMPONENT I,J APPEARS IN THE HEAT DIFFUSION
C       COEFF(I,J,2) = Q(I,J)/T(I)        
        CHI(I,J) = -COEFF(I,J,2)/(DEN(I,J)*1E19)

C       DO NOT CALCULATE A NEW FRICTION AND VISCOSITY
        NEOGEO = .FALSE. 
        NEOFRC = .TRUE.

 100  CONTINUE
      
      RETURN 
      END

