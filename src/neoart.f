      SUBROUTINE NEOART(NS,NC,NAR,NZM,ZSP,M,T,DEN,DS,RHO,EPS,
     +                  ISEL,ISHOT,NREG,SIGMA,NLEG,NENERGY,NCOF,
     +                  NEOGEO,NEOFRC,IC,EPARR,COEFF)
C--------------------------------------------------------------------
C     SUBROUTINE THIS ROUTINE CALCULATES THE TRANSPORT COEFFICIENTS
C     OF NEOCLASSICAL THEORY. 
C
C     INPUT   NS        : THE NUMBER OF SPECIES
C             NC        : ARRAY(NS) THAT GIVES THE NUMBER OF 
C                         CHARGE STATES PER SPECIES
C             NAR       : ACTUAL (FIRST) DIMENSION OF THE ZSP,DEN, 
C                         AND DS ARRAY, USED FOR CHECKING CONSISTENCY
C             NZM       : ACTUAL (SECOND) DIMENS. OF ARRAY (CONSIS.)
C             ZSP       : ARRAY(NAR,NZM) OF REAL THAT CONTAINS
C                         THE CARGE NUMBER OF EVERY COMPONENT
C             M         : ARRAY(NS) THAT CONTAINS THE MASS OF 
C                         EVERY SPECIES IN KG
C             T         : ARRAY(NS) THAT CONTAINS THE TEMPERATURE
C                         OF THE SPECIES IN KEV
C             DEN       : ARRAY(NAR,NZM) THAT CONTAINS THE DENSITY
C                         OF EVERY COMPONENT IN UNITS OF 1 10^19
C             DS        : ARRAY(NAR,NZM,2) THAT CONTAINS THE FIRST
C                         AND SECOND THERMODYNAMIC FORCE (PRESSURE
C                         GRADIENT AND TEMPERATURE GRADIENT)
C
C                                          D LN P_IJ 
C                              DS(I,J,1) = ---------   
C                                            D rho
C
C                                          D LN T_I
C                              DS(I,J,2) = --------
C                                           D rho
C
C                         WHERE RHO IS THE FLUX LABEL ACTUALLY USED
C                         THE FLUXES ARE ALSO CALCULATED IN THIS 
C                         COORDINATE. FOR THIS TO BE POSSIBLE THE 
C                         GEOMETRY ROUTINE MUST SUPPLY DPSIDR = 
C                         D PSI / D RHO WHERE PSI IS THE POLOIDAL 
C                         FLUX. 
C             RHO       : THE FLUX SURFACE LABEL
C             EPS       : REQUIRED ACCURACY.
C             ISEL      : PARAMETER THAT SELECTS HOW THE GEOMETRY
C                         DEPENDENT PARAMETERS ARE CALCULATED 
C                         1 ASSUME HAMADA COORDINATES AND USE THE
C                           ROUTINE VISGEOM
C                         2 ASSUME CIRCULAR GEOMETRY AND USE CIRC-
C                           GEOM TO DETERMINE THE PARAMETERS
C                         3 READ FROM FILE
C             ISHOT     : SHOT NUMBER, ONLY USED WHEN ISEL = 3
C             NREG      : PARAMETER THAT FORCES A CERTAIN REGIME
C                         IN THE CALCULATION OF THE VISCOSITY.
C                         0     THE NORMAL VALUE. THE BANANA AND
C                               PFIRSCH SCHLUETER CONTRIBUTIONS
C                               ARE WEIGHTED
C                         1     THE BANANA REGIME IS FORCED.
C                         OTHER THE PFIRSCH SCHLUETER REGIME IS
C                               FORCED.
C                         NOTE : FORCING THE BANANA REGIME IS NOT 
C                                THE SAME THING AS CALCULATING THE
C                                BANANA PLATEAU CONTRIBUTION. FOR
C                                HIGH COLLISION FREQ. THIS CONTRI-
C                                BUTION USUALLY DECREASES. WHEN
C                                NREG = 1 THERE IS NO SUCH DECR.
C             SIGMA(4)    THE SIGMA'S DETERMINE WHICH OF THE COUPLING 
C                         TERMS IN THE EQUATION FOR THE PFIRSCH 
C                         SCHLUETER REGIME ARE TAKEN INTO ACOUNT. 
C                         SIGMA(1)= 0,1  IS THE CROSS COUPLING BETWEEN 
C                                     THE HEAT EQUATION AND THE EQUATION
C                                     FOR THE THIRD LAGUERRE HARM.
C                         SIGMA(2)= 0.,1 IS THE CROSS COUPLING BETWEEN 
C                                     THE EQUATION FOR THE THIRD HARMONIC
C                                     AND THE TEMPERATURE PERTURB.
C                         SIGMA(3)= 0,1 IS THE ENERGY EXCHANGE BETWEEN 
C                                     THE DIFFERENT SPECIES
C                         SIGMA(4)= 0,1 IS THE COUPLING TO THE "DENSITY 
C                                     PERTURB." IN THE EQUATION FOR THE 
C                                     THRIRD LAGUERRE HARMONIC
C                         IN THE APPROXIMATIONS OF HIRSHMANN AND SIGMAR
C                         SIGMA(1)=0 SIGMA(2)=0 SIGMA(3)=1 SIGMA(4)=1
C                         THAT IS THE FIRST TWO TERMS ARE NEGLECTED.
C                         THE USUAL VALUES OF ALL THESE COEFFICIENTS IS 1.
C             NLEG      : NUMBER OF LEGENDRE POLYNOMALS IN THE 
C                         EXPANSION (MAXIMUM AND NORMAL VALUE IS 3)
C             NENERGY   : PARAMETER THAT DETERMINES WHETHER ENERGY
C                         SCATTERING IS TAKEN INTO ACCOUNT IN THE 
C                         CALCULATION OF THE VISCOSITY 
C                         0 NOT ACCOUNTED FOR 
C                         1 ACCOUNTED FOR (NORMAL VALUE)
C             NCOF        PARAMETERS THAT DETERMINES WHETHER ION-
C                         ELECTRON COLLISIONS ARE TAKEN INTO ACCOUNT
C                         0 NOT TAKEN INTO ACCOUNT
C                         1 TAKEN INTO ACCOUNT (NORMAL VALUE, 1) 
C                         0 ONLY TO OBTAIN SOME ANALYTIC RESULTS)
C             NEOGEO      LOGICAL IF TRUE THEN THE GEOMETRY DEP.
C                         PARAMETERS ARE RECALCULATED.
C             NEOFRC      LOGICAL IF TRUE THEN THE FRICTION AND 
C                         VISCOSITY MATRIX ARE NOT NEWLY CALCULATED
C             IC        : THE CONTRIBUTION FOR WHICH THE COEFFICIENTS
C                         ARE CALCULATED. 
C                         0 THE CLASSICAL PARTICLE FLUX
C                         1 BANANA PLATEAU CONTRIBUTION
C                         2 PFIRSCH SCHLUETER CONTRIBUTION
C                         3 BOTH BANANA-PLATEAU, CLASSICAL, AND P.S.
C             EPARR     : THE PARALLEL ELECTRIC FIELD TIMES THE 
C                         MAJOR RADIUS, I.E. THE LOOP VOLTAGE. 
C             COEFF     : ARRAY(NSM,NCM,4) THAT CONTAINS THE TRANSP.
C                         COEFFICIENTS OF SPECIES NSM,NCM THE LAST
C                         INDEX GIVES
C                         1 PARTICLE FLUX   \GAMMA_IJ
C                         2 ENERGY FLUX     Q_IJ / T_I
C                         3 PARALLEL FLOW   IN UNITS OF CURRENT
C                                           <J.B><B>/<B**2>
C                         4 POLOIDAL FLOW.  <U.VT>/<B.VT>
C  
C     THE ROUTINE CALLS THE FOLLOWING ROUTINES
C     PERR    :  ERROR HANDLING
C     GEOM    :  CALCULATE THE GEOMETRY DEPENDENT QUANTITIES
C     MENN    :  CALCULATE THE FRICTION COEFFICIENTS
C     PS      :  CALCULATE THE PFIRSCH SCHLUETER CONTRIBUTION
C     BP      :  CALCULATE THE BANANA PLATEAU CONTRIBUTION
C--------------------------------------------------------------------
 
      IMPLICIT NONE

      include 'elem_config.inc'
      
      integer NSM, NCM 	
      parameter(NSM = NELMAX+2)
      parameter(NCM = NIONMAX)

      INTEGER NMAXGR, NKEEP

      PARAMETER(NKEEP = NCM)
      PARAMETER(NMAXGR = 1000)

      INTEGER NZM,NAR,ISEL,NREG,NLEG,NENERGY,NCOF,IC
      REAL DS,EPARR
      LOGICAL NEOGEO, NEOFRC

      INTEGER NS,NC,MMX,IKEEP,IT,I,J,ISK,ISHOT,
     +        K,L
      REAL M,RHO,DEN,T,TAU,COEFF,ZSP,EPARN,R2I,
     +       EPS,BAV,B2AV,BGRADP,FC,FM,XI,RHOK,
     +       BI2A,RBT,DPSIDR,RNQ,UAI,SIGMA, DUM, LA, LAB, MM, NN,
     +       DENI,TA,TB,MA,MB,GCLASS,COEFFC
      DIMENSION M(NSM), ZSP(NSM,NCM), NC(NSM),
     +          DEN(NSM,NCM),T(NSM),TAU(NSM,NSM),
     +          XI(NSM,NCM),GCLASS(NKEEP), 
     +          DS(NSM,NCM,2),DENI(NSM,NCM),COEFFC(NSM,NCM,4),
     +          COEFF(NSM,NCM,4),SIGMA(4),LA(3,3,NSM),
     +          BAV(NKEEP),B2AV(NKEEP),BGRADP(NKEEP),FC(NKEEP),
     +          RHOK(NKEEP),FM(NKEEP,NMAXGR),MMX(NKEEP),R2I(NKEEP),
     +          ISK(NKEEP),BI2A(NKEEP),RBT(NKEEP),DPSIDR(NKEEP),
     +          RNQ(NKEEP),UAI(NSM,NCM,3),LAB(3,3,NSM,NSM),
     +          MM(3,3),NN(3,3)

      SAVE XI, LA, LAB, DENI, TAU

C     STORE THE GEOMETRY DEPENDENT VARIABLES
      COMMON / KGEOM / BAV, B2AV, BGRADP, FC, FM, MMX, RHOK,
     +         BI2A, RBT, DPSIDR, RNQ, ISK

      DATA IKEEP / 0  / 

      IF  ((NEOFRC).AND.(NEOGEO)) CALL PERR(16)

C     IF THE FRICTION AND VISCOSITY COEFFICIENTS ARE NEWLY CAL-
C     CULATED THEN GO THROUGH THE FOLLOWING LOOP. 
      IF (.NOT.NEOFRC) THEN

C     STORE THE DENSITY VALUES
      DO 111 J = 1, NCM
         do 111 I=1,NSM
        DENI(I,J) = DEN(I,J)
 111  CONTINUE

      IF ((NAR.NE.NSM).OR.(NZM.NE.NCM)) CALL PERR(1)
      IF (RHO.EQ.0.) CALL PERR(2)

      IF ((IC.EQ.1).OR.(IC.EQ.3)) THEN
        IF (NREG.NE.0) CALL PERR(11)
        IF (NLEG.NE.3) CALL PERR(12)
        IF (NENERGY.NE.1) CALL PERR(13)
      ENDIF
      IF ((IC.EQ.2).OR.(IC.EQ.3)) THEN
        DUM = 1.
        DO 1 I = 1, 4
          DUM = DUM*SIGMA(I)
 1      CONTINUE
        IF (DUM.NE.1.) CALL PERR(15)
      ENDIF
      IF (NCOF.NE.1) CALL PERR(14)


C     IS THE LOGICAL NEOGEO SET ????
      IF (NEOGEO) THEN
        IKEEP = 1
        IT = IKEEP
        RHOK(IT) = RHO
        ISK(IT) = ISEL
        CALL GEOM(NMAXGR,ISEL,ISHOT,RHO,EPS,BAV(IT),B2AV(IT),BI2A(IT),
     +    RBT(IT),BGRADP(IT),DPSIDR(IT),RNQ(IT),FC(IT),GCLASS(IT),
     +    FM(IT,1),MMX(IT),R2I(IT))
      ELSE
C       DOES THE VALUE OF RHO EXIST
        IT = 0
        DO 10 I = 1, IKEEP
          IF ((RHO.EQ.RHOK(I)).AND.(ISEL.EQ.ISK(I)))
     +      IT = I
 10     CONTINUE
        IF (IT.EQ.0) THEN
          IKEEP = IKEEP + 1
          IF (IKEEP.GT.NKEEP) THEN
            IKEEP = NKEEP
            CALL PERR(17)
          ENDIF
          IT = IKEEP
          RHOK(IT) = RHO
          ISK(IT) = ISEL
          CALL GEOM(NMAXGR,ISEL,ISHOT,RHO,EPS,BAV(IT),B2AV(IT),
     +              BI2A(IT),RBT(IT),BGRADP(IT),DPSIDR(IT),RNQ(IT),
     +              FC(IT),GCLASS(IT),FM(IT,1),MMX(IT),R2I(IT))
        ENDIF
      ENDIF

      CALL COLXI(NSM,NCM,NS,NC,ZSP,DEN,T,M,TAU,XI)

C     SWITCH OF ION - ELECTRON COLLISIONS? 
      IF (NCOF.EQ.0) THEN
C       YES SWITCH OF
        J = 0
        DO 100 I = 1, NS
          IF (ZSP(I,1).LT.0.) J = I
 100    CONTINUE
        IF (J.EQ.0) GOTO 200
        DO 300 I = 1, NS
          IF (J.NE.I) TAU(I,J) = 0.
 300    CONTINUE
      ENDIF
 200  CONTINUE

C     CALCULATED THE FRICTION COEFFICIENTS
      DO 20 I = 1, NS
        DO 30 K = 1, 3
          DO 30 L = 1, 3
            LA(K,L,I) = 0.
 30     CONTINUE
        DO 40 J = 1, NS
          TA = T(I)
          TB = T(J)
          MA = M(I)
          MB = M(J)
          CALL MENN(TA,TB,MA,MB,MM,NN)
          DO 50 K = 1, 3
            DO 50 L = 1, 3
               LAB(K,L,I,J) = TAU(I,J)*NN(K,L)
               LA(K,L,I) = LA(K,L,I) + TAU(I,J)*MM(K,L)
 50       CONTINUE
 40     CONTINUE
 20   CONTINUE

C     END OF THE IF STATEMENT ON NEOFRC
      ELSE
C       RESTORE THE DENSITY VALUES
         DO 222 J = 1, NCM
            do 222 I=1,NSM
               DEN(I,J) = DENI(I,J)
 222        CONTINUE
      ENDIF      


C     MAKE THE COEFFICIENTS ZERO
      DO 4000 I = 1, NSM
        DO 4000 J = 1, NCM
          DO 4000 K = 1,4
            COEFF(I,J,K) = 0.
            COEFFC(I,J,K) = 0.
 4000 CONTINUE

C     NORMALIZE THE EPARR
      EPARN = 1E-3*DPSIDR(IT)*R2I(IT)*EPARR 

C     CALCULATE THE CLASSICAL TRANSPORT CONTRIBUTION
      IF ((IC.EQ.0).OR.(IC.EQ.3)) THEN
        CALL CLASS(NSM,NCM,NS,NC,T,ZSP,XI,LA,LAB,DS,COEFFC)
C       ADD THE CORRECT NORMALIZATION FACTOR, PUT RESULT 
C       IN COEFF AND CLEAR COEFFC
        DO 3010 I = 1, NS
          DO 3010 J = 1, NC(I)
            DO 3010 K = 1, 2
              COEFF(I,J,K) = COEFF(I,J,K)+ COEFFC(I,J,K)*
     +          GCLASS(IT)/(1.6022E-22*ZSP(I,J)*DPSIDR(IT)**2)
              COEFFC(I,J,K) = 0.
 3010   CONTINUE
      ENDIF

C     CALCULATE THE PFIRSCH SCHLUETER CONTRIBUTION
      IF ((IC.EQ.2).OR.(IC.EQ.3)) THEN
        call PS(NSM,NCM,ns,nc,NEOFRC,tau,XI,LA,LAB,m,t,DENI,
     +           ZSP,sigma,ds,RNQ(it),uai,coeffc)
C       ADD THE CORRECT NORMALIZATION FACTOR, PUT RESULT IN 
C       COEFF AND CLEAR COEFFC
        do 3000 i = 1, ns
          do 3000 j = 1, nc(i)
            do 3100 k = 1, 2
              coeff(i,j,k) = coeff(i,j,k) - coeffc(i,j,k)*
     +          (1.- B2AV(it)*BI2A(it))*(RBT(it)/DPSIDR(it))**2/
     +          (1.60E-22*ZSP(I,J)*b2av(it))
              coeffc(i,j,k) = 0.
 3100       CONTINUE
 3000   CONTINUE
      ENDIF


      IF ((IC.EQ.1).OR.(IC.EQ.3)) THEN
        CALL BP(NSM,NCM,NS,NC,NEOFRC,T,M,DEN,ZSP,TAU,XI,LA,
     +   LAB,EPS,NENERGY,NLEG,NREG,BGRADP,FC,FM,MMX,DS,EPARN,
     +   COEFFC)
C       ADD THE CORRECT NORMALIZING FACTOR. PUT RESULT IN 
C       COEFF AND CLEAR COEFFC
        do 3200 i = 1, ns
          do 3200 j = 1, nc(i)
            do 3300 k = 1, 2
              coeff(i,j,k) = coeff(i,j,k) - coeffc(i,j,k)*
     +          (RBT(it)/DPSIDR(it))**2/(1.60E-22*ZSP(I,J)*
     +           b2av(it))
              coeffc(i,j,k) = 0.
 3300       CONTINUE
 3200   CONTINUE
      ENDIF 

      do 3400 i = 1, ns
        do 3400 j = 1, nc(i)
          coeff(i,j,3) = coeffc(i,j,3)*1600*den(i,j)*zsp(i,j)*
     +                   rbt(it)*bav(it)/(dpsidr(it)*b2av(it))
          coeff(i,j,4) = 1e3*coeffc(i,j,4)*rbt(it)/(b2av(it)*
     +                   dpsidr(it))
 3400 continue
      
      RETURN
      END

