
      SUBROUTINE VISCOS(NSM, NCM, NS, IS, IC, NREG, BGRADP, FC, FM, 
     +              MMX, NENERGY, EPS, TAU, M, T, XI, DEN, NCC, MU)
C--------------------------------------------------------------------
C     THIS ROUTINE CALCULATES THE VISCOSITY COEFFICIENTS USING THE 
C     EXPRESSIONS OF K.C. SHAING, PHYS. PLASMAS 3 965 (1996)
C     EXTENDED TO THE 3 LAGUERRE POLYNOMIAL EXPANSION
C
C     INPUT   NSM            : MAXIMUM NUMBER OF SPECIES. THIS PAR-
C                              AMETER IS GIVEN TO THE ROUTINE TO 
C                              CHECK INTERNAL CONSISTENCY
C             NCM            : MAXIMUM NUMBER OF CHARGE STATES, USED
C                              FOR CHECK ON INTERNAL CONSISTENCY
C             NS             : THE ACTUAL NUMBER OF SPECIES
C             IS             : THE SPECIES NUMBER FOR WHICH THE 
C                              VISCOSITY IS TO BE CALCULATED
C             IC             : THE CHARGE STATE FOR WHICH THE VISC-
C                              OSITY IS TO BE CALCULATED.
C             NREG           : INTEGER THAT DETERMINES WHAT REGIME
C                              IS TO BE FORCED IN THE CALCULATION 
C                              OF THE VISCOSITY
C                              0 = WEIGTED OVER ALL REGIMES
C                              1 = THE BANANA REGIMES IS FORCED
C                              OTHER = THE PFIRSCH SCHLUETER REGIME
C                              IS FORCED. 
c             BGRADP         : THE FLUX SURFACE QUANTITY N IN GRAD T
C                              (THE INNER PRODUCT OF THE UNIT VECTOR
C                              IN THE DIRECTION OF THE MAGNETIC FIELD
C                              AND THE GRADIENT OF THE POLOIDAL 
C                              ANGLE)
C             FC             : THE NUMBER OF PASSING PARTICLES
C             FM             : THE COEFFICIENTS NECESSARY TO CALCULATE
C                              THE VISCOSITY IN THE PFIRSCH SCHLUETER
C                              REGIME.
C             MMX            : THE NUMBER OF COEFFICIENTS FM
C             EPS            : THE ACCURACY WITH WHICH THE NUMERICAL
C                            : COEFICIENTS SHOULD BE CALCULATED. 
C                              NOTE THAT THE ACCURACY OF THE ANALYTIC
C                              EQUATIONS IMPLEMENTED HERE IS NOT 
C                              EXPECTED TO BE LARGER THAN 10% IN 
C                              CERTAIN COLLISIONAL REGIMES. WHEN EPS = 0., 
C                              THE VALUE 0.001 WILL BE ASSUMED.
C             TAU            : ARRAY(NSM,NSM) THE COLLISION FREQUENCY
C                              N_I M_I / TAU_IJ
C             M              : ARRAY(NSM) THE MASSES OF THE SPECIES 
C                              IN KG
C             T              : ARRAY(NSM) TEMPERATURE OF THE SPECIES 
C                              IN KEV
C             XI             : ARRAY(NSM,NCM) RELATIVE WEIGTH OF A 
C                              CHARGE STATE.
C             DEN            : ARRAY(NSM,NCM) DENSITY IN UNIT 10^19 M^-3
C             NCC            : ARRAY(NSM) THAT GIVES THE NUMBER OF 
C                              CHARGE STATES PER SPECIES.
C     OUTPUT  MU             : ARRAY(3,3) THE VISCOSITY COEFFICIENTS
C
C     THE ROUTINE CALLS THE FOLLOWING SUBROUTINE
C     VISCOL    : ROUTINE THAT CALCULATES THE ENERGY DEPENDENT 
C                 COLLISION AND ENERGY SCATTERING FREQUENCY
C     PERR      : ROUTINE THAT DOES THE ERROR HANDLING.
C--------------------------------------------------------------------
       
      IMPLICIT NONE

      REAL  EPS, MU, TAU, T, M, XI, DEN
      INTEGER NENERGY, NSM, NCM, NCC, IS, IC, NS, NREG
      INTEGER MMX,I,J,NV,K
      LOGICAL NLTEST, NLERR
      DIMENSION MU(3,3), TAU(NSM,NSM), T(NSM), M(NSM), XI(NSM,NCM),
     +          NCC(NSM), DEN(NSM,NCM)
      REAL BGRADP, TWOPI, FC, FT, PI, X, DUM, VTH,
     +       NUD, NUE, NUT, KB, KPS, OMMN, NTOM, KTOT, VOORF,
     +       DX, WEIGHT, MEANMUO
      REAL FM(MMX), A(6), MUO(3,3)

      NLTEST = .FALSE.

C     CALCULATE VTH
      VTH = SQRT(2.*1.6E-16*T(IS)/M(IS))

C     THE CONSTANT PI, 2 PI
      PI = 4.*ATAN(1.)
      TWOPI = 2.*PI

C     THE CONSTANTS TO APPROXIMATE THE NUT*IRM FUNCTION
      A(1) = 2. / 5.
      A(2) = - 22. / 105.
      A(3) = 6. / 35.
      A(4) = - 34. / 231.
      A(5) = 166. / 1287.
      A(6) = - 82. / 715.


      FT = 1.- FC

      NV = 100
      DO 2500 K = 1, 9
        MUO(K,1) = 0.
 2500 CONTINUE
 3000 CONTINUE

C       THE LOOP OVER THE VELOCITY
        DX = 10./REAL(NV)
        DO 3050 K = 1, 9
          MU(K,1) = 0.
 3050   CONTINUE

        DO 3100 I = 0, NV 

          X = I*DX
          OMMN = VTH * X * BGRADP
        
          IF (I.EQ.0) THEN

            VOORF = 0.
 
          ELSE 

c         NOW CALCULATE THE COLLISION AND ENERGY SCATTERING 
C         FREQUENCY. 
          CALL VISCOL(NSM, NCM, NS, IS, IC, TAU, M, T, XI, NCC, X,
     +                DEN, NUD, NUE)
          NUT = 3*NUD + REAL(NENERGY)*NUE

          KB  = FT * NUD / FC
          KPS = 0.
          DO 3200 J = 1, MMX
            NTOM = NUT / ( OMMN * REAL(J) )
            IF (NTOM.LT.10) THEN
              DUM  = -1.5*NTOM**2 - 4.5*NTOM**4 + 
     +             (0.25+(1.5+2.25*NTOM**2)*NTOM**2)*2*NTOM *
     +             ATAN(1/NTOM)
            ELSE
              NTOM = (1. / NTOM)**2 
              DUM = A(1)
              DO 3300 K = 2, 6
                DUM = DUM + A(K)*NTOM**(K-1)
 3300         CONTINUE
            ENDIF
            KPS = KPS + FM(J)*DUM
 3200     CONTINUE
          KPS = KPS * 1.5 * (VTH * X)**2/NUT
          IF (NREG.EQ.0) THEN
            KTOT = KB*KPS/(KB+KPS)
          ELSE
            IF (NREG.EQ.1) THEN
              KTOT = KB
            ELSE 
              KTOT = KPS
            ENDIF
          ENDIF
          IF ((NLTEST).AND.(NV.EQ.100)) THEN
            WRITE(*,*)NUD/(X*VTH), KTOT/(X*VTH)
          ENDIF

          VOORF = X**4 * EXP(-X**2)*KTOT*DX/3.

          ENDIF

          if ((i.eq.0).or.(i.eq.nv)) then
             weight=1.
           else
             weight=2.*2.**mod(i,2)
          endif

          MU(1,1) = MU(1,1) + weight*VOORF
          MU(1,2) = MU(1,2) + weight*VOORF*(X**2-2.5)
          MU(2,2) = MU(2,2) + weight*VOORF*(X**2-2.5)**2
          MU(1,3) = MU(1,3) + weight*VOORF*(35./8.-3.5*X**2 + 
     +              0.5*X**4)
          MU(2,3) = MU(2,3) + weight*VOORF*(35./8.-3.5*X**2 + 
     +              0.5*X**4)*(X**2-2.5)
          MU(3,3) = MU(3,3) + weight*VOORF*(35./8.-3.5*X**2 + 
     +              0.5*X**4)**2
                      
 3100   CONTINUE

        VOORF = 1.E19*DEN(IS,IC)*M(IS)* 8. / (3.*SQRT(PI))
        MU(2,1) = MU(1,2)
        MU(3,1) = MU(1,3)
        MU(3,2) = MU(2,3)
        MEANMUO = 0.
        DO 3150 K = 1, 9
          MU(K,1) = MU(K,1) * VOORF
          MEANMUO = MEANMUO + MU(K,1)
 3150   CONTINUE
        MEANMUO = MEANMUO / 9.
        NLERR = .FALSE.
        DO 3160 K = 1, 9
          IF (ABS(MU(K,1)-MUO(K,1)).GT.ABS(EPS*MU(K,1))) THEN 
            IF (ABS(MU(K,1)).GT.EPS*MEANMUO) NLERR = .TRUE.
          ENDIF
 3160   CONTINUE
        IF (NLERR) THEN
          DO 3170 K = 1, 9
            MUO(K,1) = MU(K,1)
 3170     CONTINUE
          NV = NV*2
          IF (NV.GT.4e6) CALL PERR(4)
          GOTO 3000
        ENDIF
         

      RETURN 
      END
