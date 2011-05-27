      SUBROUTINE BP(NSC,NCC,NS,NC,NEOFRC,T,M,DEN,ZSP,TAU,XI,LA,
     +     LAB,EPS,NENERGY,NLEG,NREG,BGRADP,FC,FM,MMX,DSI,EPAR,
     +     COEFF)
C--------------------------------------------------------------------
C     INPUT  NSC NCC : THE PARAMETERS NSM AND NCM USED TO CHECK
C                      CONSISTENCY
C            NS      : THE NUMBER OF SPECIES
C            NC      : ARRAY(NSM) THE NUMBER OF CHARGE STATES
C                      OF EVERY SPECIES.
C            NEOFRC  : LOGICAL IF TRUE THE MATRICES ARE NOT 
C                      NEWLY CALCULATED. 
C            T       : ARRAY(NSM) THE TEMPERATURE OF THE SPECIES
C                      IN KEV
C            M       : ARRAY(NSM) THE MASS OF THE SPECIES IN KG.
C            DEN     : ARRAY(NSM,NCM) DENSITY OF EVERY COMPONENT
C                      IN UNITS 10^19 M^-3
C            ZSP     : ARRAY(NSC,NCC) THE CHARGE OF EVERY COMP.
C                      NORMALIZED TO E 
C            TAU     : ARRAY(NSM,NSM) THE COLLISION FREQ. 
C                      M_I N_I / TAU_IJ
C            XI      : ARRAY(NSC,NCC) THE RELATIVE WEIGTH OF 
C                      EVERY CHARGE STATE. 
C            LA      : THE FRICTION COEFFICIENTS
C            LAB     : THE FRICTION COEFFICIENTS
C            EPS     : DESIRED ACCURACY
C            NENERGY : INTEGER IF 1 ENERGY SCATTERING IS USED
C                      IN THE CALCULATION OF THE VISCOSITY
C            NLEG    : THE NUMBER OF LAGUERRE POLYNOMALS 
C            NREG    : A PARAMETER THAT FORCES A CERTAIN REGIME
C                      IN THE CALCULATION OF THE VISCOSITY
C                      0 WHEIGTHED REGIME VALID FOR ALL COLLISION
C                        FREQUENCIES
C                      1 FORCE THE BANANA REGIME
C                      OTHER FORCE THE PFIRSCH SCHLUETER REGIME
C            BGRADP  : THE FLUX SURFACE QUANTITY N IN GRAD THETA
C                      THE INNER PRODUCT OF THE UNIT VECTOR IN 
C                      THE DIRECTION OF THE MAGNETIC FIELD AND 
C                      THE GRADIENT OF THE POLOIDAL ANGLE
C            FC      : THE NUMBER OF PASSING PARTICLES
C            FM      : THE COEFFICIENTS NECESSARY TO CALCULATE
C                      THE VISCOSITY IN THE PFIRSCH SCHLUETER 
C                      REGIME. 
C            MMX     : THE NUMBER OF FM COEFFICIENTS
C            DSI     : THE THERMODYNAMIC FORCE, FOR A DEFINITION
C                      SEE ROUTINE NEOART
C     OUTPUT COEFF   : NORMALIZED TRANSPORT COEFFICIENTS
C
C     THE ROUTINE CALLS THE FOLLOWING SUBROUTINES
C     PERR    : FOR ERROR HANDLING
C     VISCOS  : TO CALCULATE THE VISCOSITY
C     LUDCMP  : TO DO THE LU DECOMPOSITION
C     LUBKSB  : TO SOLVE THE MATRIX EQUATION
C--------------------------------------------------------------------
      IMPLICIT NONE

      include 'elem_config.inc'
      
C     MAXIMUM NUMBER OF SPECIES AND CHARGE STATES
      
      integer NSM, NCM 	
      parameter(NSM = NELMAX+2)
      parameter(NCM = NIONMAX)

      INTEGER  NS, NC, I, J, MT, NREG, 
     +         INDX, NSC, NCC, NLEG, MMX, K, L
      REAL T, M, DEN, TAU, AAA, D, AKP, LA, LAB,
     +       VISC, XI, FC, FM, BGRADP,ZSP, EPAR,
     +       AINV, DS, COEFF, BT, UPAR,DSI
      DIMENSION NC(NSM),T(NSM),M(NSM),DEN(NSM,NCM),UPAR(NSM,NCM,3),
     +          TAU(NSM,NSM), VISC(3,3,NSM,NCM), FM(MMX), XI(NSM,NCM),
     +          AAA(3,3,NSM,NSM), AKP(3,3),ZSP(NSM,NCM),
     +          D(3*NSM),INDX(3*NSM), LA(3,3,NSM),LAB(3,3,NSM,NSM),
     +          AINV(3*NSM,3*NSM),DS(NSM,NCM,2),COEFF(NSM,NCM,4),
     +          BT(3*NSM),DSI(NSM,NCM,2)
      LOGICAL NEOFRC


      SAVE VISC,AAA,AINV,INDX

      INTEGER NENERGY
      REAL EPS, MU
      DIMENSION MU(3,3)
      LOGICAL NLTEST

C     CHECK THE CONSISTENCY OF THE NSM NCM PARAMETERS
      IF ((NSC.NE.NSM).OR.(NCC.NE.NCM)) CALL PERR(1)

C     DEFINE THE PROPER THERMODYNAMIC FORCES FOR THE CALCULATION BELOW
      DO 111 I = 1, NS
        DO 111 J = 1, NC(I)
          DO 111 K = 1, 2
            DS(I,J,K) = T(I)*DSI(I,J,K)/ZSP(I,J)
 111  CONTINUE

C     CHECK WHETHER THE MATRICES HAVE TO BE RECALCULATED
      IF (.NOT.NEOFRC) THEN 

C     CALCULATE THE VISCOSITY
      DO 600 I = 1, NS
        DO 600 J = 1, NC(I)
          DO 610 K = 1, 9
            VISC(K,1,I,J) = 0.
 610      CONTINUE
          NLTEST = .FALSE.
          CALL VISCOS(NSM, NCM, NS, I, J, NREG, BGRADP, FC, FM, MMX, 
     +                NENERGY, EPS, TAU, M, T, XI, DEN, NC, MU)
          DO 700 K = 1, NLEG
            DO 700 L = 1, NLEG 
              VISC(K,L,I,J) = MU(K,L)
 700      CONTINUE
 600  CONTINUE

c     CALCULATE THE A MATRICES
      DO 1000 I = 1, NS 
        DO 1000 J = 1, NC(I)
          DO 1200 K = 1, 9
            AAA(K,1,I,J) = 0.
 1200     CONTINUE
          DO 1300 K = 1, 3
            AAA(K,K,I,J) = 1.
 1300     CONTINUE
          DO 1100 K = 1, NLEG
            DO 1100 L = 1, NLEG
              AAA(K,L,I,J) = VISC(K,L,I,J) - XI(I,J)*LA(K,L,I)
 1100     CONTINUE
C         DO THE LU DECOMPOSITION
          CALL LUDCMP(AAA(1,1,I,J),3,3,INDX,D)
          DO 1400 K = 1, 9
            AKP(K,1) = 0.
 1400     CONTINUE
          DO 1500 K = 1, 3
            AKP(K,K) = 1.
 1500     CONTINUE
          DO 1600 K = 1, 3
            CALL LUBKSB(AAA(1,1,I,J),3,3,INDX,AKP(1,K))
 1600     CONTINUE
          DO 1700 K = 1, 9
            AAA(K,1,I,J) = AKP(K,1)
 1700     CONTINUE

 1000 CONTINUE
 
C     NOW CALCULATE THE MATRIX FOR INVERSION
      DO 2000 I = 1, 3*NS
        DO 2000 J = 1, 3*NS
          AINV(I,J) = 0.
 2000 CONTINUE
      DO 2100 I = 1, 3*NS
        AINV(I,I) = 1.
 2100 CONTINUE
      DO 2200 I = 1, NS
        DO 2300 K = 1, 9
          AKP(K,1) = 0.
 2300   CONTINUE
        DO 2400 K = 1, 9
          DO 2400 J = 1, NC(I)
            AKP(K,1) = AKP(K,1) + XI(I,J)**2*AAA(K,1,I,J)
 2400   CONTINUE
        DO 2500 J = 1, NS
          DO 2500 K = 1, NLEG
            DO 2500 L = 1, NLEG
              DO 2500 MT = 1, 3
                AINV((I-1)*3+K,(J-1)*3+L) = AINV((I-1)*3+K,(J-1)*3+L)
     +             - AKP(K,MT)*LAB(MT,L,I,J)
 2500   CONTINUE
 2200 CONTINUE

C     CALCULATE THE LU DECOMPOSITION
      CALL LUDCMP(AINV,3*NS,3*NSM,INDX,D)  

C     END THE IF STATEMENT THAT DETERMINES WHETHER THE MATRICES ARE 
C     RECALCULATED
      ENDIF

C     NOW CALCULATE THE RIGHT HAND SIDE    
      DO 3000 I = 1, NS
        DO 3000 K = 1, 3 
          BT((I-1)*3+K) = 0.
          IF (K.GT.NLEG) GOTO 3000
          DO 3050 J = 1, NC(I)
            BT((I-1)*3+K) = BT((I-1)*3+K)+1.6*XI(I,J)*ZSP(I,J)*
     +                      DEN(I,J)*AAA(K,1,I,J)*EPAR
 3050     CONTINUE
          DO 3100 L = 1, 3
            DO 3100 MT = 1, 2
              DO 3100 J = 1, NC(I)
                BT((I-1)*3+K) = BT((I-1)*3+K) - XI(I,J)*
     +            AAA(K,L,I,J)*VISC(L,MT,I,J)*DS(I,J,MT)
 3100     CONTINUE
 3000 CONTINUE


      CALL LUBKSB(AINV,3*NS,3*NSM,INDX,BT)

C     NOW CALCULATE THE INDIVUDUAL PARTICLE VELOCITIES
      DO 4000 I = 1, NS
        DO 4100 K = 1, 3 
          AKP(K,1) = 0.
          IF (K.GT.NLEG) GOTO 4100
          DO 4200 L = 1, 3
            DO 4200 J = 1, NS
              AKP(K,1) = AKP(K,1)+LAB(K,L,I,J)*BT((J-1)*3+L)
 4200     CONTINUE
 4100   CONTINUE
        DO 4300 J = 1, NC(I)
          DO 4400 K = 1, 3
            UPAR(I,J,K) = 0.
            IF (K.GT.NLEG) GOTO 4400
            DO 4410 L = 1, 3
              AKP(L,2) = 0.
              DO 4500 MT = 1, 2
                AKP(L,2) = AKP(L,2)+VISC(L,MT,I,J)*DS(I,J,MT)
 4500         CONTINUE
              UPAR(I,J,K) = UPAR(I,J,K) + AAA(K,L,I,J)*(
     +          XI(I,J)*AKP(L,1)-AKP(L,2))
 4410       CONTINUE
            UPAR(I,J,K) = UPAR(I,J,K) + AAA(K,1,I,J)*
     +          1.6*ZSP(I,J)*DEN(I,J)*EPAR
 4400     CONTINUE
 4300   CONTINUE
 4000 CONTINUE


C     FINALLY CALCULATE THE COEFFICIENTS COEFF
      DO 5000 I = 1, NS
        DO 5000 K = 1, NC(I)

          COEFF(I,K,3) = UPAR(I,K,1)
          DO 5100 MT = 1, 3
            IF (MT.NE.3) 
     +        UPAR(I,K,MT) = UPAR(I,K,MT) + DS(I,K,MT)
            COEFF(I,K,1) = COEFF(I,K,1) + VISC(1,MT,I,K)*
     +                     UPAR(I,K,MT)
            COEFF(I,K,2) = COEFF(I,K,2) + VISC(2,MT,I,K)*
     +                     UPAR(I,K,MT)
5100      CONTINUE
          COEFF(I,K,4) = UPAR(I,K,1)

 5000 CONTINUE
 
      RETURN
      END
