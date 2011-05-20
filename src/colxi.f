
      SUBROUTINE COLXI(NSM,NCM,NS,NC,ZSP,DEN,T,M,TAU,XI)
C--------------------------------------------------------------------
C     THIS SUBROUTINE CALCULATES THE COLLISION FREQUENCIES AND 
C     THE WEIGHTING FACTORS XI
C
C     INPUT    NSM    : MAXIMUM NUMBER OF SPECIES
C              NCM    : MAXIMUM NUMBER OF CHARGE STATES
C              NS     : ACTUAL NUMBER OF SPECIES
C              NC     : ARRAY(NSM) GIVES THE ACTUAL AMOUNT OF 
C                       CHARGE STATES PER SPECIES
C              ZSP    : ARRAY(NSM,NCM) THE CARGE NUMBER
C              DEN    : ARRAY(NSM,NCM) THE DENSITY IN 10^19 M^-3
C              T      : ARRAY(NSM) TEMPERATURE OF SPECIES
C              M      : ARRAY(NSM) MASS OF SPECIES
C     OUTPUT   TAU    : ARRAY(NSM,NSM) COLLISION FREQUENCY 
C                       WEIGHTED OVER CHARGE STATES.
C              XI     : RELATIVE WEIGHT OF EVERY CHARGE STATE.
C--------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER NSM,NCM,I,J,NS,NC,K,L
      REAL*8 DEN,T,M,DENE,LNAB,TAU,XI,TCONS,ZSP
      DIMENSION ZSP(NSM,NCM),NC(NSM),DEN(NSM,NCM),T(NSM),
     +          M(NSM),XI(NSM,NCM),TAU(NSM,NSM)
 
C     CALCULATE THE COLLISION TIMES
      TCONS = 8.74202E6
C     THE COLLISION TIME BELOW IS THE DOUBLE AVERAGE (OVER CHARGE 
C     STATES) OF THE DENSITY DEVIDED BY THE COLLISION TIME (EQUATION 
C     A3 AND A4 OF HOULBERG, PHYS PLASMAS 4 3230 (1997)) 
C     FIRST DETERMINE THE ELECTRON DENSITY
      DENE = 0.
      DO 1 I = 1, NS
        IF (ZSP(I,1).LT.0.) DENE = DEN(I,1)
 1    CONTINUE
      IF (DENE.EQ.0.) THEN
      DO 2 I = 1, NS
        DO 2 J = 1, NC(I)
          DENE = DENE + DEN(I,J)*ZSP(I,J)
 2    CONTINUE
      ENDIF
   
      DO 3 I = 1, NS
        DO 3 J = 1, NS

C         CALCULATE THE COULOMB LOGARITHM 
          IF ((ZSP(I,1).LT.0.).AND.(ZSP(J,1).LT.0.)) THEN
C           ELECTRON-ELECTRON
            LNAB = 14.9 - 0.5*LOG(0.1*DENE)+LOG(T(I))
          ELSE
            IF (ZSP(I,1).LT.0.) THEN
C             ELECTRON-ION
              LNAB = 15.2-0.5*LOG(0.1*DENE) + LOG(T(I))
            ELSE
              IF (ZSP(J,1).LT.0.) THEN
C               ION-ELECTRON
                LNAB = 15.2-0.5*LOG(0.1*DENE) + LOG(T(J))
              ELSE
C               ION-ION
                LNAB = 17.3-0.5*LOG(0.1*DENE)+
     +                 1.5*LOG(0.5*(T(I)+T(J)))
              ENDIF
            ENDIF
          ENDIF

          TAU(I,J) = 0.
          DO 4 K = 1, NC(I)
            DO 4 L = 1, NC(J)
              TAU(I,J) = TAU(I,J) + DEN(I,K)*DEN(J,L)*ZSP(I,K)**2
     +                   *ZSP(J,L)**2
 4        CONTINUE

c         WRITE(*,*)I,J,NC(I),NC(J),DEN(I,1),ZSP(I,1),LNAB
          TAU(I,J) = TAU(I,J)*TCONS*LNAB*SQRT(M(I)/T(I)**3)
 3    CONTINUE
   
C     CALCULATE THE XI COEFFICIENTS
      DO 50 I = 1, NS
        XI(I,NC(I)) = 0.
        DO 51 J = 1, NC(I)
          XI(I,NC(I)) = XI(I,NC(I)) + DEN(I,J)*ZSP(I,J)**2
 51     CONTINUE
        DO 52 J = 1, NC(I)
          XI(I,J) = DEN(I,J)*ZSP(I,J)**2 / XI(I,NC(I))
 52     CONTINUE
 50   CONTINUE

      RETURN
      END
