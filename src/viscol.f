
     
      SUBROUTINE VISCOL(NSM, NCM, NS, IS, IC, TAU, M, T, XI, NC, X,
     +                  DEN, NUD, NUE)
C--------------------------------------------------------------------
C     THIS SUBROUTINE CALCULATES THE ENERGY DEPENDENT COLLISION
C     FREQUENCY AND THE ENERGY DEPENDT ENERGY SCATTERING.
C
C     INPUT  : NSM    THE MAXIMUM NUMBER OF SPECIES
C              NCM    THE MAXIMUM NUMBER OF CHARGE STATES
C              NS     THE ACTUAL NUMBER OF SPECIES
C              IS     THE SPECIES NUMBER FOR WHICH THE COLLISION 
C                     FREQUENCY IS TO BE CALCULATED
C              IC     THE CHARGE STATE OF THE SPECIES FOR WHICH 
C                     THE COLLISION FREQUENCY IS TO BE CALCULATED.
C              TAU    NORMALIZED COLLISION FREQUENCY (N_I M_I /
C                     TAU_IJ)  ARRAY(NSM,NSM)
C              M      ARRAY(NSM) THE MASS OF THE SPECIES IN KG
C              T      ARRAY(NSM) THE TEMPERATURE OF THE SPECIES 
C                     IN KEV
C              XI     ARRAY(NSM,NCM) THE RELATIVE WEIGHT OF EVERY
C                     CHARGE STATE
C              NC     ARRAY(NSM) THE NUMBER OF CHARGE STATES FOR 
C                     EVERY SPECIES
C              X      NORMALIZED (TO THERMAL) VELOCITY FOR WHICH
C                     THE ENERGY DEPENDENT COLLISION FREQUENCY IS
C                     TO BE CALCULATED
C              DEN    ARRAY(NSM,NCM) DENSITY OF EVERY COMPONENT 
C                     IN 10^19 M^-3
C     OUTPUT   NUD    THE PITCH ANGLE SCATTERING FREQUENCY
C              NUE    THE ENERGY SCATTERING FREQUENCY
C
C     THE ROUTINE CALLS THE FOLLOWING FUNCTION
C     ERF   : CALCULATES THE ERROR FUNCTION 
C
C--------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER NSM, NCM, IS, IC, NC, I, NS 
      REAL TAU,M,T,XI,X,NUD,NUE,DEN,ERF
      REAL PH, G, PI, XAB
      DIMENSION TAU(NSM,NSM), M(NSM), T(NSM), XI(NSM,NCM), NC(NSM),
     +          DEN(NSM,NCM)

      PI = 4.*ATAN(1.)

      NUD = 0.
      NUE = 0.

C     THE LOOP OVER SPECIES
      DO 100 I = 1, NS

C       CALCULATE XAB = VTHB / VTHA
        XAB = SQRT(M(IS)*T(I)/(M(I)*T(IS)))

        PH = ERF(X/XAB)
        G  = (PH - 2*X*EXP(-(X/XAB)**2)/(XAB*SQRT(PI)))
     +       /(2*(X/XAB)**2)

        NUD = NUD + TAU(IS,I) * ( PH - G )/ X**3
        NUE = NUE + TAU(IS,I) * (-2.*PH/X**3 + 4.*(T(IS)/T(I) + 
     +        1./XAB**2) *G/X)

 100  CONTINUE

      NUD = NUD * 0.75 * SQRT(PI) * XI(IS,IC) / (1.E19*DEN(IS,IC)*
     +      M(IS))
      NUE = NUE * 0.75 * SQRT(PI) * XI(IS,IC) / (1.E19*DEN(IS,IC)*
     +      M(IS))

      RETURN
      END
