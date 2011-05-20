
      SUBROUTINE VISFUS(NSM, NCM, NS, IS, IC, NREG, BGRADP, FC, FM, 
     +                  MMX, EPS, TAU, M, T, XI, DEN, NCC, MU)
C     ROUTINE FOR THE SIMPLER AND LESS ACCURATE VALUES OF THE 
C     VISCOSITY COEFFICIENTS. THIS ROUTINE IS NOT READY, IT STILL
C     NEEDS TO BE CHANGED. 
    
      IMPLICIT NONE

      REAL*8  EPS, MU, TAU, T, M, XI, DEN
      INTEGER NENERGY, NSM, NCM, NCC, IS, IC, NS, NREG
      INTEGER MMX,NV,K
      LOGICAL NLTEST
      DIMENSION MU(3,3), TAU(NSM,NSM), T(NSM), M(NSM), XI(NSM,NCM),
     +          NCC(NSM), DEN(NSM,NCM)
      REAL*8 BGRADP, TWOPI, FC, FT, PI, VTH,
     +       VOORF
      REAL*8 FM(MMX), MUO(3,3)
      REAL*8 XA, ZSP2, ERHOLON, Q, DUM1, DUM2, DUM3, R, Y

      NENERGY = 1
      NLTEST = .FALSE.
      WRITE(*,*)'NREG NENERGY NLTEST',NREG,' ',NENERGY,' ',NLTEST

C     CALCULATE VTH
      VTH = SQRT(2.*1.6E-16*T(IS)/M(IS))
      Q   = 1.
      R   = 1.65
      ERHOLON = 0.05
      ZSP2 = 18.

C     THE CONSTANT PI, 2 PI
      PI = 4.*ATAN(1.)
      TWOPI = 2.*PI

      FT = 1.- FC

      NV = 100
      DO 2500 K = 1, 9
        MUO(K,1) = 0.
 2500 CONTINUE
 3000 CONTINUE

      XA = DEN(2,1)*ZSP2**2 / DEN(1,1)
      IF (IS.EQ.2) XA = 1. / XA
      Y = TAU(IS,IS)/ (M(IS)*DEN(IS,1)*1E19) / (VTH/(Q*R))  

      write(*,*)'y = ', y, xa
      DUM1 = 0.53 + XA
      DUM2 = 1.77
      DUM3 = (3.02 + 4.25*XA)/(2.23 + 5.32*XA+2.40*XA**2)
      MU(1,1) = (VTH/(Q*R)) * (FT/FC) * DUM1 * Y / 
     +          (1.+2.92 * DUM1 * Y / (DUM2 * SQRT(ERHOLON)**3)) /
     +          (1.+ DUM2 * Y / (6. * DUM3))

      write(*,*) y*dum1/(1.+2.92 * DUM1 * Y / (DUM2 * SQRT(ERHOLON)**3))
     +           ,2.92 * DUM1 * Y / (DUM2 * SQRT(ERHOLON)**3)
      DUM1 = 0.71 + XA
      DUM2 = 5.32
      DUM3 = (12.43 + 20.13*XA)/(2.23 + 5.32*XA+2.40*XA**2)
      MU(1,2) = (VTH/(Q*R)) * (FT/FC) * DUM1 * Y / 
     +          (1.+2.92 * DUM1 * Y / (DUM2 * SQRT(ERHOLON)**3)) /
     +          (1.+ DUM2 * Y / (6. * DUM3))

      DUM1 = 1.39 + 13.*XA/4.
      DUM2 = 5.76
      DUM3 = (15.38 + 26.97*XA)/(2.23 + 5.32*XA+2.40*XA**2)
      MU(2,2) = (VTH/(Q*R)) * (FT/FC) * DUM1 * Y / 
     +          (1.+2.92 * DUM1 * Y / (DUM2 * SQRT(ERHOLON)**3)) /
     +          (1.+ DUM2 * Y / (6. * DUM3))
 
      MU(2,1) =  - 2.5*MU(1,1) +  MU(1,2)
      MU(1,2) = MU(2,1)
 
      VOORF = DEN(IS,1)*1E19*M(IS)
      DO 112 K = 1, 9
        MU(K,1) = VOORF*MU(K,1)
 112  CONTINUE
            
      RETURN 
      END
