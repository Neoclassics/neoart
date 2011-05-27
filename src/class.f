      SUBROUTINE CLASS(NSM,NCM,NS,NC,T,ZSP,XI,LA,LAB,DS,
     +                 COEFFC)
C-----------------------------------------------------------
C     THIS SUBROUTINE CALCULATES THE CLASSICAL CONTRIBUTION
C     TO THE TRANSPORT
C 
C     INPUT  NSM  : MAXIMUM NUMBER OF SPECIES
C            NCM  : MAXIMUM NUMBER OF CHARGE STATES
C            NS   : ACTUAL NUMBER OF SPECIES
C            NC   : ARRAY(NSM) ACTUAL NUMBER OF CHARGE 
C                   STATES PER SPECIES
C            T    : ARRAY(NSM) TEMPERATURE OF EVERY SPECIES
C            ZSP  : ARRAY(NSM,NCM) THE CHARGE IN UNITS OF E
C                   OF EVERY CHARGE STATE
C            XI   : ARRAY(NSM,NCM) WEIGTING FACTOR FOR THE 
C                   REDUCED CHARGE STATE METHOD
C            LA   : ARRAY(3,3,NSM) TEST PARTICLE FRICTION 
C                   COEFFICIENTS
C            LAB  : ARRAY(3,3,NSM,NSM) FRICTION COEFFICIENTS
C            DS   : ARRAY(NSM,NCM,2) THE THERMODYNAMIC
C                   FORCES. FOR DEFINITION SEE ROUTINE 
C                   NEOART
C 
C     OUTPUT COEFFC : THE TRANPORT COEFFICIENTS (SEE NEOART
C                     FOR DEFINITION
C-----------------------------------------------------------

      IMPLICIT NONE

      INTEGER NSM,NCM,NS,NC,I,J,K,L
      REAL  XI,LA,LAB,DS,COEFFC,V,ZSP,T

      DIMENSION NC(NSM), XI(NSM,NCM), LA(3,3,NSM),
     +          LAB(3,3,NSM,NSM), DS(NSM,NCM,2),
     +          COEFFC(NSM,NCM,4),ZSP(NSM,NCM),
     +          T(NSM)

      DO 100 I = 1, NS
        DO 100 J = 1, NC(I)
          
          COEFFC(I,J,1) = 0.
          COEFFC(I,J,2) = 0.
 
          DO 200 K = 1, NS
            DO 200 L = 1, NC(K)

              V = T(K)/ZSP(K,L)
              COEFFC(I,J,1) = COEFFC(I,J,1) + V*XI(I,J)*
     +         XI(K,L)*(LAB(1,1,I,K)*DS(K,L,1)+LAB(1,2,I,K)*
     +           DS(K,L,2)) 
              COEFFC(I,J,2) = COEFFC(I,J,2) + V*XI(I,J)*
     +         XI(K,L)*(LAB(2,1,I,K)*DS(K,L,1)+LAB(2,2,I,K)*
     +           DS(K,L,2)) 
              IF ((K.EQ.I).AND.(L.EQ.J)) THEN 
                COEFFC(I,J,1) = COEFFC(I,J,1) + V*XI(I,J)*(
     +            LA(1,1,I)*DS(I,J,1)+LA(1,2,I)*DS(I,J,2))
                COEFFC(I,J,2) = COEFFC(I,J,2) + V*XI(I,J)*(
     +            LA(2,1,I)*DS(I,J,1)+LA(2,2,I)*DS(I,J,2))
              ENDIF

 200      CONTINUE
 100  CONTINUE

      RETURN
      END
