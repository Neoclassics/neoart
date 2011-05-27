

      SUBROUTINE CIRCGEOM(IT,RHO,RN,E,Q,BN)
C--------------------------------------------------------------------
C     THIS ROUTINE SUPLIES THE PARAMETERS THAT DESCRIBE THE CIRCU-
C     LAR EQUILIBRIUM. IN PRINCIPLE IT IS THE USER RESPONCIBILITY
C     TO GIVE THE CORRECT VALUES OF ALL THE QUANTITIES. IT IS HERE
C     USED IN A SIMPLE WAY. THE DESIRED PARAMETERS ARE SAVED AND 
C     THE CALLS FROM GEOM THAN READ OUT THE QUANTITIES
C   
C     INPUT  :  IT       : CONTROL PARAMETER (INTEGER) IF 
C                          1 = QUANTITIES ARE SAVED
C                          2 = READ OUT OF SAVED QUANTITIES
C                          THE ROUTINE SHOULD BE CHANGED SUCH
C                          THAT ALL QUANTITIES ARE CALCULATED 
C                          AS A FUNCTION OF RHO
C               RHO      : FLUX SURFACE LABEL.
C     OUTPUT :  RN       : THE MAJOR RADIUS OF THE MAGNETIC AXIS
C               E        : THE INVERSE ASPECT RATIO OF THE SURFACE
C               Q        : THE SAFETY FACTOR OF THE SURFACE
C               BN       : THE MAGNETIC FIELD STRENGTH IN CHI = 
C                          PI / 2. WITH CHI BEING THE POLOIDAL 
C                          ANGLE
C--------------------------------------------------------------------

      IMPLICIT NONE

      REAL E, RN, Q, BN, RHO, EK, RNK, QK, BNK
      INTEGER IT

      SAVE EK, RNK, QK, BNK

      IF (IT .EQ.1) THEN
        EK  = E
        RNK = RN
        QK  = Q
        BNK = BN
      ELSE
        E  = EK
        RN = RNK
        Q  = QK
        BN = BNK
      ENDIF

      RETURN
      END
