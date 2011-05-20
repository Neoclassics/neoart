     



      SUBROUTINE PERR(N)
C--------------------------------------------------------------------
C     THIS SUBROUTINE DOES THE ERROR HANDLING
C--------------------------------------------------------------------

      IMPLICIT NONE
 
      INTEGER N,IFILE, IWARN1, IWARN2, IWARN3, IWARN4, IWARN5, 
     +        IWARN6, IWARN7
      LOGICAL NLWARN
      
      DATA IFILE / 0 /
      DATA IWARN1, IWARN2, IWARN3, IWARN4, IWARN5, IWARN6, IWARN7
     +    /     0,      0,      0,      0,      0,      0,      0/
      
      NLWARN = .FALSE. 
      IF (N.LE.10) THEN
        GOTO (1,2,3,4,5) N
      ELSE
C       TO SWITCH OFF ALL WARNINGS ACTIVATE THE LINE BELOW
        IF (.NOT.NLWARN) RETURN
        IF (IFILE.EQ.0) THEN
          OPEN(11,FILE = 'NEOWARN.TXT')
          IFILE = 1
        ENDIF
        GOTO (100,101,102,103,104,105,106)N-10
      ENDIF
      RETURN

C-------------------------------------------------------------------
C     THE ERROR MESSAGES
C-------------------------------------------------------------------
 1    CONTINUE
      WRITE(*,*)'THE PARAMETERS NSM AND/OR NCM ARE NOT CONSISTENTLY'
      WRITE(*,*)'SET THROUGHOUT THE CODE. THEY APPEAR IN DIFFERENT'
      WRITE(*,*)'SUBROUTINES BUT MUST ALL HAVE THE SAME VALUE'
      STOP

 2    CONTINUE
      WRITE(*,*)'SORRY, YOU CAN NOT RUN THE CODE FOR RHO = 0.'
      STOP

 3    CONTINUE 
      WRITE(*,*)'THE PARAMETER NMAXGR IS NOT CONSISTENTLY'
      WRITE(*,*)'SET THROUGHOUT THE CODE. THEY APPEAR IN DIFFERENT'
      WRITE(*,*)'SUBROUTINES BUT MUST HAVE THE SAME VALUE'
      STOP

 4    CONTINUE
      WRITE(*,*)'UNABLE TO REACH CONVERGENCE IN VISCOS'
      STOP

 5    CONTINUE
      WRITE(*,*)'UNABLE TO REACH CONVERGENCE IN THE ROUTINE GEOM'
      STOP
      
C-------------------------------------------------------------------
C     THE WARNING MESSAGES
C-------------------------------------------------------------------
 100  CONTINUE
      IF (IWARN1.EQ.0) THEN
        WRITE(11,*)'WARNING'
        WRITE(11,*)'YOU HAVE USED NREG.NE.0. THEREFORE A CERTAIN',
     +            ' COLLISION'
        WRITE(11,*)'REGIME WILL BE FORCED. THE RESULTS MIGHT BE ',
     +            'WRONG'
        WRITE(11,*)
        IWARN1 = 1
      ENDIF 
      RETURN

 101  CONTINUE
      IF (IWARN2.EQ.0) THEN
        WRITE(11,*)'WARNING'
        WRITE(11,*)'YOU HAVE USED NLEG.NE.3. THE NUMBER OF LEG',
     +            'ENDRE HAR-'
        WRITE(11,*)'MONICS IS NOT EQUAL 3. THIS CAN MAKE THE ',
     +            'RESULTS INACCURATE'
        WRITE(11,*)
        IWARN2 = 1
      ENDIF
      RETURN

 102  CONTINUE
      IF (IWARN3.EQ.0) THEN
        WRITE(11,*)'WARNING'
        WRITE(11,*)'YOU HAVE USED NENERGY.NE.1 THE ENERGY SCAT',
     +            'TERING IN'
        WRITE(11,*)'THE CALCULATION OF THE VISCOSITY IS NOT TA',
     +            'KEN INTO ACCOUNT'
        WRITE(11,*)
        IWARN3 = 1
      ENDIF
      RETURN

 103  CONTINUE
      IF (IWARN4.EQ.0) THEN
        WRITE(11,*)'WARNING'
        WRITE(11,*)'YOU HAVE USED NCOFC.NE.1. ION-ELECTRON COL',
     +            'LISIONS ARE'
        WRITE(11,*)'NOT PROPERLY ACCOUNTED FOR'
        WRITE(11,*)
        IWARN4 = 1
      ENDIF
      RETURN

 104  CONTINUE
      IF (IWARN5.EQ.0) THEN
        WRITE(11,*)'WARNING'
        WRITE(11,*)'ONE OF THE SIGMAS DOES NOT HAVE THE VALUE ',
     +             '1. AND ONE'
        WRITE(11,*)'OF THE COUPLINGS IN THE PFIRSCH SCHLUETER ',
     +             'CALCULATION IS'
        WRITE(11,*)'THEREFORE TAKEN INTO ACCOUNT'
        WRITE(11,*)
        IWARN5 = 1
      ENDIF
      RETURN

 105  CONTINUE
      IF (IWARN6.EQ.0) THEN
        WRITE(11,*)'WARNING'
        WRITE(11,*)'BOTH NEOGEO AND NEOFRC ARE SET TRUE. HOWEVER',
     +             ', THE VISCOSITY'
        WRITE(11,*)'COEFFICIENTS CONTAIN GEOMETRIC DEPENDENCIES',
     +             ' AND MUST BE'
        WRITE(11,*)'NEWLY CALCULATED'
        WRITE(11,*)
        IWARN6 = 1
      ENDIF
      RETURN
 
 106  CONTINUE
      IF (IWARN7.EQ.0) THEN
        WRITE(11,*)'WARNING'
        WRITE(11,*)'THE ROUTINE NEOART HAS BEEN CALLED WITH DIF',
     +             'FERENT RHO/ISEL'
        WRITE(11,*)'MORE THAN NKEEP TIMES. THE VALUES OF THE ',
     +             'GEOMETRY DEPENDENT'
        WRITE(11,*)'QUANTITIES CAN NO LONGER BE KEPT IN CORE'
        WRITE(11,*)
        IWARN7 = 1
      ENDIF
      RETURN
  
      END
