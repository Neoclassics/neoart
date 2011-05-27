      PROGRAM TEST
C-----------------------------------------------------------------
C     CALCULATE THE TRANSPORT COEFFICIENTS IN THE BANANA REGIME 
C     FOR A SIMPLE HYDROGEN PLASMA AND COMPARE THEM WITH THE 
C     VALUES GIVEN IN THE LITERATURE.
C-----------------------------------------------------------------
      WRITE(*,*)'START TEST1   '
      CALL TEST1

C-----------------------------------------------------------------
C     CALCULATE THE TRANSPORT COEFFICIENTS IN THE P.S. REGIME 
C     FOR A SIMPLE HYDROGEN PLASMA AND COMPARE THEM WITH THE 
C     VALUES GIVEN IN THE LITERATURE.
C-----------------------------------------------------------------
      WRITE(*,*)'START TEST2   '
      CALL TEST2
 
C-----------------------------------------------------------------
C     TEST THE REDUCED CHARGE STATE METHOD
C-----------------------------------------------------------------
      WRITE(*,*)'START TEST3   '
      CALL TEST3

C-----------------------------------------------------------------
C     CALCULATE THE TRANSPORT COEFFICIENTS OF AN PLASMA CONTAINING
C     OF HYDROGEN AND OXYGEN, AND COMPARE THEM WITH THE ANALYTIC 
C     EXPRESSIONS OF HIRSHMANN AND SIGMAR.
C-----------------------------------------------------------------
      WRITE(*,*)'START TEST4   '
      CALL TEST4

C-----------------------------------------------------------------
C     CALCULATE THE HEAT TRANSPORT OF A PLASMA CONTAINING
C     OF HYDROGEN AND OXYGEN, AND COMPARE THEM WITH THE ANALYTIC 
C     EXPRESSIONS OF HIRSHMANN AND SIGMAR, AND THE DERIVED 
C     EXPRESSION KEEPING FINITE  MASS RATIO EFFECTS.
C-----------------------------------------------------------------
      WRITE(*,*)'START TEST5   '
      CALL TEST5

C-----------------------------------------------------------------
C     EXAMPLE: THE ION HEAT CONDUCTIVITY OF A PURE PLASMA: THE 
C     TRANSITION FROM THE COLLISIONAL CASE TO THE COLLISIONLESS
C     CASE.
C-----------------------------------------------------------------
      WRITE(*,*)'START TEST6   '
      CALL TEST6

C-----------------------------------------------------------------
C     EXAMPLE: THE SCREENING OF CARBON IN AN HYDROGEN PLASMA: 
C     THE TRANSITION FROM THE COLLISIONAL TO THE COLLISIONLESS
C     CASE. 
C-----------------------------------------------------------------
      WRITE(*,*)'START TEST7   '
      CALL TEST7

C-----------------------------------------------------------------
C     CALCULATE THE HEAT TRANSPORT OF A PLASMA CONTAINING
C     OF HYDROGEN AND A HEAVY IMPURITY, AND COMPARE THEM WITH 
C     THE ANALYTIC EXPRESSIONS OF HIRSHMANN AND SIGMAR.
C     HERE ALSO ENERGY COUPLING IS CONSIDERED. 
C-----------------------------------------------------------------
      WRITE(*,*)'START TEST8   '
      CALL TEST8

C-----------------------------------------------------------------
C     CALCULATE THE HEAT TRANSPORT OF A PLASMA CONTAINING
C     OF HYDROGEN AND OXYGEN, AND COMPARE THEM WITH THE ANALYTIC 
C     EXPRESSIONS OF HIRSHMANN AND SIGMAR.
C     ALSO ENERGY COUPLING IS CONSIDERED HERE.
C-----------------------------------------------------------------
      WRITE(*,*)'START TEST9   '
      CALL TEST9
 
C-----------------------------------------------------------------
C     CALCULATE THE HEAT TRANSPORT OF A PLASMA CONTAINING
C     OF HYDROGEN AND OXYGEN, AND COMPARE THEM WITH THE ANALYTIC 
C     EXPRESSIONS OF HIRSHMANN AND SIGMAR.
C     HERE THE HYDROGEN HEAT FLUX IS CALCULATED AS A FUNCTION 
C     OF THE IMPURITY CONTENT (GIVEN BY ALPHA) FOR DIFFERENT 
C     CHARGE STATES OF THE OXYGEN.
C-----------------------------------------------------------------
      WRITE(*,*)'START TEST10  '
      CALL TEST10

C-----------------------------------------------------------------
C     CALCULATE THE HEAT TRANSPORT OF A PLASMA CONTAINING
C     OF HYDROGEN AND OXYGEN, AND COMPARE THEM WITH THE ANALYTIC 
C     EXPRESSIONS OF HIRSHMANN AND SIGMAR.
C     HERE THE OXYGEN HEAT FLUX IS CALCULATED AS A FUNCTION 
C     OF THE IMPURITY CONTENT (GIVEN BY ALPHA) FOR DIFFERENT 
C     CHARGE STATES OF THE OXYGEN.
C-----------------------------------------------------------------
      WRITE(*,*)'START TEST11  '
      CALL TEST11

C-----------------------------------------------------------------
C     CALCULATE THE EPSILON (INVERSE ASPECT RATIO) DEPENDENCE 
C     OF THE ION HEAT FLUX IN A PURE PLASMA
C-----------------------------------------------------------------
      WRITE(*,*)'START TEST12  '
      CALL TEST12

C-----------------------------------------------------------------
C     COMPARE THE 2 AND 3 LAGUERRE POLYNOMAL EXPANSIONS OF THE 
C     ION HEAT FLUX AT FINITE INVERSE ASPECT RATIO.
C-----------------------------------------------------------------
      WRITE(*,*)'START TEST13  '
      CALL TEST13
 
C-----------------------------------------------------------------
C     THE TOTAL OXYGEN PARTICLE FLUX AS A FUNCTION OF NORMALIZED
C     COLLISIONALITY
C-----------------------------------------------------------------
      WRITE(*,*)'START TEST14  '
      CALL TEST14  

C-----------------------------------------------------------------
C     THE PARTICLE FLUX OF A TRACE IMPURITY IN AN HYDROGEN 
C     CARBON PLASMA
C-----------------------------------------------------------------
C      WRITE(*,*)'START TEST15  '
C      CALL TEST15

C-----------------------------------------------------------------
C     TEST THE DANDV ROUTINE. A PURE PLASMA IN THE BANANA 
C     PLATEAU REGIME 
C-----------------------------------------------------------------
      WRITE(*,*)'START TEST16  '
      CALL TEST16

C-----------------------------------------------------------------
C     TEST THE HAMADA COORDINATE VERSION OF THE GEOMETRY DEPEN-
C     DENT QUANTITIES
C-----------------------------------------------------------------
C      WRITE(*,*)'START TEST17  '
C      CALL TEST17

c-----------------------------------------------------------------
C     TEST THE CONTRIBUTION FROM THE LOOP VOTAGE 
C-----------------------------------------------------------------
      WRITE(*,*)'START TEST18'
      CALL TEST18

      END
