      PROGRAM TEST
C-----------------------------------------------------------------
C     CALCULATE THE TRANSPORT COEFFICIENTS IN THE BANANA REGIME 
C     FOR A SIMPLE HYDROGEN PLASMA AND COMPARE THEM WITH THE 
C     VALUES GIVEN IN THE LITERATURE.
C-----------------------------------------------------------------
c     CALL TEST1

C-----------------------------------------------------------------
C     CALCULATE THE TRANSPORT COEFFICIENTS IN THE P.S. REGIME 
C     FOR A SIMPLE HYDROGEN PLASMA AND COMPARE THEM WITH THE 
C     VALUES GIVEN IN THE LITERATURE.
C-----------------------------------------------------------------
c     CALL TEST2
 
C-----------------------------------------------------------------
C     TEST THE REDUCED CHARGE STATE METHOD
C-----------------------------------------------------------------
C     CALL TEST3

C-----------------------------------------------------------------
C     CALCULATE THE TRANSPORT COEFFICIENTS OF AN PLASMA CONTAINING
C     OF HYDROGEN AND OXYGEN, AND COMPARE THEM WITH THE ANALYTIC 
C     EXPRESSIONS OF HIRSHMANN AND SIGMAR.
C-----------------------------------------------------------------
C     CALL TEST4

C-----------------------------------------------------------------
C     CALCULATE THE HEAT TRANSPORT OF A PLASMA CONTAINING
C     OF HYDROGEN AND OXYGEN, AND COMPARE THEM WITH THE ANALYTIC 
C     EXPRESSIONS OF HIRSHMANN AND SIGMAR, AND THE DERIVED 
C     EXPRESSION KEEPING FINITE  MASS RATIO EFFECTS.
C-----------------------------------------------------------------
C     CALL TEST5

C-----------------------------------------------------------------
C     EXAMPLE: THE ION HEAT CONDUCTIVITY OF A PURE PLASMA: THE 
C     TRANSITION FROM THE COLLISIONAL CASE TO THE COLLISIONLESS
C     CASE.
C-----------------------------------------------------------------
C     CALL TEST6

C-----------------------------------------------------------------
C     EXAMPLE: THE SCREENING OF CARBON IN AN HYDROGEN PLASMA: 
C     THE TRANSITION FROM THE COLLISIONAL TO THE COLLISIONLESS
C     CASE. 
C-----------------------------------------------------------------
C     CALL TEST7

C-----------------------------------------------------------------
C     CALCULATE THE HEAT TRANSPORT OF A PLASMA CONTAINING
C     OF HYDROGEN AND A HEAVY IMPURITY, AND COMPARE THEM WITH 
C     THE ANALYTIC EXPRESSIONS OF HIRSHMANN AND SIGMAR.
C     HERE ALSO ENERGY COUPLING IS CONSIDERED. 
C-----------------------------------------------------------------
C     CALL TEST8

C-----------------------------------------------------------------
C     CALCULATE THE HEAT TRANSPORT OF A PLASMA CONTAINING
C     OF HYDROGEN AND OXYGEN, AND COMPARE THEM WITH THE ANALYTIC 
C     EXPRESSIONS OF HIRSHMANN AND SIGMAR.
C     ALSO ENERGY COUPLING IS CONSIDERED HERE.
C-----------------------------------------------------------------
c     CALL TEST9
 
C-----------------------------------------------------------------
C     CALCULATE THE HEAT TRANSPORT OF A PLASMA CONTAINING
C     OF HYDROGEN AND OXYGEN, AND COMPARE THEM WITH THE ANALYTIC 
C     EXPRESSIONS OF HIRSHMANN AND SIGMAR.
C     HERE THE HYDROGEN HEAT FLUX IS CALCULATED AS A FUNCTION 
C     OF THE IMPURITY CONTENT (GIVEN BY ALPHA) FOR DIFFERENT 
C     CHARGE STATES OF THE OXYGEN.
C-----------------------------------------------------------------
C     CALL TEST10

C-----------------------------------------------------------------
C     CALCULATE THE HEAT TRANSPORT OF A PLASMA CONTAINING
C     OF HYDROGEN AND OXYGEN, AND COMPARE THEM WITH THE ANALYTIC 
C     EXPRESSIONS OF HIRSHMANN AND SIGMAR.
C     HERE THE OXYGEN HEAT FLUX IS CALCULATED AS A FUNCTION 
C     OF THE IMPURITY CONTENT (GIVEN BY ALPHA) FOR DIFFERENT 
C     CHARGE STATES OF THE OXYGEN.
C-----------------------------------------------------------------
c     CALL TEST11

C-----------------------------------------------------------------
C     CALCULATE THE EPSILON (INVERSE ASPECT RATIO) DEPENDENCE 
C     OF THE ION HEAT FLUX IN A PURE PLASMA
C-----------------------------------------------------------------
c     CALL TEST12

C-----------------------------------------------------------------
C     COMPARE THE 2 AND 3 LAGUERRE POLYNOMAL EXPANSIONS OF THE 
C     ION HEAT FLUX AT FINITE INVERSE ASPECT RATIO.
C-----------------------------------------------------------------
c     CALL TEST13
 
C-----------------------------------------------------------------
C     THE TOTAL OXYGEN PARTICLE FLUX AS A FUNCTION OF NORMALIZED
C     COLLISIONALITY
C-----------------------------------------------------------------
c     CALL TEST14  

C-----------------------------------------------------------------
C     THE PARTICLE FLUX OF A TRACE IMPURITY IN AN HYDROGEN 
C     CARBON PLASMA
C-----------------------------------------------------------------
c     CALL TEST15

C-----------------------------------------------------------------
C     TEST THE DANDV ROUTINE. A PURE PLASMA IN THE BANANA 
C     PLATEAU REGIME 
C-----------------------------------------------------------------
C     CALL TEST16

C-----------------------------------------------------------------
C     TEST THE HAMADA COORDINATE VERSION OF THE GEOMETRY DEPEN-
C     DENT QUANTITIES
C-----------------------------------------------------------------
      CALL TEST17

C-----------------------------------------------------------------
C     THE ASDEX UPGRADE IMPURITY TRANSPORT EXPERIMENTS. 
C-----------------------------------------------------------------
      do 1 i = 1, 100
c     CALL RALF
1     continue
 
      
      END
