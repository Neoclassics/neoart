
      SUBROUTINE TEST15
C-----------------------------------------------------------------
C     THE PARTICLE FLUX OF A TRACE IMPURITY IN AN HYDROGEN 
C     CARBON PLASMA
C-----------------------------------------------------------------
      implicit none

      integer nsm, ncm
      INTEGER NMAXGR, ISEL

      PARAMETER(NMAXGR = 100000)
      parameter(nsm = 5)
      parameter(ncm = 30)

      integer i, j, k, itel, ns, nc, 
     +        it, IC, NLEG, NREG, NCOF, NENERGY
      real*8 tau, xi, sigma,  mas, temp, ZSP,
     +       den, rho, coeff, EPS, Q, RN, BN, norm, alpha,
     +       l1,l2,l3,l11tt,l11th,l12th,DS,e
      LOGICAL NEOFRC, NEOGEO

      dimension tau(nsm,nsm), xi(nsm,ncm), sigma(4),
     +          mas(nsm), temp(nsm),
     +          den(nsm,ncm), DS(nsm,ncm,2),
     +          zsp(nsm,ncm), coeff(nsm,ncm,4)
      dimension nc(nsm)

      data itel / -1 /

C     THE GEOMETRY IS DETERMINED BY THE FOLLOWING QUANTITIES
C     RHO IS THE INVERSE ASPECT RATIO 
      RHO = 0.3
      E = 0.3
C     Q IS THE SAFETY FACTOR
      Q = 4.
C     BN IS THE TOTAL MAGNETIC FIELD AT THE LOCATION CHI = PI/2 
      BN = 2.5 
C     RN IS THE MAJOR RADIUS OF THE MAGNETIC AXIS
      RN = 1.65
C     COPY THEM INTO THE VALUES USED BY THE CODE
      CALL CIRCGEOM(1,RHO,RN,E,Q,BN)

C     SET THE COUPLING
      sigma(1) = 0
      sigma(2) = 0
      sigma(3) = 1
      sigma(4) = 1

C     THE NUMBER OF SPECIES
      NS = 3
C     THE NUMBER OF CHARGE STATES PER SPECIES
      NC(1) = 1
      NC(2) = 1
      NC(3) = 1

      do 1111 it = 1, 21
        alpha = exp((-3.+3.*(it-1)/20.)*log(10.))


C     THE MASS OF THE SPECIES (IN KG) 
      MAS(1) = 1.6727E-27
c      MAS(1) = 1e-30*1.6727E-27
      MAS(2) = 12*1.6727E-27
      MAS(3) = 20*1.6727E-27

C     THE CHARGE 
      ZSP(1,1) = 1.
      ZSP(2,1) = 4.
      ZSP(3,1) = 11.

C     THE DENSITY IN UNITS OF 10^19 M^-3
      DEN(1,1) = 7
      DEN(2,1) = alpha*DEN(1,1)
      DEN(3,1) = 1E-8


C     THE TEMPERATURE IN UNITS OF KEV
      TEMP(1) = 0.001
      TEMP(2) = 0.001
      TEMP(3) = 0.001

C     USE ELECTRON-ION COLLISIONS
      NCOF = 1

C     CALCULATE THE PFIRSCH SCHLUETER CONTRIBUTION
      IC = 2

      NEOGEO = .TRUE.
      NEOFRC = .FALSE. 

C     SET THE THERMODYNAMIC FORCES
      do 2000 i = 1, ns
        do 2000 j = 1, nc(i)
          do 2000 k = 1, 2 
            DS(I,J,K) = 0.
 2000 continue
      DS(2,1,2) = 1.
      
      ISEL = 2
      EPS = 1E-5

      CALL NEOART(NS,NC,NSM,NCM,ZSP,MAS,TEMP,DEN,DS,RHO,EPS,
     +            ISEL,NREG,SIGMA,NLEG,NENERGY,NCOF,
     +            NEOGEO,NEOFRC,IC,COEFF)


C     NOW CALCULATE THE NORMALIZED COLLISION FREQUENCIES
      CALL COLXI(NSM,NCM,NS,NC,ZSP,DEN,TEMP,MAS,TAU,XI)

      NORM = 1.6E-22*ZSP(3,1)*BN**2/(2.*Q**2*TEMP(1))

      call psflux(temp(1),nsm,ncm,mas,den,zsp,tau,
     +                  l1,l2,l3,l11tt,l11th,l12th)


c      write(*,301)alpha,-coeff(3,1,1)/temp(1),
c     +          (l1+l2+l3)/zsp(2,1)-l11th+l12th
c      write(*,301)den(2,1)/den(1,1),-zsp(2,1)*coeff(3,1,1)*norm,
c     +            l3
c      write(*,301)temp(1),-zsp(2,1)*coeff(3,1,1)*norm
c     +            l3
c      write(*,301)den(2,1)/den(1,1),zsp(2,1)*coeff(3,1,1)*norm,
c     +            -l1-l2
c      write(*,301)den(2,1)/den(1,1),coeff(3,1,1)*norm,
c     +             l11th
c      write(*,301)den(2,1)/den(1,1),-zsp(3,1)*coeff(3,1,1)*norm,
c     +             l11tt
      write(*,301)den(2,1)/den(1,1),-coeff(3,1,1)*norm,
     +             l12th


 301  format(3(1x,1pe13.5))
 1111 continue
      return
      end
