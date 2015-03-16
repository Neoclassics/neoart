# Introduction #

Neoart computes the neoclassical fluxes as outputs given the thermodynamic forces as inputs.  In order to compute neoclassical fluxes the following parameters are necessary
  * temperature and density values and their gradient length scales for all species
  * magnetic geometry
  * loop voltage

The most convenient way to define inputs is to use the inverse aspect ratio as the flux surface label.  Then the gradient length scales are:
```
1 / L = 1/ A d(A)/d(eps) 
```

# Simple example #

Two species (electron and proton) with equal temperature and density. The
only non-zero thermodynamic force is the ion temperature gradient.  The  fluxes are evaluated in the banana-plateau regime.  Circular geometry is used ([issue 11](https://code.google.com/p/neoart/issues/detail?id=11)).

## Fortran 90 ##

```
! set the parameters for the circular geometry
rho = 0.05
e = 1e-4
q = 2.
rn = 1.65
bn = 2.5
!copy them into the values used by the code
call circgeom(1,rho,rn,e,q,bn)
!set the electric field to zero 
eparr = 0. 
!use the circular geometry approximation
isel = 2
!set the accuracy
eps = 1e-5
!force the banana regime in the calculation of the transport
!fluxes (in the calculation of the viscosity)
nreg = 1
!the number of legendre harmonics
nleg = 3
!use energy scattering in the calc. of viscosity
nenergy = 1
!switch off ion-electron collisions
ncof = 0
!in the first call the matrices have to be calculated
neofrc = .false.
!calculate the banana plateau contribution
ic = 1

!the number of species is 2
ns = 2
!ions and electrons have only one charge 
nc(1) = 1
nc(2) = 1
!the mass of the electron and proton
m(1) = 9.1096e-31
m(2) = 1.6727e-27
!the charge of the electron and ion
zsp(1,1) = -1
zsp(2,1) = 1
!the density of the species in 10^19 m^-3
den(1,1) = 5.
den(2,1) = 5. 
!the temperature in kev
t(1) = 10.
t(2) = 10.

!the thermodynamic forces      
do i = 1, ns
  do j = 1, nc(i)
    ds(i,j,1) = 0.
    ds(i,j,2) = 0.
  enddo
enddo
ds(1,1,1) = 1.

call neoart(ns,nc,nar,nzm,zsp,m,t,den,ds,rho,eps, &
            isel,ishot,nreg,sigma,nleg,nenergy,ncof, &
            neogeo,neofrc,ic,eparr,cff1)
```

## Python ##

```
>>> nc = [1, 1] # number of ionisation stages per species
>>> zsp = [[-1], [1]] # charge for each ionisation stage
>>> mass = [9.1096e-31, 1.6727e-27] # kg
>>> density = [[5], [5]] # 1e19 m^-3
>>> temperature = [10, 10] # keV
>>> forces = np.zeros((2,1,2))
>>> forces[1,0,1] = 1 # ion temperature gradient
>>> g = dict(rho=0.000165, Rmag=1.65, q=2, bn=2.5)
>>> coeff = neoart.neoart(nc, zsp, mass, temperature, density, forces,
        geometry=g, contribution='bp')
```

## Matlab ##

```
ns = 2; % number of species
nc = [1, 1]; % number of ionisation stages per species
zsp = [-1; 1]; % charge for each ionisation stage
mass = [9.1096e-31, 1.6727e-27]; % kg
density = [5; 5]; % 1e19 m^-3
temperature = [10, 10];
forces = zeros(2,1,2);
forces(2,1,2) = 1; % ion temperature  gradient
isel = 2;
rho = 0.000165;
RN = 1.65;
BN = 2.5;
q = 2;
ic = 2;
eparr = 0;
[c1, c2, c3, c4] = neoart(ns, nc, mass, temperature, zsp, density, ...
                    forces, isel, rho, RN, BN, q, ic, eparr);
```


# Contribution #

If you want to contribute to Neoart make sure to conform the following coding standards.

## Code layout ##
  1. Use Fortran 90. (see [issue 4](https://code.google.com/p/neoart/issues/detail?id=4))
  1. Name source files lowercase with an extension `.f90`. Example:
```
   foo.f90
   foo_bar.f90
```
  1. Use 2 spaces for indention. Do not mix tabs and spaces.
  1. Limit the lines to 69 characters.
  1. Do not leave whitespaces at the end of the lines.
  1. Use lowercase variable names.

## Subroutine headers ##

Comment every subroutine, explain briefly what it does, what are the input and output variables. Example from
[menn.f90](http://code.google.com/p/neoart/source/browse/trunk/src/menn.f90):
```
!--------------------------------------------------------------------
!> subroutine that calculates the friction coefficients m^ab and 
!> n^ab 
!>
!> input    ta   : the temperature of species a
!>          tb   : the temperature of species b
!>          ma   : the mass of species a
!>          mb   : the mass of species b
!> output   mm   : array(3,3) containing the test particle friction
!>                 coefficients
!>          nn   : array(3,3) containing the field particle friction
!>                 coefficients
!--------------------------------------------------------------------
subroutine menn(ta, tb, ma, mb, mm, nn)
```