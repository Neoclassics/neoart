!--------------------------------------------------------------------
!> subroutine that calculates the heat diffusion coefficient
!>
!> input    : nar   the maximum number of species (used for 
!>                  consistency check.
!>            nzm   the maximum number of charge states (used
!>                  for consistency check)
!>            ns    the actual number of species
!>            nc    array(ns) the number of charges states per
!>                  species
!>            zsp   array(nar,nzm) real.  the charge of every
!>                  components (normalized to e)
!>            isel  integer that selects the representation of 
!>                  the geometry 1 = hamada coordinates 2 = 
!>                  circular approximation, 3 = read from file
!>            file  geoemetry file, only used when isel = 3
!>            ic    the contribution for which the coefficients
!>                  are calculated.
!>                  0 the classical particle flux            
!>                  1 banana plateau contribution
!>                  2 pfirsch schlueter contribution
!>                  3 banana-plateau, classical and p.s.
!>            m     array(ns) the mass of the species in kg
!>            t     array(ns) the temperature of the species in
!>                  kev
!>            den   array(nar,nzm) the density of every component
!>                  in units of 1.10^-19 m^-3
!>            ds    array(nar,nzm,2) the thermodynamic forces 
!>                  see the comments of the routine neoart to 
!>                  find the exact definition
!>            rho   the flux surface label     
!>            eps   the desired accuracy
!> output     chi   array(nar,nzm) of the heat diffusion coefficients 
!>                  of every components in units of rho^2/sec    
!>
!> the routine calls the following subroutines
!> perr    : for error handling
!> neoart  : to calculate the fluxes
--------------------------------------------------------------------
subroutine chis(nar,nzm,ns,nc,zsp,isel,file,ic,m,t,den,ds,rho,eps,chi)

  implicit none

  include 'elem_config.inc'
      
  integer nsm, parameter :: nsm = nelmax+2 
  integer ncm, parameter :: ncm = nionmax 

  integer  ns, nc, isel, nreg, nleg, nenergy, ncof, ic
  integer  i, j, k, l, mt, nar, nzm
  real     m, t, den, ds, dsl, rho, eps, sigma, coeff, chi
  real     zsp, eparn
  character*(*) file
  logical  neofrc, neogeo

  dimension m(nsm), t(nsm), den(nsm,ncm), ds(nsm,ncm,2)
  dimension dsl(nsm,ncm,2), sigma(4), coeff(nsm,ncm,4)
  dimension chi(nsm,ncm), nc(nsm), zsp(nsm,ncm)
         
  ! check consitency of settings
  if ((nar.ne.nsm).or.(nzm.ne.ncm)) call perr(1)

    ! always use all the full expression of the viscosity
    nreg = 0
    ! use 3 legendre polynomals
    nleg = 3
    ! use energy scattering in the collision operator
    nenergy = 1
    ! use ion - electron collisions
    ncof = 1
    ! use all couplings in the pfirsch schlueter regime
    sigma(1) = 1
    sigma(2) = 1
    sigma(3) = 1
    sigma(4) = 1
    ! on the first call recaclulate the friction/viscosity and
    ! geometry parameters
    neofrc = .false.
    neogeo = .true. 

    ! first initialize chi to zero
    do i = 1, ns 
      do j = 1, nc(i)
          chi(i,j) = 0.
      end do 
    end do 
      
    ! set the epar to zero for the time being
    eparn = 0.

    ! loop over the species 
    do i = 1, ns
      do j = 1, nc(i)

        ! determine d
        do k = 1, ns ; do l = 1, nc(k) ; do mt = 1, 2
           dsl(k,l,mt) = 0.
        end do ; end do ; end do 
        dsl(i,j,1) = 1.  !only temperature in pressure term
        dsl(i,j,2) = 1.  !temperature term

        call neoart(ns,nc,nar,nzm,zsp,m,t,den,dsl,rho,eps, &
                  & isel,file,nreg,sigma,nleg,nenergy,ncof,&
                  & neogeo,neofrc,ic,eparn,coeff)

        ! the flux of component i,j appears in the heat diffusion
        ! coeff(i,j,2) = q(i,j)/t(i)        
        chi(i,j) = -coeff(i,j,2)/(den(i,j)*1e19)

        ! do not calculate a new friction and viscosity
        neogeo = .false. 
        neofrc = .true.

      end do 
    end do 
      
return 
end subroutine chis

