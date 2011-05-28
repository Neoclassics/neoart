!--------------------------------------------------------------------
!>   subroutine that calculates the diffusion and pinch coefficients
!>
!>   input    : nar   the maximum number of species (used for 
!>                    consistency check.
!>              nzm   the maximum number of charge states (used
!>                    for consistency check)
!>              ns    the actual number of species
!>              nc    array(ns) the number of charges states per
!>                    species
!>              zsp   array(nar,nzm) real.  the charge of every
!>                    components (normalized to e)
!>              isel  integer that selects the representation of 
!>                    the geometry 1 = hamada coordinates 2 = 
!>                    circular approximation, 3 = read from file
!>              ishot shot number, only used when isel = 3
!>              ic    the contribution for which the coefficients
!>                    are calculated. 
!>                    1 banana plateau contribution
!>                    2 pfirsch schlueter contribution
!>                    3 both banana-plateau and p.s.
!>              id    the drift contribution
!>                    0 do not give temperature drift contributions vt
!>                    1 give temperature drift contribution vt
!>              m     array(ns) the mass of the species in kg
!>              t     array(ns) the temperature of the species in
!>                    kev
!>              den   array(nar,nzm) the density of every component
!>                    in units of 1.10^-19 m^-3
!>              ds    array(nar,nzm,2) the thermodynamic forces 
!>                    see the comments of the routine neoart to 
!>                    find the exact definition
!>              rho   the flux surface label
!>              eparr loop voltage 
!>              eps   the desired accuracy
!>   output     dd    array(nar,nzm) the diffusion coefficients 
!>                    of every components in units of rho^2/sec
!>              vv    array(nar,nzm) the total pinch coefficients of 
!>                    every component in units of rho/sec
!>              vt    array(nar,nzm) temperature dependent pinch      
!>                    every component in units of rho/sec
!>                    (only calculated for id=1)      
!>
!>   the routine calls the following subroutines
!>   perr    : for error handling
!>   neoart  : to calculate the fluxes
!--------------------------------------------------------------------
subroutine dandv(nar,nzm,ns,nc,zsp,isel,ishot, ic, id,    &
                 & m,t,den,ds,rho,eparr,eps,dd,vv,vt)

  implicit none

  include 'elem_config90.inc'
      
  integer, parameter :: nsm = nelmax + 2 
  integer, parameter :: ncm = nionmax 

  integer ns, nc, isel, nreg, nleg, nenergy, ncof, ic, id
  integer i, j, k, l, mt, nar, nzm, ishot
  real    m, t, den, ds, dsl, rho, eps, sigma, coeff, dd, vv
  real    zsp, vt, eparr, eparn
  logical neofrc, neogeo

  dimension m(nsm), t(nsm), den(nsm,ncm), ds(nsm,ncm,2) 
  dimension dsl(nsm,ncm,2), sigma(4), coeff(nsm,ncm,4)
  dimension dd(nsm,ncm), vv(nsm,ncm), nc(nsm), zsp(nsm,ncm)
  dimension vt(nsm,ncm) 

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

  ! first initialize d and v to zero
  do i = 1, ns ; do j = 1, nc(i)
    dd(i,j) = 0.
    vv(i,j) = 0.
  end do ; end do 

  if (id.eq.1) then
     do i=1,ns ; do j=1,nc(i)
       vt(i,j) = 0.
     end do ; end do
  endif

  ! set the epar to zero for the time being
  eparn = 0.

  ! loop over the species 
  do i = 1, ns ; do j = 1, nc(i)

    ! determine d
    do k = 1, ns ; do l = 1, nc(k) ; do mt = 1, 2
      dsl(k,l,mt) = 0.E0
    end do ; end do ; end do 
    dsl(i,j,1) = 1.E0

    call neoart(ns,nc,nar,nzm,zsp,m,t,den,dsl,rho,eps,      &
           &    isel,ishot,nreg,sigma,nleg,nenergy,ncof,     &
           &    neogeo,neofrc,ic,eparn,coeff)

    ! the flux of component i,j appears in the diffusion
    dd(i,j) = -coeff(i,j,1)/(den(i,j)*1e19)

    ! the contribution to all other fluxes appear as pinch
    do 300 k = 1, ns ; do 300 l = 1, nc(k)
      if ((k.ne.i).and.(l.ne.j)) then 
        vv(k,l) = vv(k,l) + coeff(k,l,1)*                  &
                &   (ds(i,j,1)-ds(i,j,2))/(den(k,l)*1e19)
      endif 
    end do ; end do 

    ! do not calculate a new friction and viscosity
    neogeo = .false. 
    neofrc = .true.

    ! now calculate the pinch contribution due to the temperature
    ! gradient.
    do k = 1, ns ; do l = 1, nc(k) ; do mt = 1, 2
      dsl(k,l,mt) = 0.
    end do ; end do ; end do 
    dsl(i,j,1) = ds(i,j,2)
    dsl(i,j,2) = ds(i,j,2)

    call neoart(ns,nc,nar,nzm,zsp,m,t,den,dsl,rho,eps,      &
            &       isel,ishot,nreg,sigma,nleg,nenergy,ncof, &
            &       neogeo,neofrc,ic,eparn,coeff)

    ! all these contributions appear in the pinch
    do k = 1, ns ; do l = 1, nc(k)
      vv(k,l) = vv(k,l) + coeff(k,l,1)/(den(k,l)*1e19)
    end do ; end do 

    ! the temperature dependent part is stored seperately in vt
    if (id.eq.1) then
      do  k = 1, ns
        do  l = 1, nc(k)
          vt(k,l) = vt(k,l) + coeff(k,l,1)/(den(k,l)*1e19)
        end do
      end do
    endif

  end do ; end do 

  ! add the ware pinch
  do k = 1, ns ;  do l = 1, nc(k) ; do mt = 1, 2
     dsl(k,l,mt) = 0.
  end do ; end do ; end do 


  call neoart(ns,nc,nar,nzm,zsp,m,t,den,dsl,rho,eps,       &
          &       isel,ishot,nreg,sigma,nleg,nenergy,ncof,  &
          &       neogeo,neofrc,ic,eparr,coeff)

  ! add to the pinch
  do k = 1, ns ; do l = 1, nc(k)
    vv(k,l) = vv(k,l) + coeff(k,l,1)/(den(k,l)*1e19)
  end do ; end do 


return 
end subroutine dandv 
 






