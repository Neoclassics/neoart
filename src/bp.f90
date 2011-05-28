!--------------------------------------------------------------------
!> input  nsc ncc : the parameters nsm and ncm used to check
!>                  consistency
!>        ns      : the number of species
!>        nc      : array(nsm) the number of charge states
!>                  of every species.
!>        neofrc  : logical if true the matrices are not 
!>                  newly calculated. 
!>        t       : array(nsm) the temperature of the species
!>                  in kev
!>        m       : array(nsm) the mass of the species in kg.
!>        den     : array(nsm,ncm) density of every component
!>                  in units 10^19 m^-3
!>        zsp     : array(nsc,ncc) the charge of every comp.
!>                  normalized to e 
!>        tau     : array(nsm,nsm) the collision freq. 
!>                  m_i n_i / tau_ij
!>        xi      : array(nsc,ncc) the relative weigth of 
!>                  every charge state. 
!>        la      : the friction coefficients
!>        lab     : the friction coefficients
!>        eps     : desired accuracy
!>        nenergy : integer if 1 energy scattering is used
!>                  in the calculation of the viscosity
!>        nleg    : the number of laguerre polynomals 
!>        nreg    : a parameter that forces a certain regime
!>                  in the calculation of the viscosity
!>                  0 wheigthed regime valid for all collision
!>                    frequencies
!>                  1 force the banana regime
!>                    other force the pfirsch schlueter regime
!>        bgradp  : the flux surface quantity n in grad theta
!>                  the inner product of the unit vector in 
!>                  the direction of the magnetic field and 
!>                  the gradient of the poloidal angle
!>        fc      : the number of passing particles
!>        fm      : the coefficients necessary to calculate
!>                  the viscosity in the pfirsch schlueter 
!>                  regime. 
!>        mmx     : the number of fm coefficients
!>        dsi     : the thermodynamic force, for a definition
!>                  see routine neoart
!> output coeff   : normalized transport coefficients
!>
!> the routine calls the following subroutines
!> perr    : for error handling
!> viscos  : to calculate the viscosity
!> ludcmp  : to do the lu decomposition
!> lubksb  : to solve the matrix equation
!-------------------------------------------------------------------
subroutine bp(nsc,ncc,ns,nc,neofrc,t,m,den,zsp,tau,xi,la,      & 
     &     lab,eps,nenergy,nleg,nreg,bgradp,fc,fm,mmx,dsi,epar, &
     &     coeff)

  implicit none

  include 'elem_config90.inc'
      
  ! maximum number of species and charge states
      
  integer, parameter :: nsm = nelmax + 2 
  integer, parameter :: ncm = nionmax 

  integer  ns, nc, i, j, mt, nreg, indx, nsc, ncc, nleg, mmx, k, l
  real t, m, den, tau, aaa, d, akp, la, lab, visc, xi, fc, fm 
  real bgradp,zsp, epar, ainv, ds, coeff, bt, upar,dsi

  dimension nc(nsm),t(nsm),m(nsm),den(nsm,ncm),upar(nsm,ncm,3)
  dimension tau(nsm,nsm), visc(3,3,nsm,ncm), fm(mmx), xi(nsm,ncm)
  dimension aaa(3,3,nsm,nsm), akp(3,3),zsp(nsm,ncm)
  dimension d(3*nsm),indx(3*nsm), la(3,3,nsm),lab(3,3,nsm,nsm)
  dimension ainv(3*nsm,3*nsm),ds(nsm,ncm,2),coeff(nsm,ncm,4)
  dimension bt(3*nsm),dsi(nsm,ncm,2)

  logical neofrc


  save visc,aaa,ainv,indx

  integer nenergy
  real eps, mu

  dimension mu(3,3)
  logical nltest

  ! check the consistency of the nsm ncm parameters
  if ((nsc.ne.nsm).or.(ncc.ne.ncm)) call perr(1)

  ! define the proper thermodynamic forces for the calculation below
  do i = 1, ns ; do j = 1, nc(i) ; do k = 1, 2
    ds(i,j,k) = t(i)*dsi(i,j,k)/zsp(i,j)
  end do ; end do ; end do 

  ! check whether the matrices have to be recalculated
  if (.not.neofrc) then 

    ! calculate the viscosity
    do i = 1, ns ;  do j = 1, nc(i)

      do k = 1, 3 ; do l = 1, 3
        visc(k,l,i,j) = 0.
      end do ; end do 

      nltest = .false.
      call viscos(nsm, ncm, ns, i, j, nreg, bgradp, fc, fm, mmx, & 
              &    nenergy, eps, tau, m, t, xi, den, nc, mu)

      do k = 1, nleg ; do l = 1, nleg 
        visc(k,l,i,j) = mu(k,l)
      end do ; end do 

    end do ; end do 

    ! calculate the a matrices
    do i = 1, ns ; do j = 1, nc(i)

      do k = 1, 3 ; do l = 1, 3
        aaa(k,l,i,j) = 0.
      end do ; end do 
      do k = 1, 3
        aaa(k,k,i,j) = 1.
      end do 

      do k = 1, nleg ; do l = 1, nleg
        aaa(k,l,i,j) = visc(k,l,i,j) - xi(i,j)*la(k,l,i)
      end do ; end do 

      ! do the lu decomposition
      call ludcmp(aaa(1,1,i,j),3,3,indx,d)

      do k = 1, 3 ; do l = 1, 3
        akp(k,l) = 0.
      end do ; end do 
      do k = 1, 3
        akp(k,k) = 1.
      end do 

      do k = 1, 3
        call lubksb(aaa(1,1,i,j),3,3,indx,akp(1,k))
      end do 

      do k = 1, 3 ; do l = 1, 3
        aaa(k,l,i,j) = akp(k,l)
      end do ; end do 

    end do ; end do 
 
    ! now calculate the matrix for inversion
    do i = 1, 3*ns ;  do j = 1, 3*ns
      ainv(i,j) = 0.
    end do ; end do 
    do i = 1, 3*ns
      ainv(i,i) = 1.
    end do 

    do i = 1, ns

      do k = 1, 3 ; do l = 1, 3
        akp(k,l) = 0.
      end do ; end do 

      do k = 1, 3 ; do l = 1, 3 ;  do j = 1, nc(i)
         akp(k,l) = akp(k,l) + xi(i,j)**2*aaa(k,l,i,j)
      end do ; end do ; end do 

      do j = 1, ns ; do k = 1, nleg ; do l = 1, nleg ; do mt = 1, 3
        ainv((i-1)*3+k,(j-1)*3+l) =   ainv((i-1)*3+k,(j-1)*3+l) & 
                                  & - akp(k,mt)*lab(mt,l,i,j)
      end do ; end do ; end do ; end do 

    end do 

    ! calculate the lu decomposition
    call ludcmp(ainv,3*ns,3*nsm,indx,d)  

    ! end the if statement that determines whether the matrices are 
    ! recalculated
  endif

  ! now calculate the right hand side    
  do i = 1, ns ; do k = 1, 3 
    bt((i-1)*3+k) = 0.
    if (k.le.nleg) then 
      do j = 1, nc(i)
        bt((i-1)*3+k) = bt((i-1)*3+k)+1.6*xi(i,j)*zsp(i,j)* &
                      & den(i,j)*aaa(k,1,i,j)*epar
      end do 
      do l = 1, 3 ; do mt = 1, 2 ; do j = 1, nc(i)
         bt((i-1)*3+k) = bt((i-1)*3+k) - xi(i,j)*       &
            &   aaa(k,l,i,j)*visc(l,mt,i,j)*ds(i,j,mt)
       end do ; end do ; end do 
    end if 
  end do ; end do 


  call lubksb(ainv,3*ns,3*nsm,indx,bt)

  ! now calculate the indivudual particle velocities
  do i = 1, ns 

    do k = 1, 3 
      akp(k,1) = 0.
      if (k.le.nleg) then 
        do l = 1, 3 ; do j = 1, ns
          akp(k,1) = akp(k,1)+lab(k,l,i,j)*bt((j-1)*3+l)
        end do ; end do 
      endif 
    end do 
    do j = 1, nc(i) ;  do k = 1, 3
      upar(i,j,k) = 0.
      if (k.le.nleg) then 
        do l = 1, 3
          akp(l,2) = 0.
          do mt = 1, 2
            akp(l,2) = akp(l,2)+visc(l,mt,i,j)*ds(i,j,mt)
          end do
          upar(i,j,k) = upar(i,j,k) + aaa(k,l,i,j)*(   &
                      & xi(i,j)*akp(l,1)-akp(l,2))
        end do 
        upar(i,j,k) = upar(i,j,k) + aaa(k,1,i,j)*  & 
                    & 1.6E0*zsp(i,j)*den(i,j)*epar
      end if 
    end do ; end do 

  end do 


  ! finally calculate the coefficients coeff
  do i = 1, ns ; do k = 1, nc(i)

    coeff(i,k,3) = upar(i,k,1)
    do mt = 1, 3
      if (mt.ne.3) upar(i,k,mt) = upar(i,k,mt) + ds(i,k,mt)
      coeff(i,k,1) = coeff(i,k,1) + visc(1,mt,i,k)*upar(i,k,mt)
      coeff(i,k,2) = coeff(i,k,2) + visc(2,mt,i,k)*upar(i,k,mt)
    end do 
    coeff(i,k,4) = upar(i,k,1)

  end do ; end do 
 
return
end subroutine bp 
