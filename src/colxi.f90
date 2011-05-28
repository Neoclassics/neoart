!--------------------------------------------------------------------
!>   this subroutine calculates the collision frequencies and 
!>   the weighting factors xi
!>
!>   input    nsm    : maximum number of species
!>            ncm    : maximum number of charge states
!>            ns     : actual number of species
!>            nc     : array(nsm) gives the actual amount of 
!>                     charge states per species
!>            zsp    : array(nsm,ncm) the carge number
!>            den    : array(nsm,ncm) the density in 10^19 m^-3
!>            t      : array(nsm) temperature of species
!>            m      : array(nsm) mass of species
!>   output   tau    : array(nsm,nsm) collision frequency 
!>                     weighted over charge states.
!>            xi     : relative weight of every charge state.
!--------------------------------------------------------------------
subroutine colxi(nsm,ncm,ns,nc,zsp,den,t,m,tau,xi)

  implicit none

  integer nsm,ncm,i,j,ns,nc,k,l
  real den,t,m,dene,lnab,tau,xi,tcons,zsp
  dimension zsp(nsm,ncm),nc(nsm),den(nsm,ncm),t(nsm)
  dimension m(nsm),xi(nsm,ncm),tau(nsm,nsm)
 
  ! calculate the collision times
  tcons = 8.74202e6
  ! the collision time below is the double average (over charge 
  ! states) of the density devided by the collision time (equation 
  ! a3 and a4 of houlberg, phys plasmas 4 3230 (1997)) 
  ! first determine the electron density
  dene = 0.
  do i = 1, ns
    if (zsp(i,1).lt.0.) dene = den(i,1)
  end do 
  if (dene.eq.0.) then
    do i = 1, ns ; do j = 1, nc(i)
      dene = dene + den(i,j)*zsp(i,j)
    end do ; end do 
  endif
   
  do i = 1, ns ; do j = 1, ns

    ! calculate the coulomb logarithm 
    if ((zsp(i,1).lt.0.).and.(zsp(j,1).lt.0.)) then

      ! electron-electron
      lnab = 14.9 - 0.5*log(0.1*dene)+log(t(i))

    else

      if (zsp(i,1).lt.0.) then
        ! electron-ion
        lnab = 15.2-0.5*log(0.1*dene) + log(t(i))
      else
        if (zsp(j,1).lt.0.) then
          ! ion-electron
          lnab = 15.2-0.5*log(0.1*dene) + log(t(j))
        else
          ! ion-ion
          lnab = 17.3-0.5*log(0.1*dene)+1.5*log(0.5*(t(i)+t(j)))
        endif
      endif

    endif

    tau(i,j) = 0.
    do k = 1, nc(i) ; do l = 1, nc(j)
      tau(i,j) = tau(i,j) + den(i,k)*den(j,l)*zsp(i,k)**2*zsp(j,l)**2
    end do ; end do 

    tau(i,j) = tau(i,j)*tcons*lnab*sqrt(m(i)/t(i)**3)

  end do ; end do 
   
  ! calculate the xi coefficients
  do i = 1, ns
    xi(i,nc(i)) = 0.
    do j = 1, nc(i)
      xi(i,nc(i)) = xi(i,nc(i)) + den(i,j)*zsp(i,j)**2
    end do 
    do j = 1, nc(i)
      xi(i,j) = den(i,j)*zsp(i,j)**2 / xi(i,nc(i))
    end do 
  end do 

return
end subroutine colxi 
