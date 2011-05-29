!--------------------------------------------------------------------
! this subroutine calculates all the quantities related with the
!>  geometry. several possibilities to supply input to this routine
!>  exist. one can either use analytic expressions for a circular
!>  equilibrium, read the coefficients from file (interface to 
!>  asdex upgrade, or call a routine that suplies the magnetic field
!>  information in hamada coordinates. these options are selected 
!>  through the parameter isel.
!>
!>  input   nmg    :  maximum number of quantities in fm used for
!>                    consistency check.
!>          isel   :  parameter that selects the various options
!>                    1 use hamada coordinates and call visgeom
!>                    2 use circular geometry
!>                    3 read the parameters from file
!>          ishot  :  shot number, only used if isel = 3
!>          rho    :  the surface label
!>          eps    :  accuracy required.
!>  output  bav    :  average of the magnetic field strength
!>          b2av   :  average of the magnetic field strength squared.
!>          bi2a   :  average of the inverse magnetic field strength
!>                    squared. 
!>          rbt    :  the product of major radius and toroidal magn.
!>                    field
!>          bgradp :  inner product of direction of magnetic field
!>                    and gradient of theta, a flux function 
!>          dpsidr :  gradient of poloidal flux towards dimensional 
!>                    radial coordinate
!>          rnq    :  field line length (approx r * q)
!>          fc     :  number of passing particles
!>          gclass :  geometry dependent coefficient for the cal-
!>                    culation of classical transport. gclass = 
!>                    < r^2 b_p^2 / b^2 >
!>          fm     :  array(nmaxgr) coefficients needed to 
!>                    calculate the viscosity in the pfirsch 
!>                    schlueter regime.
!>          mmx    :  the actual number of fourier coefficients
!>          r2i    :  flux surface average of r^-2
!>
!>  the routine call the folowing subroutines
!>
!>  visgeom :  if isel =1 the geometry should be suplied in hamada 
!>             coordinates
!>  circgeom:  if isel =2 the parameters of circular geometry are 
!>             suplied through this routine
!--------------------------------------------------------------------
subroutine geom(nmg,isel,ishot,rho,eps,bav,b2av,bi2a,rbt,bgradp, &
                 & dpsidr,rnq,fc,gclass,fm,mmx,r2i)

  implicit none

  real rho, eps
  integer nmaxgr, isel, nmg

  parameter(nmaxgr = 1000)

  real bmax, tmax, dthet, thet, bgradp, twopi, bdum, bav
  real b2av, err, error, lam, dlam, fdenom, fc, fco, pi
  real rn, e, bn, q, dum, bi2a, rbt, dpsidr, rnq, gclass
  real r2i
  real b(nmaxgr), theta(nmaxgr),sin1(nmaxgr),sin2(nmaxgr)
  real cos1(nmaxgr), cos2(nmaxgr), fm(nmaxgr), a(5)

  integer ngr,mmx,i,j,nc,new_read,ishot
  real raxis, rgeo, idlp, cosb2, sinb2, cosbl, sinbl

  dimension cosb2(3), sinb2(3), cosbl(3), sinbl(3)
 
  data new_read / 1 / 

  !>  the constant pi, 2 pi
  pi = 4.*atan(1.)
  twopi = 2.*pi

  if (nmg.ne.nmaxgr) call perr(3)
  
  ! The different options 
  select case(isel) 
     
  ! in this part the coefficients are calculated through calls
  ! to the routine visgeom which should give the magnetic field
  ! strength in hamada coordinates.
  !
  ! ***** warning: r2i is not calculated in this part therefore 
  !            no ware pinch can be caculated but the rest 
  !            is o.k.
  case(1) 

    ! first calculate the coefficients f_m, the number of passing
    ! particles f_c, and the quantity (n grad theta) in the coordinates
    ! of shaing.
      
    ! the initial choise of the number of grid points in the theta 
    ! direction and the number of terms in the fourier expansion in 
    ! the poloidal angle
    ngr = 200
    mmx = 20

    ! make sure the calculation is done once 
    error = 2.*eps
 
    ! repeat until accurate enough 
    do while (error.gt.eps)  

      ! distance between the grid in theta direction
      dthet = 1 / real(ngr)

      ! calculate the magnetic field on the grid points and determine 
      ! the maximum of the field. (the minimum is assumed to be at 
      ! thet = 0)
      bmax = 0.
      tmax = 0.
      bav  = 0.
      b2av = 0.
      bi2a = 0.
      do i = 1, ngr
        thet = (i-1.)*dthet
        call visgeom(rho,thet,bgradp,rbt,b(i),dpsidr)
        bav  = bav  + b(i)
        b2av = b2av + b(i)**2
        bi2a = bi2a + 1/b(i)**2
        if (b(i).gt.bmax) then
          bmax = b(i)
          tmax = thet
        endif
      end do 
      bav  = bav  * dthet
      b2av = b2av * dthet
      bi2a = bi2a * dthet
      ! estimate the error in the integrals
      error = 0. 
      err = 0.
      do i = 2, ngr-1
        err = max(err, abs( (b(i-1) - 2.*b(i) + b(i+1))/(dthet**2) ))
      end do 
      error = max(dthet**2*err/(abs(bav)), error)

      ngr = ngr * 2
      if (ngr.gt.nmaxgr) call perr(5)

    end do 

    ! calculate the poloidal angle in the coordinates of shaing
    ! (this coordinate is denoted theta instead of thet for the 
    ! hamada coordinates
    theta(1) = 0.
    do i = 2, ngr
      thet = (i-1.5)*dthet
      call visgeom(rho,thet,bgradp,rbt,bdum,dpsidr)
      theta(i) = theta(i-1) + bdum*dthet
    end do 
    do i = 1, ngr
      theta(i) = theta(i) * twopi / bav
    end do 

    ! calculate the coefficients fm
    ! first do the integrals
    do i = 1, mmx
      sin1(i) = 0.
      sin2(i) = 0.
      cos1(i) = 0.
      cos2(i) = 0.
      do j = 1, ngr
        sin1(i) = sin1(i) + cos(i*theta(j))*b(j)*log(b(j))
        sin2(i) = sin2(i) + cos(i*theta(j))*b(j)**2
        cos1(i) = cos1(i) + sin(i*theta(j))*b(j)*log(b(j))
        cos2(i) = cos2(i) + sin(i*theta(j))*b(j)**2
      end do 
      fm(i) = 2. * ( i * twopi * bgradp * dthet )**2 /     &
            & (bav**3 * b2av) * (sin1(i)*sin2(i)+cos1(i)*cos2(i))
    end do 

    bgradp = twopi * bgradp / bav
      
    rnq = 1./bgradp 

    ! now determine the number of passing particles
    nc = 100
    fco = 0.
    fc  = 0. 
    do while ((abs(fco-fc).gt.eps*abs(fc)).or.(fco.eq.0.)) 

      fco = fc 
      if (nc.gt.nmaxgr) stop 'err 2 no convergence in viscos'
      dlam = 1./ (bmax*real(nc))
      fc  = 0.
      do i = 1, nc
        lam = (i-0.5)*dlam
        fdenom = 0.
        do j = 1, ngr
          fdenom = fdenom + sqrt(1.- lam*b(j))*dthet
        end do 
        fc = fc + lam * dlam / fdenom
      end do 
      nc = nc*2 

    end do 
    ! the number of trapped and passing particles is  
    fc = 0.75*b2av*fc 

    ! the quantity gclass is not yet programmed for hamada 
    ! coordinates
    gclass = 0.

    return 

  ! in this part analytic expressions for an circular geometry 
  ! are used. 
  !10000 continue
  case(2) 

    ! obtain the values of the major radius inverse aspect ratio
    ! safety factor and magnetic field strength.
    call circgeom(2,rho,rn,e,q,bn)
 
    ! for circular surface the flux label r is used. 
    dpsidr = rn * bn / sqrt(1 + q**2 * (1/e**2 - 1.))

    ! averages of the magnetic field strengths
    bav = bn
    b2av = bn**2/sqrt(1-e**2)
    r2i  = 1./(rn**2*sqrt(1.-e**2))
    bi2a = (1 + 1.5 * e**2 ) / bn**2

    ! the product of major radius and toroidal magnetic field
    rbt = sqrt(1 - e**2) * rn * bn / sqrt(1 + e**2 *(1/q**2 - 1))

    ! measure of the field line length
    rnq = rn*q

    bgradp = 1/(rn*sqrt(e**2 + q**2*(1-e**2)))

    a(1) = -1.46655
    a(2) = 1.0241
    a(3) = -1.20107
    a(4) = 1.356234
    a(5) = -0.662881

    dum =  1 - e
    fc = dum
    do i = 1, 5
      dum = dum*sqrt(e)
      fc = fc + a(i)*dum
    end do 

    i = 0
    do while ((abs(dum).gt.0.01*eps).or.(i.lt.15)) 
      i = i + 1  
      if (i.gt.nmaxgr) call perr(5)
      fm(i) = 2*i*bgradp**2*((-1.+sqrt(1-e**2))/e)**(2*i)
      dum = (2 - e**2 - 2*sqrt(1-e**2))/e**2
      dum = i*dum**(i-1) + (1-i)*dum**i
    end do 
    mmx = i

    ! the geometry dependent quantity for the calculation 
    ! of classical transport
    gclass = (rn*e)**2*(1 + 3.*e**2/2.)/(q**2 + e**2 * (1 - q**2))      

    return

  ! in this part the coefficients are read from file
  !20000 continue
  case(3) 
      
    call get_geom(ishot, new_read, rho, raxis, idlp, bav, b2av, &
            &      bi2a, rbt, fc, mmx, cosb2, sinb2, cosbl, sinbl,r2i)   

    if (mmx.ne.3) stop 'inconsistent settings with nf.ne.3'

    ! warning rgeo is set to a fixed value
    rgeo = 1.65
      
    dpsidr = twopi*rho*rgeo / idlp
    bgradp = twopi / (bav*idlp)

    do i = 1, mmx
      fm(i) = 2.*(i*bgradp)**2/(b2av*bav)*(cosb2(i)*cosbl(i)+  &
            &     sinb2(i)*sinbl(i))
    end do 

    rnq = 1. / bgradp
  
    if (new_read.eq.1) new_read = 0
     
  case default 
    return 

  end select

return
end subroutine geom 
