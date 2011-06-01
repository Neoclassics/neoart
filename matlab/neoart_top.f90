      subroutine mexfunction(nlhs, plhs, nrhs, prhs)
      implicit none

!   the mex function should have the name of the file. In this case neoart
!   Then typically in matlab one calls:
!   >> [a,b,c,..]=neoart(arg1,arg2,...)
!
!   nlhs, nrhs: number of input arguments on left and right-hand side respectively
!   thus one can test and depending on number of arguments choose different options
!
!   plhs, prhs: pointer to the different arguments: plhs(i), i=1,nlhs
!   (plhs and prhs are integer arrays)
!
!   fmex function: mxGetPr gets the pointer value in prhs, plhs
!   .              mxCopyPtrToReal8: copy pointed values to local array
!   .              mxCopyReal8ToPtr: copy local array values to pointed array
!   .              mxGetM: get nb of rows in pointed array
!   .              mxGetN: get nb of columns in pointed array
!   .              mxCreateDoubleMatrix: creates a matrix for the return argument
!
!   In this example, matlab call expected:
!
!   >> [c1,c2,c3,c4] = neoart(ns,nc,M,T,zsp,den,ds,isel,rho{,RN,BN,q},ic,eparr{,eps,nreg,sigma,nleg,nenergy,ncof,neogeo,neofrc})



! variables for the mexfunction I/O pointers arrays
      integer*8 plhs(*), prhs(*)
      integer nlhs, nrhs
      integer*8 mxGetPr
      integer*8 mxGetM, mxGetN, mxCreateDoubleMatrix, mxCreateNumericArray, mxDOUBLE_CLASS
! Pointers and variables for matlab input arguments
      integer*8 pns, pnc, pM, pT, pzsp, pden, pds, pisel, prho, pRN, pBN, pq, &
	 &  pic, peparr, peps, pnreg, psigma, pnleg, pnenergy, pncof, pneogeo, pneofrc
!
      integer neogeo, neofrc
      integer ns, isel, ic, nreg, nleg, nenergy, ncof
      REAL*8 zns, zisel, zic, znreg, znleg, znenergy, zncof, zneogeo, zneofrc
      integer, pointer :: nc(:)
      REAL*8 :: rho, RN, BN, q, eparr, eps, sigma(4)
      REAL*8,pointer :: znc(:), M(:), T(:)
      REAL*8,pointer :: zsp(:,:), den(:,:), ds(:,:,:)
      REAL*8,pointer :: zsp_tmp(:,:), den_tmp(:,:), ds_tmp(:,:,:)
! Pointers and variables for matlab output arguments
      integer*8 pcoeff_out1, pcoeff_out2, pcoeff_out3, pcoeff_out4
      REAL*8, pointer :: coeff(:,:,:)
      REAL*8, pointer :: c1(:,:), c2(:,:), c3(:,:), c4(:,:)
! Additional variables
      integer ishot, nzm, ninrow, nincol, nin, ipr, i, j, k
      real*8 :: e
      integer*8 ninint8, nelint8, ninrowint8, nincolint8
      integer*8 noutint8, noutrowint8, noutcolint8
! dimension for the arrays used in the NEOART call (arrays must have the size set in elem_config.inc)
      include '../src/elem_config90.inc'
      integer NSM, NCM 	
      parameter(NSM = NELMAX+2)
      parameter(NCM = NIONMAX)

!
! 1. Default values
      ishot = 0
      eps = 1.E-5
	  nreg = 0
	  sigma(1) = 1
	  sigma(2) = 2
	  sigma(3) = 3
	  sigma(4) = 4
	  nleg = 3
	  nenergy = 1
	  ncof = 1
	  neogeo = 1
	  neofrc = 0
!
! 2. Checks
      if (nrhs .ne. 22 .and. nrhs .ne. 19 .and. nrhs .ne. 14 .and. nrhs .ne. 11) then
          print *
          print *,' Number of inputs must be 22, 19, 14 or 11'
          print *,' Type ''help neoart'' for details'
          return
      end if
!      
! 3. Read inputs
! 3.1 number of species: ns
      pns = mxGetPr(prhs(1))
!   call mxCopyPtrToInteger4(pns, ns, 1)  ! does not work => use real8
      nelint8 = int(1,8)
      call mxCopyPtrToReal8(pns, zns, nelint8)
	  ns = zns
	  if (ns .gt. NSM) then
          print *
          print *,' Error: the maximum number of species in NEOART is ', NSM
          print *,' Can be changed in the file elem_config.inc (need to recompile NEOART)'
          return
      end if	  
! 3.2 charge state of each species	  
      pnc = mxGetPr(prhs(2))
      ninrowint8 = mxGetM(prhs(2))
      nincolint8 = mxGetN(prhs(2))
      ninrow = int(ninrowint8)
      nincol = int(nincolint8)
	  nin = max(ninrow,nincol)
	  if (nin .ne. ns) then
          print *
          print *,' Error: nc must have ns elements'
          print *
          return
      end if
	  allocate(nc(NSM))
 	  allocate(znc(nin))
      ninint8 = int(nin,8)
      call mxCopyPtrToReal8(pnc, znc, ninint8)
	  nc = znc
	  nzm = maxval(znc(:))
	  if (nzm .gt. NCM) then
          print *
          print *,' Error: the max. number of charge states in NEOART is ', NCM
          print *,' Can be changed in the file elem_config.inc (need to recompile NEOART)'
          return
      end if	  
! 3.3 species mass: M
      pM = mxGetPr(prhs(3))
      ninrowint8 = mxGetM(prhs(3))
      nincolint8 = mxGetN(prhs(3))
      ninrow = int(ninrowint8)
      nincol = int(nincolint8)
	  nin = max(ninrow,nincol)
	  if (nin .ne. ns) then
          print *
          print *,' Error: M must have ns elements'
          print *
          return
      end if
	  allocate(M(NSM))
      ninint8 = int(nin,8)
      call mxCopyPtrToReal8(pM, M, ninint8)
! 3.4 species temperature: T
      pT = mxGetPr(prhs(4))
      ninrowint8 = mxGetM(prhs(4))
      nincolint8 = mxGetN(prhs(4))
      ninrow = int(ninrowint8)
      nincol = int(nincolint8)
	  nin = max(ninrow,nincol)
	  if (nin .ne. ns) then
          print *
          print *,' Error: T must have ns elements'
          print *
          return
      end if
	  allocate(T(NSM))
      ninint8 = int(nin,8)
      call mxCopyPtrToReal8(pT, T, ninint8)
! 3.5 charge states: zsp
      pzsp = mxGetPr(prhs(5))
      ninrowint8 = mxGetM(prhs(5))
      nincolint8 = mxGetN(prhs(5))
      ninrow = int(ninrowint8)
      nincol = int(nincolint8)
	  nin = ninrow*nincol
	  if (ninrow .ne. ns .or. nincol .ne. nzm) then
          print *
          print *,' Error: zsp must be of size (ns,max(nc))'
          print *
          return
      end if
	  allocate(zsp_tmp(ninrow,nincol))
	  allocate(zsp(NSM,NCM))
      ninint8 = int(nin,8)
      call mxCopyPtrToReal8(pzsp, zsp_tmp, ninint8)
	  do i = 1, ninrow
	   do j = 1, nincol
          zsp(i,j)=zsp_tmp(i,j)		
       end do
	  end do	  
! 3.6 densities: den
      pden = mxGetPr(prhs(6))
      ninrowint8 = mxGetM(prhs(6))
      nincolint8 = mxGetN(prhs(6))
      ninrow = int(ninrowint8)
      nincol = int(nincolint8)
	  nin = ninrow*nincol
	  if (ninrow .ne. ns .or. nincol .ne. nzm) then
          print *
          print *,' Error: den must be of size (ns,max(nc))'
          print *
          return
      end if
	  allocate(den_tmp(ninrow,nincol))
	  allocate(den(NSM,NCM))
      ninint8 = int(nin,8)
      call mxCopyPtrToReal8(pden, den_tmp, ninint8)
	  do i = 1, ninrow
	   do j = 1, nincol
          den(i,j)=den_tmp(i,j)		
       end do
	  end do	  
! 3.7 thermodynamic forces: ds
      pds = mxGetPr(prhs(7))
      ninrowint8 = mxGetM(prhs(7))
      nincolint8 = mxGetN(prhs(7))
      ninrow = int(ninrowint8)
      nincol = int(nincolint8)
	  nin = ninrow*nincol
	  if (ninrow .ne. ns .or. nincol .ne. 2*nzm) then
          print *
          print *,' Error: ds must be of size (ns,max(nc),2)'
          print *
          return
      end if
	  allocate(ds_tmp(ninrow,nincol/2,2))
      allocate(ds(NSM,NCM,2))
	  ninint8 = int(nin,8)
      call mxCopyPtrToReal8(pds, ds_tmp, ninint8)
	  do i = 1, ninrow
	   do j = 1, nincol/2
	    do k = 1, 2
          ds(i,j,k)=ds_tmp(i,j,k)		
		end do
       end do
	  end do
! 3.8 geometry selection: isel
      pisel = mxGetPr(prhs(8))
      nelint8 = int(1,8)
      call mxCopyPtrToReal8(pisel, zisel, nelint8)
	  isel = zisel
! 3.9 flux surface label: rho
      prho = mxGetPr(prhs(9))
      nelint8 = int(1,8)
      call mxCopyPtrToReal8(prho, rho, nelint8)
! 3.9.1 to 3.9.3: reads the values of RN, BN and q if needed
      select case(isel)
	    case(2)
		  if (nrhs .ne. 22 .and. nrhs .ne. 14) then
            print *
            print *,'  Error: isel=2, number of inputs must be 22, or 14'
            print *
          return
          end if
		  ipr = 1
          pRN = mxGetPr(prhs(9+ipr))
          nelint8 = int(1,8)
          call mxCopyPtrToReal8(pRN, RN, nelint8)
		  ipr=ipr+1
          pBN = mxGetPr(prhs(9+ipr))
          nelint8 = int(1,8)
          call mxCopyPtrToReal8(pBN, BN, nelint8)
		  ipr=ipr+1
          pq = mxGetPr(prhs(9+ipr))
          nelint8 = int(1,8)
          call mxCopyPtrToReal8(pq, q, nelint8)
		  e = rho/RN	  
		case(4)
		  if (nrhs .ne. 19 .and. nrhs .ne. 11) then
            print *
            print *,'  Error: isel=4, number of inputs must be 19, or 11'
            print *
          return
          end if
		  ipr = 0
		case default
          print *
          print *,' Error: isel=', isel, ' only values allowed are 2 and 4'
          print *
          return
	  end select
! 3.10 type of flux contribution: ic
      pic = mxGetPr(prhs(10+ipr))
      nelint8 = int(1,8)
      call mxCopyPtrToReal8(pic, zic, nelint8)
	  ic = zic
! 3.11 parallel electric field: eparr
      peparr = mxGetPr(prhs(11+ipr))
      nelint8 = int(1,8)
      call mxCopyPtrToReal8(peparr, eparr, nelint8)
! optional
      if (nrhs .gt. 14) then
! 3.12 accuracy: eps
      peps = mxGetPr(prhs(12+ipr))
      nelint8 = int(1,8)
      call mxCopyPtrToReal8(peps, eps, nelint8)
! 3.13 regime for viscosity: nreg
      pnreg = mxGetPr(prhs(13+ipr))
      nelint8 = int(1,8)
      call mxCopyPtrToReal8(pnreg, znreg, nelint8)
	  nreg = znreg
! 3.14 coupling: sigma	  
      psigma = mxGetPr(prhs(14+ipr))
      ninrowint8 = mxGetM(prhs(14+ipr))
      nincolint8 = mxGetN(prhs(14+ipr))
      ninrow = int(ninrowint8)
      nincol = int(nincolint8)
	  nin = max(ninrow,nincol)
	  if (nin .ne. 4) then
          print *
          print *,' Error: sigma must have 4 elements'
          print *
          return
      end if
      ninint8 = int(nin,8)
      call mxCopyPtrToReal8(psigma, sigma, ninint8)
! 3.15 number of Legendre polynomial: nleg
      pnleg = mxGetPr(prhs(15+ipr))
      nelint8 = int(1,8)
      call mxCopyPtrToReal8(pnleg, znleg, nelint8)
	  nleg = znleg
! 3.16 energy selection (viscosity): nenergy
      pnenergy = mxGetPr(prhs(16+ipr))
      nelint8 = int(1,8)
      call mxCopyPtrToReal8(pnenergy, znenergy, nelint8)
	  nenergy = znenergy
! 3.17 ion-electron collisions: ncof
      pncof = mxGetPr(prhs(17+ipr))
      nelint8 = int(1,8)
      call mxCopyPtrToReal8(pncof, zncof, nelint8)
	  ncof = zncof
! 3.18 recalculate geometry: neogeo
      pneogeo = mxGetPr(prhs(18+ipr))
      nelint8 = int(1,8)
      call mxCopyPtrToReal8(pneogeo, zneogeo, nelint8)
	  neogeo = zneogeo
! 3.19 recalculate viscosity and friction matrixes: neofrc
      pneofrc = mxGetPr(prhs(19+ipr))
      nelint8 = int(1,8)
      call mxCopyPtrToReal8(pneofrc, zneofrc, nelint8)
	  neofrc = zneofrc
	  end if
!
! 4. NEOART
      allocate(coeff(NSM,NCM,4))
! 4.1 Store the geometry values if isel=2
      if (isel .eq. 2) then
        CALL CIRCGEOM(1,RHO,RN,E,Q,BN)
      end if
! 4.2 Call the main routine	  
      CALL NEOART(ns,nc,NSM,NCM,zsp,M,T,den,ds,rho,eps, &
     &            isel,ishot,nreg,sigma,nleg,nenergy,ncof, &
     &            neogeo,neofrc,ic,eparr,coeff)
!
! 5. Export outputs
	  allocate(c1(ns,nzm))
	  allocate(c2(ns,nzm))
	  allocate(c3(ns,nzm))
	  allocate(c4(ns,nzm))
	  do i = 1, ns
	   do j = 1, nzm
          c1(i,j)=-coeff(i,j,1)		
          c2(i,j)=-coeff(i,j,2)		
          c3(i,j)=coeff(i,j,3)		
          c4(i,j)=coeff(i,j,4)		
       end do
	  end do
      if (nlhs .ge. 4) then
        noutrowint8=int(ns,8)
		noutcolint8=int(nzm,8)
        nelint8=int(0,8)
!        plhs(1) = mxCreateFull(noutrowint8,noutcolint8,nelint8)
        plhs(1) = mxCreateDoubleMatrix(noutrowint8,noutcolint8,0)
        pcoeff_out1 = mxGetPr(plhs(1))
!        plhs(2) = mxCreateFull(noutrowint8,noutcolint8,nelint8)
        plhs(2) = mxCreateDoubleMatrix(noutrowint8,noutcolint8,0)
        pcoeff_out2 = mxGetPr(plhs(2))
!        plhs(3) = mxCreateFull(noutrowint8,noutcolint8,nelint8)
        plhs(3) = mxCreateDoubleMatrix(noutrowint8,noutcolint8,0)
        pcoeff_out3 = mxGetPr(plhs(3))
!        plhs(4) = mxCreateFull(noutrowint8,noutcolint8,nelint8)
        plhs(4) = mxCreateDoubleMatrix(noutrowint8,noutcolint8,0)
        pcoeff_out4 = mxGetPr(plhs(4))
		noutint8 = int(ns*nzm,8)
        call mxCopyReal8ToPtr(c1(:,:), pcoeff_out1, noutint8)
        call mxCopyReal8ToPtr(c2(:,:), pcoeff_out2, noutint8)
        call mxCopyReal8ToPtr(c3(:,:), pcoeff_out3, noutint8)
        call mxCopyReal8ToPtr(c4(:,:), pcoeff_out4, noutint8)
      end if
!
! 6. Deallocate
      deallocate(nc)
      deallocate(znc)
      deallocate(M)
      deallocate(T)
      deallocate(zsp)
      deallocate(zsp_tmp)
      deallocate(den)
      deallocate(den_tmp)
      deallocate(ds)
      deallocate(ds_tmp)
      deallocate(coeff)
      return
	  end
	  
