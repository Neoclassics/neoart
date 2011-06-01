
      subroutine get_geom_ch
     >     (new_read, rho_in, nmaxgr,
     >      bav_out, b2av_out, bi2a_out, rbt_out, bgradp_out,
     >      dpsidr_out, rnq_out, fc_out, gclass_out,
     >      fm_out, mmx_out, r2i_out)
c
      implicit none
c
c     input variables
c      
      integer
     >     new_read,             !set 1 to force reading from file
     >     nmaxgr
c
      real
     >     rho_in               ! radius of interest : rho_in=(Rmax-Rmin)/2/R0EXP (same def. as in GKW)
c
c     output variables
c
      integer
     >     mmx_out              !number of Fourier coefficients
      
c     magnetic fields and averages
      
      real
     >     bav_out,             !<B> 
     >     b2av_out,            !<B^2> 
     >     bi2a_out,            !<1/B^2>
     >     rbt_out,             !< R B_T>
     >     bgradp_out,
     >     dpsidr_out,
     >     rnq_out,
     >     fc_out,              !fraction of circulating particles
     >     gclass_out,
     >     r2i_out
      real fm_out(NMAXGR)
c
c     local variables
c
c
      integer npsi ! number of points for the radial grid
      integer mmx_l ! number of fourier coefficients
      real  r0exp  ! normalising radius used in CHEASE
      character (len=20) tdum
      real  dum
      integer
     >     i,k,npoly, ierr
c     functions
      real
     >     interp
c     rho
      real, allocatable ::
     >     rho_l(:)              !rho=(Rmax-Rmin)/2/R0EXP (same def. as in GKW)
c     magnetic fields and averages
      real, allocatable :: 
     >     bav_l(:),
     >     b2av_l(:),
     >     bi2a_l(:),
     >     rbt_l(:),
     >     bgradp_l(:),
     >     dpsidr_l(:),
     >     rnq_l(:),
     >     fc_l(:),
     >     gclass_l(:),
     >     fm_l(:,:),
     >     r2i_l(:)

c     file
      character*2 GRIDDIR
      parameter (GRIDDIR='./')
      character*15 file
      integer finum, readfile, extfile

c     save variables
      save
     >     npsi, mmx_l,
     >     rho_l, 
     >     bav_l, b2av_l, bi2a_l, rbt_l, bgradp_l, dpsidr_l, rnq_l,
     >     fc_l, gclass_l, fm_l, r2i_l
      
c*************************************************************      
c     Read file
c*************************************************************

      if (new_read.eq.1) then 
         extfile =0
         finum=10
         file = 'neoart_geom.dat'
         
         open(finum,file=GRIDDIR//file, status='old', err=9000)
         rewind finum
         
         ! grid size
         ! number of points for psi-grid and fourier coefficients
         read(finum,*) tdum, npsi, tdum, mmx_l
         mmx_out=mmx_l
         ! allocate arrays
         allocate(rho_l(1:npsi),stat=ierr)
         if (ierr.ne.0) then
          stop 'Could not allocate rho_l in get_geom_ch'
         endif
         allocate(bav_l(1:npsi),stat=ierr)
         if (ierr.ne.0) then
          stop 'Could not allocate bav_l in get_geom_ch'
         endif
         allocate(b2av_l(1:npsi),stat=ierr)
         if (ierr.ne.0) then
          stop 'Could not allocate b2av_l in get_geom_ch'
         endif
         allocate(bi2a_l(1:npsi),stat=ierr)
         if (ierr.ne.0) then
          stop 'Could not allocate bi2a_l in get_geom_ch'
         endif
         allocate(rbt_l(1:npsi),stat=ierr)
         if (ierr.ne.0) then
          stop 'Could not allocate rbt_l in get_geom_ch'
         endif
         allocate(bgradp_l(1:npsi),stat=ierr)
         if (ierr.ne.0) then
          stop 'Could not allocate bgradp_l in get_geom_ch'
         endif
         allocate(dpsidr_l(1:npsi),stat=ierr)
         if (ierr.ne.0) then
          stop 'Could not allocate dpsidr_l in get_geom_ch'
         endif
         allocate(rnq_l(1:npsi),stat=ierr)
         if (ierr.ne.0) then
          stop 'Could not allocate rnq_l in get_geom_ch'
         endif
         allocate(fc_l(1:npsi),stat=ierr)
         if (ierr.ne.0) then
          stop 'Could not allocate fc_l in get_geom_ch'
         endif
         allocate(gclass_l(1:npsi),stat=ierr)
         if (ierr.ne.0) then
          stop 'Could not allocate gclass_l in get_geom_ch'
         endif
         allocate(fm_l(1:npsi,1:mmx_l),stat=ierr)
         if (ierr.ne.0) then
          stop 'Could not allocate fm_l in get_geom_ch'
         endif
         allocate(r2i_l(1:npsi),stat=ierr)
         if (ierr.ne.0) then
          stop 'Could not allocate r2i_l in get_geom_ch'
         endif
         
         ! read data
         ! normalising radius
         read(finum,*) tdum, r0exp
         ! not used
         read(finum,*) tdum, dum, tdum, dum
         read(finum,*) tdum, (dum, i=1,npsi)
         read(finum,*) tdum, (dum, i=1,npsi)
         ! amin -> rho
         read(finum,*) tdum, (rho_l(i), i=1,npsi)
         do i=1, npsi
            rho_l(i) = rho_l(i)/r0exp
         enddo
         ! damindpsi -> dpsidr
         read(finum,*) tdum, (dpsidr_l(i), i=1,npsi)
         do i=1, npsi
            dpsidr_l(i) = r0exp/dpsidr_l(i)
         enddo
         ! not used
         read(finum,*) tdum, (dum, i=1,npsi)
         read(finum,*) tdum, (dum, i=1,npsi)
         read(finum,*) tdum, (dum, i=1,npsi)
         read(finum,*) tdum, (dum, i=1,npsi)
         read(finum,*) tdum, (dum, i=1,npsi)
         read(finum,*) tdum, (dum, i=1,npsi)
         ! rbt
         read(finum,*) tdum, (rbt_l(i), i=1,npsi)
         ! not used
         read(finum,*) tdum, (dum, i=1,npsi)
         ! bav
         read(finum,*) tdum, (bav_l(i), i=1,npsi)
         ! b2av
         read(finum,*) tdum, (b2av_l(i), i=1,npsi)
         ! bi2a
         read(finum,*) tdum, (bi2a_l(i), i=1,npsi)
         ! not used
         read(finum,*) tdum, (dum, i=1,npsi)
         ! r2iav
         read(finum,*) tdum, (r2i_l(i), i=1,npsi)
         ! gclass
         read(finum,*) tdum, (gclass_l(i), i=1,npsi)
         ! bgradp
         read(finum,*) tdum, (bgradp_l(i), i=1,npsi)
         ! fc
         read(finum,*) tdum, (fc_l(i), i=1,npsi)
         ! fm
         do k=1,mmx_l
           read(finum,*) tdum, (fm_l(i,k), i=1,npsi)
         enddo
         
         close(finum)
      endif
c
      do i=1, npsi
         rnq_l(i) = 1./bgradp_l(i)
      enddo
c
      if (rho_in.gt.rho_l(npsi)) goto 9010
c
c     intperpolate
c
      npoly=2
      
      bav_out = interp(rho_l,bav_l,npsi,npoly,rho_in)
      b2av_out = interp(rho_l,b2av_l,npsi,npoly,rho_in)
      bi2a_out = interp(rho_l,bi2a_l,npsi,npoly,rho_in)
      rbt_out = interp(rho_l,rbt_l,npsi,npoly,rho_in)
      bgradp_out = interp(rho_l,bgradp_l,npsi,npoly,rho_in)
      dpsidr_out = interp(rho_l,dpsidr_l,npsi,npoly,rho_in)
      rnq_out = interp(rho_l,rnq_l,npsi,npoly,rho_in)
      fc_out = interp(rho_l,fc_l,npsi,npoly,rho_in)
      gclass_out = interp(rho_l,gclass_l,npsi,npoly,rho_in)
      r2i_out = interp(rho_l,r2i_l,npsi,npoly,rho_in)
      do k=1,mmx_l
       fm_out(k) = interp(rho_l,fm_l(1:npsi,k),npsi,npoly,rho_in)
      enddo
      return
c==================================================================
c
c     write errors and stop
c
9000  write(*,*) 'geometry file: '//GRIDDIR//file//' does not exist.'
      stop
9010  write(*,*) 'rho-grid in '//file//' too large !'
      stop
      end
