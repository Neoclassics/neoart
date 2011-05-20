
      subroutine get_geom
     >     (shot, new_read, rho_in,
     >      r0_out, bp_invers_out, b_out, b_sq_out,
     >      b_sq_inv_out, r_bt_out, fc_out,
     >      n_f_out, 
     >      cos_b_sq_out, sin_b_sq_out, cos_b_lnB_out, sin_b_lnB_out)   
c
      implicit none

c
c      number of rho_poloidal values
c
       integer NR
       parameter (NR = 20)
c
c      number of rho_poloidal values up to separatrix
c
       integer NSEP
       parameter (NSEP = 13)
c
c      number of fourier coeficients
c
       integer NF
       parameter (NF = 3)
c
c      path to profile data
c
       character GRIDDIR*9
       parameter (GRIDDIR='./augint/')

c
c     input variables
c      
      integer
     >     shot,                !shot number
     >     new_read             !set 1 to force reading from file

      real
     >     rho_in               ! sqrt(V/2/pi**2/Rgeo) = rho_vol
      
c
c     output variables
c
      integer
     >     n_f_out              !number of Fourier coefficients
      
c     magnetic fields and averages
      
      real
     >     r0_out,              !R of magnetic axis 
     >     bp_invers_out,       !< 1/B_poloidal > 
     >     b_out,               !< total magnetic field>
     >     b_sq_out,            !< B**2>
     >     b_sq_inv_out,        !< 1./B**2>
     >     r_bt_out,            !< R B_T>
     >     fc_out,              !fraction of circulating particles
     >     cos_b_sq_out(NF),    !cosine fourier coefficient of B^2
     >     sin_b_sq_out(NF),    !sine fourier coefficient of B^2
     >     cos_b_lnB_out(NF),   !cosine fourier coefficient of BlnB
     >     sin_b_lnB_out(NF)    !sine fourier coefficient of BlnB     
      
c
c     local variables
c
      integer
     >     i,k,npoly,
     >     n_r, n_sep, n_f 

c     functions
      real
     >     interp

c     rho_pol ....
      real
     >     rgeo,                !major radius of geometrical axis
     >     r0,                  !r0 of magnetic axis
     >     time,  
     >     rho(NR),             !rho_poloidal grid
     >     rho_vol(NR),         !rho_volume grid
     >     rho_vol_sep,         !rho_volume at separatrix 
     >     rmaj_out(NR),        !major radius low field side
     >     rmaj_in(NR),         !major radius high field side
     >     grad(NR),            ! < gradient(rho_vol) >
     >     grad_sq(NR)          ! < gradient(rho_vol)**2 > 
      
c     magnetic fields and averages
      real
     >     bp_invers(NSEP),     !< 1/B_poloidal > 
     >     b(NSEP),             !< total magnetic field>
     >     b_sq(NSEP),          !< B**2>
     >     b_sq_inv(NSEP),      !< 1./B**2>
     >     r_bt(NSEP),          !< R B_T>
     >     fc(NSEP),            !fraction of circulating particles
     >     cos_b_sq(NSEP,NF),   !cosine fourier coefficient of B^2
     >     sin_b_sq(NSEP,NF),   !sine fourier coefficient of B^2
     >     cos_b_lnB(NSEP,NF),  !cosine fourier coefficient of BlnB
     >     sin_b_lnB(NSEP,NF),  !sine fourier coefficient of BlnB
     >     q(NSEP)              !safety factor       

c     file
      character*12 file
      integer index, finum, readfile, extfile

c     save variables
      save
     >     n_sep, n_f,
     >     r0, rho_vol,  
     >     bp_invers, b, b_sq, b_sq_inv, r_bt, fc,
     >     cos_b_sq, sin_b_sq, cos_b_lnB, sin_b_lnB
      
c*************************************************************      
c     Read file
c*************************************************************

      if (new_read.eq.1) then 
         index = 0
         extfile =0
         finum=10
         file = 'grid_00000.0'
         write(file(6:10),'(i5.5)') shot
         write(file(12:12),'(i1)')  index
         
         open(finum,file=GRIDDIR//file, status='old', err=9000)
         
         call advance(finum,'cv ',readfile,extfile)	
         read(readfile,*) rho_vol_sep, rgeo, r0, time
         
         call advance(finum,'cv ',readfile,extfile)
         read(readfile,*) n_r,n_sep,n_f

         if (n_r.gt.NR) goto 9010
         if (n_f.gt.NF) goto 9020
         
         call advance(finum,'cv ',readfile,extfile)
         read(readfile,*) (rho(i),i=1,n_r)
         
         call advance(finum,'cv ',readfile,extfile)
         read(readfile,*) (rho_vol(i),i=1,n_r)
         
         call advance(finum,'cv ',readfile,extfile)
         read(readfile,*) (rmaj_out(i),i=1,n_r)
         
         call advance(finum,'cv ',readfile,extfile)
         read(readfile,*) (rmaj_in(i),i=1,n_r)
         
         call advance(finum,'cv ',readfile,extfile)
         read(readfile,*) (grad(i),i=1,n_r)
         
         call advance(finum,'cv ',readfile,extfile)
         read(readfile,*) (grad_sq(i),i=1,n_r)
         
         call advance(finum,'cv ',readfile,extfile)
         read(readfile,*) (bp_invers(i),i=1,n_sep)
         
         call advance(finum,'cv ',readfile,extfile)
         read(readfile,*) (b(i),i=1,n_sep)
         
         call advance(finum,'cv ',readfile,extfile)
         read(readfile,*) (b_sq(i),i=1,n_sep)
         
         call advance(finum,'cv ',readfile,extfile)
         read(readfile,*) (b_sq_inv(i),i=1,n_sep)
         
         call advance(finum,'cv ',readfile,extfile)
         read(readfile,*) (r_bt(i),i=1,n_sep)
         
         call advance(finum,'cv ',readfile,extfile)
         read(readfile,*) (fc(i),i=1,n_sep)
      
         call advance(finum,'cv ',readfile,extfile)
         do k=1,n_f
            read(readfile,*) (cos_b_sq(i,k),i=1,n_sep)
         enddo
         
         call advance(finum,'cv ',readfile,extfile)
         do k=1,n_f
            read(readfile,*) (sin_b_sq(i,k),i=1,n_sep)
         enddo
         
         call advance(finum,'cv ',readfile,extfile)
         do k=1,n_f
            read(readfile,*) (cos_b_lnb(i,k),i=1,n_sep)
         enddo
         
         call advance(finum,'cv ',readfile,extfile)
         do k=1,n_f
            read(readfile,*) (sin_b_lnb(i,k),i=1,n_sep)
         enddo
         
         call advance(finum,'cv ',readfile,extfile)
         read(readfile,*) (q(i),i=1,n_sep)    
         
         close(finum)

c     rho_vol in meter

         do i=1, n_sep
            rho_vol(i) = rho_vol(i)*rho_vol_sep/100.
         enddo
      endif
c
c     intperpolate
c
      npoly=2
      
      bp_invers_out = interp(rho_vol,bp_invers,n_sep,npoly,rho_in)
      b_out = interp(rho_vol,b,n_sep,npoly,rho_in)
      b_sq_out = interp(rho_vol,b_sq,n_sep,npoly,rho_in)
      b_sq_inv_out = interp(rho_vol,b_sq_inv,n_sep,npoly,rho_in)
      r_bt_out = interp(rho_vol,r_bt,n_sep,npoly,rho_in)
      fc_out = interp(rho_vol,fc,n_sep,npoly,rho_in)
      do k=1,n_f
       cos_b_sq_out(k) = 
     +   interp(rho_vol,cos_b_sq(1,k),n_sep,npoly,rho_in)
      enddo
      do k=1,n_f
       sin_b_sq_out(k) = 
     +   interp(rho_vol,sin_b_sq(1,k),n_sep,npoly,rho_in)
      enddo
      do k=1,n_f
       cos_b_lnb_out(k) = 
     +   interp(rho_vol,cos_b_lnb(1,k),n_sep,npoly,rho_in)
      enddo      
      do k=1,n_f
       sin_b_lnb_out(k) = 
     +   interp(rho_vol,sin_b_lnb(1,k),n_sep,npoly,rho_in)
      enddo
      n_f_out = n_f
      r0_out = r0/100.
      
      return
c==================================================================
c
c     write errors and stop
c
9000  write(*,*) 'geometry file: '//GRIDDIR//file//' does not exist.'
      stop
9010  write(*,*) 'rho-grid in '//file//' to large !'
      stop
9020  write(*,*) 'To many fourier components in '//file//' !'
      stop        
      end
