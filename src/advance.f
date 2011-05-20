
c
c=========================================
c
 
      subroutine advance(filenum,flag,readfile,extfile)
 
      implicit none
 
      character*3 flag,dummy
 
      integer filenum,readfile,extfile,counter

      dummy = '   '
      counter = 1
       
1     read(filenum,*,end=99) dummy
      if (dummy(1:2).ne.flag(1:2)) goto 1
      if (dummy(3:3).eq.' ') readfile=filenum
      if (dummy(3:3).eq.'#'.and.extfile.ne.0) readfile=extfile
      return

99    close(filenum) 
      write(*,*) 'End of file unit: ',filenum,' reached while reading!'
      stop
      end
