   
      function erf(xxx)
c------------------------------------------------------
c     this routine computes the error function.
c------------------------------------------------------
      IMPLICIT NONE

      REAL ERF, XXX, A, PI, TUSQPI, SIGN, XCG, X2, 
     +       SUM, TERM, T1, T2, P, ETA, PHIP
      INTEGER K, KMAX

      dimension a(5)
      pi = 4.*atan(1.)
      tusqpi = 2 / sqrt(pi)
      sign=1.
      if (xxx .lt. 0.) sign=-1.
      xcg=sign*xxx
      x2=xcg*xcg
      if (xcg .ge. .6) go to 1
      sum=xcg
      term=xcg
      kmax=6
      do 2 k=1,kmax
      t1=float(k)
      t2=float(2*k+1)/float(2*k-1)
      term=-term*x2/(t1*t2)
      sum=sum+term
    2 continue
      erf=tusqpi*sum
      erf=sign*erf
      return
    1 continue
      p=.3275911
      a(1)=.225836846
      a(2)=-.252128668
      a(3)=1.25969513
      a(4)=-1.287822453
      a(5)=.94064607
      eta=1./(1.+p*xcg)
      phip=tusqpi*exp(-x2)
      term=(((a(5)*eta+a(4))*eta+a(3))*eta+a(2))*eta+a(1)
      erf=1.-term*eta*phip
      erf=sign*erf
      return
      end
