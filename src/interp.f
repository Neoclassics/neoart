
      real function interp(x,y,npts,npoly,xin)
c
c     interpolation and extrapolation routine from BEVINGTON
c     for montonically increasing x-grid
c
c     x - independent variable of data array
c     y - dependent varaible of data array
c     npts - number of data points
c     npoly - degree of fitting polynomial (up to 9)
c     xin - input value of x
c

      implicit none

      integer npts,npoly

      real
     >x(npts),y(npts),
     >xin,denom,
     >delta(10),a(10),deltax,prod,sum

      integer
     >nterms,
     >i,k,j,i1,i2,ix,
     >imax,ixmax

      nterms = npoly+1
c
c     search for appropriate value of x(npts)
c

11    do 19 i=1,npts
      if (xin - x(i)) 13,17,19
13    i1 = i - nterms/2
      if (i1) 15,15,21
15    i1=1
      goto 21
17    interp=y(i)
18    goto 61
19    continue

      i1 = npts - nterms + 1
21    i2 = i1 + nterms -1
      if (npts - i2) 23,31,31
23    i2 = npts
      i1 = i2 - nterms + 1
25    if (i1) 26,26,31
26    i1 = 1
27    nterms = i2 - i1 + 1

c
c     evaluate deviations delta
c

31    denom = x(i1+1) - x(i1)
      deltax =  (xin - x(i1)) / denom
      do 35 i=1,nterms
      ix = i1 + i - 1
35    delta(i) = (x(ix) - x(i1)) / denom

c
c     accumulate coefficients
c

      a(1) =  y(i1) 
      do 50 k=2,nterms
      prod=1.
      sum=0.
      imax = k - 1
      ixmax = i1 + imax
      do 49 i=1,imax
      j = k - i
      prod = prod * (delta(k) - delta(j))
49    sum = sum - a(j)/prod
50    a(k) = sum + y(ixmax)/prod

c
c     accumulate sum of expansion
c

51    sum= a(1)
      do 57 j=2,nterms
      prod=1.
      imax = j-1
      do 56 i=1,imax
56    prod = prod*(deltax - delta(i))
57    sum = sum + a(j)*prod
60    interp = sum
61    return
      end      
