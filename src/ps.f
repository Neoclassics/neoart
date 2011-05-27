
      subroutine PS(NSC,NCC,ns,nc,NEOFRC,tau,XI,LA,LAB,mas,temp,den,
     +              ZSP,sigma,uaioo,RNQ,uai,coeff)
C--------------------------------------------------------------------
C     THE ROUTINE THAT CALCULATES THE PFIRSCH SCHLUETER CONTRIBUTION.
C
C     INPUT   NSC    : MAXIMUM NUMBER OF SPECIES (USED TO CHECK 
C                      CONSISTENCY
C             NCC    : MAXIMUM NUMBER OF CHARGE STATES (USED TO 
C                      CHECK CONSISTENCY
C             NS     : NUMBER OF SPECIES
C             NC     : ARRAY(NS) NUMBER OF CHARGE STATES PER SPECIES
C             NEOFRC : LOGICAL IF TRUE THEN THE FRICTION MATRICES
C                      ARE NOT NEWLY CALCULATED
C             TAU    : ARRAY(NSC,NSC) NORMALIZED COLLISION FREQUEN-
C                      CIES (M_A N_A / TAU_AB)
C             XI     : ARRAY(NSC,NCC) RELATIVE WEIGHT OF THE CHARGE
C                      STATES
C             LA     : FRICTION COEFFICIENTS
C             LAB    : FRICTION COEFFICIENTS
C             MAS    : ARRAY(NSC) MASS OF THE SPECIES (KG)
C             TEMP   : ARRAY(NSC) TEMPERATURE IN KEV
C             DEN    : ARRAY(NSC,NCC) DENSITY IN 10^19 M^-3
C             ZSP    : CHARGE NORMALIZED TO E, ARRAY(NSM,NCM) 
C             SIGMA  : CONTROL PARAMETERS THAT DETERMINE WHICH OF 
C                      THE COUPLING TERMS ARE TAKEN INTO ACCOUNT. 
C                      FOR ALL COUPLING TERMS ALL SIGMAS ARE 1
C                      SEE ROUTINE NEOART FOR DEFINITION
C             UAIOO  : THE THERMODYNAMIC FORCES. ARRAY(NSM,NCM,3)
C                      SEE ROUTINE NEOART FOR DEFINITION
C             RNQ    : LENGTH OF THE FIELD LINE. IN CIRCULAR GEOMETRY
C                      THIS CAN BE APPROXIMATED BY MAJOR RADIUS 
C                      TIMES SAFETY FACTOR. 
C     OUTPUT  UAI    : THE VELOCITIES
C             COEFF  : NORMALIZED TRANSPORT COEFFICIENTS. 
C
C     THE ROUTINE CALLS THE FOLLOWING SUBROUTINES
C     PERR   : ERROR HANDLING
C     PENQ   : TO CALCULATE THE COUPLING COEFFICIENTS
C     LUDCMP : TO CALCULATE THE LU DECOMPOSITION
C     LUBKSB : TO CALCULATE THE SOLUTION OF THE MATRIX EQUATION   
C--------------------------------------------------------------------
      IMPLICIT NONE

      include 'elem_config.inc'
      
      integer NSM, NCM 	
      parameter(NSM = NELMAX+2)
      parameter(NCM = NIONMAX)

      integer i,j,k,l,indx, m, n, ns, nc, mt, nt, kt, lt, ms, 
     +        NSC, NCC
      real pab,qab,la,lab,tau,xi,kab,oa,ha,gk,
     +       aaa,a1,a2,a3,LG,da,dnom,rm,rmo,rmoi,
     +       dv,ev,fv,d,b,uo,uoo,sor,soro,soroo,
     +       uaio,uai,sigma,mas,temp,den,uaioo,
     +       dum1,dum2,pp,qq,ta,tb,ma,mb,lga,
     +       akeep,bkeep,sol,btus,f1,f2,RNQ,
     +       coeff,usolo,rmoo,rmooo,zsp
      logical nliter, nlchan, NEOFRC
      
      dimension pab(2,2,nsm,nsm),qab(2,2,nsm,nsm),
     +          la(3,3,nsm),lab(3,3,nsm,nsm),
     +          tau(nsm,nsm),xi(nsm,ncm),sigma(4),
     +          kab(2,2,nsm,nsm), oa(2,2,nsm),
     +          ha(2,3,nsm), gk(nsm), dv(2,3,nsm,nsm),
     +          aaa(2*nsm,2*nsm), ZSP(NSM,NCM),
     +          ev(2,3,nsm,nsm),fv(2,3,nsm,nsm),b(2*nsm),
     +          uaio(nsm,ncm,3),uo(nsm,3),uoo(nsm,3),
     +          sor(2,nsm,ncm),soro(2,nsm),uai(nsm,ncm,3),
     +          soroo(2,nsm),mas(nsm),uaioo(nsm,ncm,2),
     +          temp(nsm),den(nsm,ncm),
     +          pp(2,2),qq(2,2),coeff(nsm,ncm,4),
     +          akeep(2*nsm,2*nsm), bkeep(2*nsm), sol(2*nsm),
     +          btus(2*nsm),usolo(nsm,3),
     +          rm(2,2,nsm,ncm),rmo(2,2,nsm),rmoo(2,2,nsm),
     +          rmoi(2,2,nsm), rmooo(2,2,nsm)
      dimension indx(2*nsm), nc(nsm)

      SAVE AAA,INDX,HA,FV,DV,EV,RM,RMO,RMOO,RMOI

      nliter = .FALSE.
      nlchan = .false.
            
C     SOME COEFFICIENTS THAT APPEAR IN THE EQUATIONS
      a1 = 4./25.
      a2 = 16./175.
      a3 = (8./35.)**2

      IF ((NSM.NE.NSC).OR.(NCC.NE.NCM)) CALL PERR(1)

C     SINCE THE THERMODYNAMIC FORCES ARE PROPORTIONAL TO T/ Z 
C     THE UAIO QUANTITIES ARE MULTIPLIED BY THIS NUMBER.
      DO 8 I = 1, NS
        DO 8 J = 1, NC(I) 
          DO 9 K = 1, 2
            UAIO(I,J,K) = TEMP(I)*UAIOO(I,J,K) /(ZSP(I,J))
 9        CONTINUE
          UAIO(I,J,3) = 0.
 8    CONTINUE

C     IF THE ROUTINE IS CALLED WITH THE SAME DENSITIES AND TEMPERATURES
C     THEN DO NOT REPEAT THE CALCULATION OF THE MATRIX
      IF (.NOT.NEOFRC) THEN

C     THE COEFFICIENT THAT DETERMINES THE INFLUENCE OF THE FIELD
C     LINE LENGTH
      LG = -1./(RNQ)**2

C     CALCULATED THE FRICTION COEFFICIENTS
      DO 2 I = 1, NS
        DO 4 J = 1, NS
          TA = TEMP(I)
          TB = TEMP(J)
          MA = MAS(I)
          MB = MAS(J)       
          CALL PENQ(TA,TB,MA,MB,PP,QQ)
          DO 6 K = 1, 2
            DO 6 L = 1, 2
              PAB(K,L,I,J) = PP(K,L)
              QAB(K,L,I,J) = QQ(K,L)
 6          CONTINUE
 4      CONTINUE
 2    CONTINUE

c     calculate the kab and oa matrices
      do 100 i = 1, ns
        dum1 = 1./(temp(i))
        oa(1,1,i) = 0.
        oa(1,2,i) = 0.
        oa(2,1,i) = 0.
        oa(2,2,i) = 0.
        do 100 j = 1, ns
          dum2 = 1./(temp(j))
          kab(1,1,i,j) = -a1*dum1*sigma(3)*tau(i,j)*
     +                   pab(1,1,i,j)
          kab(1,2,i,j) = - a2*dum2*sigma(1)*tau(i,j)*
     +                   qab(1,2,i,j)
          kab(2,1,i,j) = a2*dum1*sigma(2)*tau(i,j)*
     +                   (pab(2,1,i,j)+pab(1,1,i,j))
          kab(2,2,i,j) = a3*dum2*tau(i,j)*sigma(4)*
     +                   (qab(2,2,i,j)+sigma(1)*qab(1,2,i,j))
          oa(1,1,i) = oa(1,1,i) + a1*dum1*sigma(3)*
     +                tau(i,j)*pab(1,1,i,j)
          oa(1,2,i) = oa(1,2,i) - a2*dum1*sigma(1)*
     +                tau(i,j)*pab(1,2,i,j)
          oa(2,1,i) = oa(2,1,i) - a2*dum1*sigma(2)*
     +                tau(i,j)*(pab(2,1,i,j) +
     +                pab(1,1,i,j))
          oa(2,2,i) = oa(2,2,i) + a3*dum1*tau(i,j)*sigma(4)*
     +                (pab(2,2,i,j)+ sigma(1)*pab(1,2,i,j))
100   continue

c     now calculate the ha coefficients
      do 200 i = 1, ns
        do 200 k = 1, 2
          do 200 l = 1, 3
            ha(k,l,i) = 0.
            do 300 j = 1, 2
              ha(k,l,i) = ha(k,l,i) - oa(k,j,i)*
     +                    la(j+1,l,i)
300         continue
200   continue

      do 400 i = 1, ns
        da = ha(1,2,i)*ha(2,3,i)-ha(1,3,i)*ha(2,2,i)
        ta = ha(1,2,i)+ha(2,3,i)
        lga = LG*mas(i)*1.6e22
        do 500 j = 1, nc(i)
          dnom = (xi(i,j)/den(i,j))/ (LGA**2 + (xi(i,j)/
     +       (den(i,j)))**2*LGA*ta+(xi(i,j)/(den(i,j)))**4*da)
          dum1 = (xi(i,j)/den(i,j))**2
          rm(1,1,i,j) =  dnom*(lga+dum1*ha(2,3,i))
          rm(1,2,i,j) = -dnom*dum1*ha(1,3,i)
          rm(2,1,i,j) = -dnom*dum1*ha(2,2,i)
          rm(2,2,i,j) =  dnom*(lga+dum1*ha(1,2,i))
500     continue
400   continue
 
      do 800 i = 1, ns
        gk(i) = 0.
        do 700 k = 1, 2
        do 701 l = 1, 2
          rmo(k,l,i) = 0.
          rmoo(k,l,i) = 0.
          rmooo(k,l,i) = 0.
 701    continue
 700    continue
        do 800 j = 1, nc(i)
          dum1 = xi(i,j)**2/den(i,j)
          dum2 = xi(i,j)*dum1/den(i,j)
          gk(i) = gk(i) + dum1
          do 750 k = 1, 2
          do 751 l = 1, 2
            rmo(k,l,i) = rmo(k,l,i)+ xi(i,j)*rm(k,l,i,j)
            rmoo(k,l,i) = rmoo(k,l,i)+dum1*rm(k,l,i,j)
            rmooo(k,l,i) = rmooo(k,l,i)+dum2*rm(k,l,i,j)
751       continue
750       continue
800   continue
 
      do 850 i = 1, ns
        dnom = 1./(rmo(1,1,i)*rmo(2,2,i)-
     +         rmo(1,2,i)*rmo(2,1,i))
        rmoi(1,1,i) =  dnom*rmo(2,2,i)
        rmoi(1,2,i) = -dnom*rmo(1,2,i)
        rmoi(2,1,i) = -dnom*rmo(2,1,i)
        rmoi(2,2,i) =  dnom*rmo(1,1,i)
 850  continue

      do 900 i = 1, ns
        do 900 j = 1, ns

          do 1000 k = 1, 2
            do 1000 l = 1, 3
              dv(k,l,i,j) = 0.
              ev(k,l,i,j) = 0.
              fv(k,l,i,j) = 0.
              do 1000 m = 1, 2
                dv(k,l,i,j) = dv(k,l,i,j) + 
     +                        kab(k,m,i,j)*la(m+1,l,j)
                ev(k,l,i,j) = ev(k,l,i,j) + 
     +                        oa(k,m,i)*lab(m+1,l,i,j)
                do 1000 n = 1, ns
                  fv(k,l,i,j) = fv(k,l,i,j) + 
     +               gk(n)*kab(k,m,i,n)*lab(m+1,l,n,j)
 1000     continue
 900  continue

c     build the matrix
      do 1100 i = 1, ns
        do 1100 j = 1, ns
          do 1100 k = 1, 2
            do 1100 l = 1, 2
              dum1 = 0.
              dum2 = 0.
              do 1120 n = 1, 2
                dum1 = dum1 + rmo(k,n,i)*fv(n,l+1,i,j)
     +                 +rmoo(k,n,i)*ev(n,l+1,i,j)
                do 1120 mt = 1, 2
                  do 1120 nt = 1, 2
                     dum2 = dum2 + rmo(k,n,i)*dv(n,mt+1,i,j)*
     +                 rmoo(mt,nt,j)*rmoi(nt,l,j)
 1120         continue
              aaa(2*(i-1)+k,2*(j-1)+l) = - dum1 - dum2

              dum1 = 0.
              do 1140 n = 1, ns
                do 1160 kt = 1, 2
                  do 1160 lt = 1, 2 
                    do 1160 mt = 1, 2
                      dum2 = 0.
                      do 1180 nt = 1, 2
                        do 1180 ms = 1, 2
                          dum2 = dum2 + rmoo(lt,nt,n)*
     +                      rmoi(nt,ms,n)*rmoo(ms,mt,n)
 1180                 continue
                      dum2 = rmooo(lt,mt,n)-dum2
                      dum1 = dum1 + rmo(k,kt,i)*dv(kt,lt+1,i,n)*
     +                    dum2*ev(mt,l+1,n,j)
 1160           continue 
 1140         continue
              aaa(2*(i-1)+k,2*(j-1)+l) = aaa(2*(i-1)+k,2*(j-1)+l)
     +           - dum1
 1100 continue

      do 1810 i = 1, 2*ns
        aaa(i,i) = aaa(i,i) + 1.
 1810 continue

      if (nlchan) then
      endif

      if (nliter) then
        do 1900 i = 1, 2*ns
          do 1900 j = 1, 2*ns
            akeep(i,j) = aaa(i,j)
 1900   continue
      endif

c     now do the lU decomposition 
      call ludcmp(aaa,2*ns,2*nsm,indx,d)

C     END THE IF STATEMENT THAT DETERMINES WHETER THE SAME 
C     MATRIX IS TO BE USED. 
      ENDIF

      do 2100 i = 1, ns
        do 2100 k = 1, 3
          uo(i,k) = 0.
          uoo(i,k) = 0.
          do 2120 n = 1, nc(i)
            uo(i,k) = uo(i,k) + xi(i,n)*uaio(i,n,k)
            uoo(i,k)=uoo(i,k)+xi(i,n)**2*uaio(i,n,k)/
     +         (den(i,n))
 2120     continue
 2100 continue

c     calculate the right hand side solutions
      do 2210 i = 1, ns
        do 2210 j = 1, nc(i)
          do 2210 k = 1, 2
            sor(k,i,j) = 0.
            do 2220 l = 1, 3
              do 2220 mt = 1, 2
                sor(k,i,j) = sor(k,i,j) - (xi(i,j)/den(i,j))*
     +            rm(k,mt,i,j)*ha(mt,l,i)*uaio(i,j,l)
                do 2230 n = 1, ns
                  sor(k,i,j) = sor(k,i,j) + rm(k,mt,i,j)*
     +            (fv(mt,l,i,n)*uo(n,l)+dv(mt,l,i,n)*uoo(n,l)+
     +            (xi(i,j)/den(i,j))*ev(mt,l,i,n)*uo(n,l))
 2230           continue
 2220       continue
 2210 continue

      do 2310 i= 1, ns
        do 2310 k = 1, 2
          soro(k,i) = 0.
          soroo(k,i) = 0.
          do 2320 j = 1, nc(i)
            soro(k,i) = soro(k,i) + xi(i,j)*sor(k,i,j)
            soroo(k,i) = soroo(k,i) + xi(i,j)**2*sor(k,i,j)/
     +                   den(i,j)
 2320     continue
 2310 continue

      do 2200 i = 1, ns
        do 2200 k = 1, 2
          b(2*(i-1)+k) = soro(k,i)
          do 2400 n = 1, ns
            do 2400 mt = 1, 2
              do 2400 nt = 1, 2
                b(2*(i-1)+k) = b(2*(i-1)+k) + rmo(k,mt,i)*
     +            dv(mt,nt+1,i,n)*soroo(nt,n)
                do 2450 kt = 1, 2
                  do 2450 lt = 1, 2
                    b(2*(i-1)+k) = b(2*(i-1)+k) - rmo(k,mt,i)*
     +                dv(mt,nt+1,i,n)*rmoo(nt,kt,n)*
     +                rmoi(kt,lt,n)*soro(lt,n)
 2450           continue
 2400     continue
 2200 continue
           

      if (nlchan) then
      endif

      if (nliter) then
        do 2411 i = 1, 2*ns
          bkeep(i) = b(i)
 2411   continue
      endif
      
c     calculate the solution through backsubstition   
      call lubksb(aaa,2*ns,2*nsm,indx,b)
 
      if (nliter) then
        do 1915 i = 1, 2*ns
          sol(i) = b(i)
 1915   continue
        do 1920 k=1, 5
          do 1930 i = 1, 2*ns
            btus(i) = bkeep(i)
            do 1930 j = 1, 2*ns
              btus(i) = btus(i)-akeep(i,j)*sol(j)
 1930     continue
          call lubksb(aaa,2*ns,2*nsm,indx,btus)
          do 1940 i = 1, 2*ns
            sol(i) = sol(i)+btus(i)
 1940     continue
 1920   continue
        do 1950 i = 1, 2*ns
          b(i) = sol(i)
 1950   continue
      endif

c     now calculate the velocities
      do 3000 i = 1, ns
        do 3000 j = 1, nc(i)
          do 3000 k = 1, 2
            uai(i,j,k+1) = sor(k,i,j)
            do 3100 mt = 1, 2
              do 3100 nt = 1, 2
                uai(i,j,k+1) = uai(i,j,k+1) + rm(k,mt,i,j)*
     +             rmoi(mt,nt,i) * (b(2*(i-1)+nt) - soro(nt,i))
 3100       continue
            do 3200 kt = 1, 2
              do 3200 lt = 1, 2
                dum1 = 0.
                if (kt.eq.lt) dum1 = dum1 + xi(i,j)/den(i,j)
                do 3300 mt = 1, 2
                  dum1 = dum1 - rmoi(kt,mt,i)*rmoo(mt,lt,i)
 3300           continue
                do 3400 mt = 1, 2 
                  do 3400 n = 1, ns
                    uai(i,j,k+1) = uai(i,j,k+1) + rm(k,kt,i,j)*
     +                dum1*ev(lt,mt+1,i,n)*b(2*(n-1)+mt)
 3400           continue
 3200       continue
 3000 continue
                  
  
      do 3500 i = 1, ns
        do 3500 j = 1, nc(i)
          uai(i,j,1) = uaio(i,j,1)
          uai(i,j,2) = uai(i,j,2)+uaio(i,j,2)
          usolo(i,1) = uo(i,1)
          usolo(i,2) = uo(i,2) + b(2*(i-1)+1)
          usolo(i,3) = b(2*(i-1)+2)
 3500 continue
 
c     calculate the friction 
      do 6000 i = 1, ns
        f1 = 0.
        f2 = 0.
        do 6100 j = 1, ns
          do 6100 k = 1, 3
            f1 = f1 + lab(1,k,i,j)*usolo(j,k)
            f2 = f2 + lab(2,k,i,j)*usolo(j,k)
 6100   continue
        do 6210 j = 1, nc(i)
          coeff(i,j,1) = 0.
          coeff(i,j,2) = 0.
          do 6200 k = 1, 3
            coeff(i,j,1) = coeff(i,j,1) + xi(i,j)*(la(1,k,i)*
     +                     uai(i,j,k))
            coeff(i,j,2) = coeff(i,j,2) + xi(i,j)*(la(2,k,i)*
     +                     uai(i,j,k))
 6200     continue
          coeff(i,j,1) = coeff(i,j,1) + xi(i,j)*f1
          coeff(i,j,2) = coeff(i,j,2) + xi(i,j)*f2
 6210   continue
 6000 continue
 
      return
      end
