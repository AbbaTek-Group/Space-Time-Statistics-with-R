      Subroutine DA(m,n,nn,XIN,din,RIN,ipick,XPOUT,nndepth,nntPTF)
	Dimension X(n,nn), d(m), R(m,m), 
     &          DD(m,nn), A(n,nn),C(n,n),
     &          V(n),AT(nn,n),HA(m,nn),HX(m,nn),w(m),
     $          HAT(nn,m), P(m,m),XP(n,nn),del(m,nn),PIN(m,m),
     &          AUX1(m,nn), AUX2(nn,nn), AUX3(n,nn), PCOPY(m,m)
      Dimension XIN(nndepth,nntPTF),din(nndepth),RIN(nndepth,nndepth)
	Dimension XPOUT(nndepth,nntPTF),    
     &          RINCOPY(m,m),diag(n) !FP for more PTFs

      Dimension indx(m)
      Dimension ipick(nndepth)
c============================
c Internal variables
c XIN -> X,  din->d,  RIN-> R
c============================
	Do i=1,n
	Do ii=1,nn
	  X(i,ii)=XIN(i,ii)
	Enddo
	Enddo
      Do j=1,m
	  d(j)=din(j)
	Enddo
      Do j1=1,m
      Do j2=1,m
        R(j1,j2)=RIN(j1,j2)
	Enddo
      Enddo
c;=================================================================
c     X is an n x NN matrix whose columns are the ensemble members, 
c     and it is called the prior ensemble.
c     n the number of state variables 
c     The data d is assumed to have Gaussian pdf with covariance R and mean Hx, 
c     where H is the so-called the observation matrix. 
c     The covariance matrix R describes the estimate of the error of the data; 
c     If the random errors in the entries of the data vector d are independent, 
c     R is diagonal and its diagonal entries are the squares of the standard deviation 
c     (“error size”) of the error of the corresponding entries of the data vector d. 
c     The value Hx is what the value of the data would be for the state x 
c     in the absence of data errors. 


c;=====================================================================
c;Replicate the data d into an m x N matrix DD
c;DD=[d_1,D-2,...D_NN], d_i=d+eps_i, eps_i is from N(0,R)
c; This is how one samples the multivariate distribution:
c; http://home.jesus.ox.ac.uk/~clifford/a5/chap1/node16.html
c==========================
c Copy for Choleski decomp
c===========================
      Do k1=1,m
	Do k2=1,m
        RINCOPY(k1,k2)=RIN(k1,k2)
      Enddo
	Enddo
	call choldc(RINCOPY,m,m,diag)
	idum = 951
      Do kk=1,NN
        Do j=1,m
	    DD(j,kk)=d(j)+diag(j)*gasdev(idum)
	    Do k1=1,j-1
		  DD(j,kk)=DD(j,kk)+RINCOPY(j,k1)*gasdev(idum)
	    Enddo
        Enddo
      Enddo

c========================================
c     Matrix A is defined as X-E(X)
c     A=X-E(x)=X-(X e_NN x 1)e_1 x NN /NN
c     Matrix C is AA^T/(NN-1)
c     Essentially covariance matrix
c========================================
      Do i=1,n 
        s=0.
        Do kk=1,NN 
          s=s+X(i,kk)
        Enddo
        V(i)=s/float(NN)
      Enddo

      Do kk=1,NN
        Do i=1,n
          A(i,kk)=X(i,kk)-V(i)
        Enddo
      Enddo
c
      Do i=1,n 
	  Do kk=1,NN 
           AT(kk,i)=A(i,kk)
        Enddo
      Enddo
c
      Do i1=1,n
        Do i2=1,n
          C(i1,i2)=0.
          Do kk=1,NN
             C(i1,i2)=C(i1,i2) + A(i1,kk)*AT(kk,i2)
          Enddo
          C(i1,i2)=C(i1,i2)/float(NN-1)
        Enddo
      Enddo

c==========================================================================================
c     Observation matrix-free implementation
c
c     Function h(x) of the form h(x)=H(x) is more natural to compute
c     The function h is called the observation function or, in the inverse problems context,
c     the forward operator. 
c     The value of h(x) is what the value of the data would be for the state x 
c     assuming the measurement is exact.
c===========================================================================================
	mm=1
      Do j=1,n
        if(ipick(j).GT.0) then
c-----h(x)=x exact observation
        Do kk=1,NN
           HX(mm,kk)=X(j,kk)
        Enddo
 	  mm=mm+1
       Endif
      Enddo
c 
      Do j=1,m
        s=0.
        Do kk=1,NN
          s=s+HX(j,kk)
        Enddo
        W(j)=s/float(NN)
      Enddo
c
      Do kk=1,NN
        Do j=1,m
          HA(j,kk)=HX(j,kk)-W(j)
        Enddo
      Enddo

      Do j=1,m
        Do kk=1,NN
          HAT(kk,j)=HA(j,kk)
        Enddo
      Enddo


      Do j1=1,m
        Do j2=1,m
          P(j1,j2)=0.
          Do kk=1,NN
            P(j1,j2)=P(j1,j2) + HA(j1,kk)*HAT(kk,j2)
          Enddo
          P(j1,j2)=P(j1,j2)/float(NN-1) + R(j1,j2)
        Enddo
      Enddo
c;==========
c; Posterior
c;==========

      Do j=1,m
        Do kk=1,NN
          del(j,kk)=DD(j,kk)-HX(j,kk)
        Enddo
      Enddo
c
	Do j1=1,m
	Do j2=1,m
	 PCOPY(j1,j2)=P(j1,j2)
	Enddo
	Enddo

	Do j1=1,m
	  Do j2=1,m
	    PIN(j1,j2)=0.
	  Enddo
	  PIN(j1,j1)=1.0
	Enddo
	CALL ludcmp(PCOPY,m, m, indx, vd)
	Do i2=1,m
        CALL lubksb(PCOPY,m,m,indx,PIN(1,i2))
	Enddo
      Do j=1,m
        Do kk=1,NN
           AUX1(j,kk)=0.
           Do j3=1,m
              AUX1(j,kk)=AUX1(j,kk)+PIN(j,j3)*del(j3,kk)
           Enddo
        Enddo
      Enddo
c
      Do kk1=1,NN
        Do kk2=1,NN
           AUX2(kk1,kk2)=0.
           Do j=1,m
              AUX2(kk1,kk2)= AUX2(kk1,kk2) + HAT(kk1,j) * AUX1(j,kk2)
           Enddo
        Enddo
      Enddo
c
       Do i=1,n
        Do kk=1,NN
           AUX3(i,kk)=0.
           Do kk3=1,NN
              AUX3(i,kk) = AUX3(i,kk) + A(i,kk3)*AUX2(kk3,kk)
           Enddo
        Enddo
      Enddo

c======================
c Finally the posterior
c======================

      Do i=1,n
        Do kk=1,NN
           XP(i,kk)=X(i,kk)+AUX3(i,kk)/float(NN-1)
	     XPOUT(i,kk)=XP(i,kk)
        Enddo
      Enddo

      return
      end





      FUNCTION gasdev(idum)
      INTEGER idum
      REAL gasdev
CU    USES ran1
      INTEGER iset
      REAL fac,gset,rsq,v1,v2,ran1
      SAVE iset,gset
      DATA iset/0/
      if (iset.eq.0) then
1       v1=2.*ran1(idum)-1.
        v2=2.*ran1(idum)-1.
        rsq=v1**2+v2**2
        if(rsq.ge.1..or.rsq.eq.0.)goto 1
        fac=sqrt(-2.*log(rsq)/rsq)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
      else
        gasdev=gset
        iset=0
      endif
      return
      END

	FUNCTION ran1(idum)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL ran1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     *NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM

      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
      END

	SUBROUTINE ludcmp(a,n,np,indx,d)
      INTEGER n,np,indx(n),NMAX
      REAL d,a(np,np),TINY
      PARAMETER (NMAX=500,TINY=1.0e-20)
      INTEGER i,imax,j,k
      REAL aamax,dum,sum,vv(NMAX)
      d=1.
      do 12 i=1,n
        aamax=0.
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
        if (aamax.eq.0.) pause 'singular matrix in ludcmp'
        vv(i)=1./aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.
        do 16 i=j,n

          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.)a(j,j)=TINY
        if(j.ne.n)then
          dum=1./a(j,j)

          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      return
      END

      SUBROUTINE lubksb(a,n,np,indx,b)
      INTEGER n,np,indx(n)
      REAL a(np,np),b(n)
      INTEGER i,ii,j,ll
      REAL sum
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      END
c==========================
c Choleski decomposition
c==========================
      SUBROUTINE choldc(a,n,np,p)
      INTEGER n,np
      REAL a(np,np),p(n)
      INTEGER i,j,k
      REAL sum
      do 13 i=1,n
        do 12 j=i,n
          sum=a(i,j)
          do 11 k=i-1,1,-1
            sum=sum-a(i,k)*a(j,k)
11        continue
          if(i.eq.j)then
            if(sum.le.0.)pause 'choldc failed'
            p(i)=sqrt(sum)
          else
            a(j,i)=sum/p(i)
          endif
12      continue
13    continue
      return
      END