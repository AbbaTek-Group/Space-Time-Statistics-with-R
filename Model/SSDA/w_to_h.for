      SUBROUTINE W_TO_H(n,x,h,mat,k,xobs,wobs,thr,ths,alpha,vgn,imodel)
	Dimension x(n),h(n),xobs(k),wobs(k),
     &          thr(k),ths(k),alpha(k),vgn(k),mat(n)
	Dimension hlgobs(k),m(k)
c==========================================================
c     This routine is used only once when no information outside
c     the range of depths [xobs(1), xobs(n)] is available
c==========================================================
        
c======================================
c  Nodal numbers for reference depths
c======================================
      Do j=1,k
        difmax=1000.
        Do i=1,n 
        delt=abs(x(i)-xobs(j))
        If(delt.LT.difmax) then
          difmax=delt
          m(j)=i
        Endif
        Enddo
      Enddo
c=============================================================================
c Get logarithms of pressure head - they will be interpolated and extrapolated
c=============================================================================
 	if (imodel.eq.0) then !Van genuchten model
	Do jj=1,k
	  j=mat(m(jj))
c======================================================
c Trim observed water contents to stay within the expected range
c======================================================
        if(wobs(jj).GT.ths(j)) wobs(jj)=0.999*ths(j)
	  if(wobs(jj).LT.thr(j)) wobs(jj)=1.001*thr(j)
	  os=(ths(j)-thr(j))/(wobs(jj)-thr(j))
	  hobspt=(os**(1.0/(1.0-1.0/vgn(j)))-1.0)**(1.0/vgn(j))
	  hobspt=hobspt/alpha(j)
	  hlgobs(jj)=alog10(hobspt)
	Enddo
	elseif(imodel.eq.2) then !Brooks-Corey model
	 Do jj =1,k
c======================================================
c Trim observed water contents to stay within the expected range
c======================================================
        if(wobs(jj).GT.ths(j)) wobs(jj)=0.999*ths(j)
	  if(wobs(jj).LT.thr(j)) wobs(jj)=1.001*thr(j)
	  j=mat(m(jj))
	  os=(ths(j)-thr(j))/(wobs(jj)-thr(j))
	  hobspt=os**(1.0/vgn(j))
	  hobspt=hobspt/alpha(j)
	  hlgobs(jj)=alog10(hobspt) 
	 enddo
	endif !FP    
      Do i=1,n
c==============
c Extrapolation
c==============
      if(i.LT.m(1)) then
          hlg=hlgobs(1)+(x(i)-x(m(1)))*
     &    (hlgobs(2)-hlgobs(1))/(x(m(2))-x(m(1)))
          goto 2
        endif
        if(i.GT.m(k)) then
          hlg=hlgobs(k-1)+(x(i)-x(m(k-1)))*
     &    (hlgobs(k)-hlgobs(k-1))/(x(m(k))-x(m(k-1)))
        goto 2
      endif
c==============
c Interpolation
c==============
      Do j=1,k-1
        if(i.GE.m(j).AND.i.LE.m(j+1)) then
        hlg=hlgobs(j)+(x(i)-x(m(j)))*(hlgobs(j+1)-hlgobs(j))
     &      /(x(m(j+1))-x(m(j)))
        goto 2
        endif
      Enddo
 2      h(i)=-10.0**hlg
	Enddo
      Return
	End

