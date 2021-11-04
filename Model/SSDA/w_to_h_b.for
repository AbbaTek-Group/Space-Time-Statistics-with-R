      SUBROUTINE W_TO_H_B(n,x,h,mat,k,xobs,wobs,thr,ths,alpha,vgn,
     CThetasav,imodel)
	Dimension x(n),h(n),xobs(k),wobs(k),
     &          thr(k),ths(k),alpha(k),vgn(k),mat(n),Thetasav(n),
     &          hlgThetasav(n)
	Dimension hlgobs(k),m(k)
c===========================================================================
c the purpose of having thetasav is to make the dependence 
c of the log pressure head on the distance outside of the interpolation range
c similar to that in the Thetasav
c============================================================================ 
c=======================================
c Nodal numbers for the reference depths
c=======================================
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
 	if (imodel.eq.0) then !van Genuchten's model
	Do jj=1,k
	  j=mat(m(jj))
        if(wobs(jj).GE.ths(j)) wobs(jj)=0.999*ths(j)
	  if(wobs(jj).LE.thr(j)) wobs(jj)=1.001*thr(j)
	  os=(ths(j)-thr(j))/(wobs(jj)-thr(j))
	  hobspt=(os**(1.0/(1.0-1.0/vgn(j)))-1.0)**(1.0/vgn(j))
	  hobspt=hobspt/alpha(j)
	  hlgobs(jj)=alog10(hobspt)
	Enddo
	elseif(imodel.eq.2) then !Brooks-Corey model
	 Do jj =1,k
	  j=mat(m(jj))
        if(wobs(jj).GE.ths(j)) wobs(jj)=0.999*ths(j)
	  if(wobs(jj).LE.thr(j)) wobs(jj)=1.001*thr(j)
	  os=(ths(j)-thr(j))/(wobs(jj)-thr(j))
	  hobspt=os**(1.0/vgn(j))
	  hobspt=hobspt/alpha(j)
	  hlgobs(jj)=alog10(hobspt) 
	 enddo
	endif  
c==================================================
c logpressure heads for the Thetasav water contents
c==================================================   
	Do i=1,n
	  j=mat(i)
       if (imodel.eq.0) then 
       if(Thetasav(i).GE.ths(j)) Thetasav(i)=0.999*ths(j)
	  if(Thetasav(i).LE.thr(j)) Thetasav(i)=1.001*thr(j)
	  os=(ths(j)-thr(j))/(Thetasav(i)-thr(j))
	  hobspt=(os**(1.0/(1.0-1.0/vgn(j)))-1.0)**(1.0/vgn(j))
	  hobspt=hobspt/alpha(j)
	  hlgThetasav(i)=alog10(hobspt)
	elseif(imodel.eq.2) then 
	  j=mat(i)
        if(Thetasav(i).GT.ths(j)) Thetasav(i)=0.999*ths(j)
	  if(Thetasav(i).LT.thr(j)) Thetasav(i)=1.001*thr(j)
	  os=(ths(j)-thr(j))/(Thetasav(i)-thr(j))
	  hobspt=os**(1.0/vgn(j))
	  hobspt=hobspt/alpha(j)
	  hlgThetasav(i)=alog10(hobspt) 
	endif !FP 
	enddo
c
      Do i=1,n
c=============================================================
c Pressure head on soil surface is set the same as in Thetasav
c=============================================================
        if(i.LT.m(1)) then
        hlg=hlgThetasav(i)+
     &         (hlgobs(1)-hlgThetasav(i))*(x(i)-x(1))/(x(m(1))-x(1))
          goto 2
        endif
c=============================================================
c Pressure head on st the bottom of the profile is set the same 
c as in Thetasav
c=============================================================
        if(i.GT.m(k)) then
        hlg=hlgThetasav(n)+
     &         (hlgobs(k)-hlgThetasav(n))*(x(i)-x(n))/(x(m(k))-x(n))
        goto 2
        endif
c===============
c Interpolation
c===============
       Do j=1,k-1
          if(i.GE.m(j).AND.i.LE.m(j+1)) then
        if(idep.EQ.0) then
          hlg=hlgobs(j)+float((i-m(j)))*(hlgobs(j+1)-hlgobs(j))
     &      /float((m(j+1)-m(j)))
	  else
          hlg=hlgobs(j)+(x(i)-x(m(j)))*(hlgobs(j+1)-hlgobs(j))
     &      /(x(m(j+1))-x(m(j)))
	  endif
          goto 2
          endif
        Enddo

 2      h(i)=-10.0**hlg
	Enddo
      Return
	End

