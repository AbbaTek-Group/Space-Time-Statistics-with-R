* Source file SINK.FOR |||||||||||||||||||||||||||||||||||||||||||||||||

      subroutine SetSnk(N,NMat,MatNum,x,hRoot,vRoot,Sink,TPot,hNew,
     !                  lMoSink,lSolRed,lSolAdd,P0,POptm,P2H,P2L,P3,r2H,
     !                  r2L,aOsm,c50,P3c,Beta,lChem,NS,NSD,Conc,cRoot)

      logical lChem,lMoSink,lSolRed,lSolAdd
      dimension x(N),MatNum(N),hNew(N),POptm(NMat),Beta(N),Sink(N),
     !          Conc(NSD,N),cRoot(NS),aOsm(NS)

      vRoot=0.
      hRoot=0.
      ARoot=0.
      do 11 ii=1,NS
        cRoot(ii)=0.
11    continue
      do 13 i=2,N
        if(Beta(i).gt.0.) then
          if(i.eq.N) then
            dxM=(x(i)-x(i-1))/2.
          else
            dxM=(x(i+1)-x(i-1))/2.
          end if
          M=MatNum(i)
          hRed=hNew(i)
          SAlfa=1.
          if(lChem.and.lSolRed) then
            cRed=0.
            do 15 j=1,NS
              cRed=cRed+aOsm(j)*Conc(j,i)
15          continue
            if(lSolAdd) then
              hRed=hRed+cRed
            else
              SAlfa=FSAlfa(cRed,c50,P3c)
            end if
          end if
          Alfa=FAlfa(lMoSink,TPot,hRed,P0,POptm(M),P2H,P2L,P3,r2H,r2L)
          Sink(i)=Alfa*SAlfa*Beta(i)*TPot
          vRoot=vRoot+Sink(i)*dxM
          hRoot=hRoot+hNew(i)*dxM
          do 12 ii=1,NS
            if(lChem) cRoot(ii)=cRoot(ii)+Conc(ii,i)*dxM
12        continue
          ARoot=ARoot+dxM
        else
          Sink(i)=0.
        end if
13    continue
      if(ARoot.gt.0.001) then
        hRoot=hRoot/ARoot
        do 14 ii=1,NS
          cRoot(ii)=cRoot(ii)/ARoot
14      continue
      end if
      return
      end

************************************************************************

      real function FSAlfa(cRed,c50,P3c)

      FSAlfa=1./(1.+(cRed/c50)**P3c)
      return
      end

************************************************************************

      real function FAlfa(lMoSink,TPot,h,P0,P1,P2H,P2L,P3,r2H,r2L)

      logical lMoSink


      if(lMoSink) then
        if(TPot.lt.r2L) P2=P2L
        if(TPot.gt.r2H) P2=P2H
        if((TPot.ge.r2L).and.(TPot.le.r2H))
     !    P2=P2H+(r2H-TPot)/(r2H-r2L)*(P2L-P2H)
        FAlfa=0.0
        if((h.gt.P3).and.(h.lt.P2)) FAlfa=(h-P3)/(P2-P3)
        if((h.ge.P2).and.(h.le.P1)) FAlfa=1.0
        if((h.gt.P1).and.(h.lt.P0)) FAlfa=(h-P0)/(P1-P0)
      else
        FAlfa=1./(1.+(h/P0)**P3)
      end if
      return
      end

************************************************************************

      subroutine SetRG(NumNP,x,Beta,t,tRMin,tRHarv,xRMin,xRMax,RGR)

      dimension x(NumNP),Beta(NumNP)

      if(t.lt.tRMin.or.t.gt.tRHarv) then
        do 11 i=1,NumNP
          Beta(i)=0.
11      continue
        return
      end if
      xR=xRMax
      if(xRMin.le.1.) xRMin=1.
      tt=t-tRMin
      xR=(xRMax*xRMin)/(xRMin+(xRMax-xRMin)*exp(-RGR*tt))
      SBeta=0.
      do 12 i=2,NumNP-1
        if(x(i).lt.x(NumNP)-xR) then
          Beta(i)=0.
        else if(x(i).lt.x(NumNP)-0.2*xR) then
          Beta(i)=2.08333/xR*(1-(x(NumNP)-x(i))/xR)
        else
          Beta(i)=1.66667/xR
        end if
        if(i.ne.NumNP) then
          SBeta=SBeta+Beta(i)*(x(i+1)-x(i-1))/2.
        else
          SBeta=SBeta+Beta(i)*(x(i)-x(i-1))/2.
        end if
12    continue
      if(SBeta.lt.0.0001) then
        Beta(NumNP-1)=1.
      else
        do 13 i=2,NumNP-1
          Beta(i)=Beta(i)/SBeta
13      continue
      end if
      return
      end

* ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||