* Source file WATFLOW.FOR ||||||||||||||||||||||||||||||||||||||||||||||

      subroutine WatFlow(NumNP,NTab,NTabD,NMat,hTab,ConTab,CapTab,hNew,
     !                   hOld,MatNum,ParD,ParW,Con,Cap,ConSat,Ah,AK,ATh,
     !                   hSat,hTemp,KodTop,KodBot,rTop,rBot,CosAlf,t,dt,
     !                   x,Sink,P,R,S,FreeD,SeepF,qGWLF,Aqh,Bqh,GWL0L,
     !                   hTop,hBot,hCritA,hCritS,WLayer,Iter,ItCum,
     !                   TopInf,KTOld,KBOld,TolTh,TolH,MaxIt,dtMin,tOld,
     !                   dtOpt,ConvgF,TheTab,ThNew,ThOld,thr,ths,lWDep,
     !                   TempN,Kappa,KappaO,AThS,ThRR,ConO,ConR,AKS,AhW,
     !                   AThW,AKW,iHyst,iModel)

      logical ConvgF,ItCrit,FreeD,qGWLF,TopInf,WLayer,SeepF,lWDep
      double precision P,R,S,PB,RB,SB,PT,RT,ST,rMin
      dimension x(NumNP),hNew(NumNP),hOld(NumNP),hTemp(NumNP),thr(NMat),
     !          MatNum(NumNP),ParD(9,NMat),Sink(NumNP),hTab(NTabD),
     !          ConTab(NTabD,NMat),CapTab(NTabD,NMat),Con(NumNP),
     !          Cap(NumNP),ConSat(NMat),P(NumNP),R(NumNP),S(NumNP),
     !          hSat(NMat),Ah(NumNP),AK(NumNP),ATh(NumNP),ths(NMat),
     !          TheTab(NTabD,NMat),ThNew(NumNP),ThOld(NumNP),
     !          TempN(NumNP),ParW(9,NMat),ConO(NumNP),Kappa(NumNP),
     !          AThS(NumNP),ThRR(NUmNP),ConR(NumNP),AKS(NumNP),
     !          AhW(NMat),AThW(NMat),KappaO(NumNP),AKW(NMat)

      rMax=1.e+10
      rMin=1.d-100
11    continue

      Iter=0
      ConvgF=.true.

12    continue

*     Generate terms of matrix equation and solve by Gauss elimination
      call SetMat(NumNP,NTab,NTabD,NMat,hTab,ConTab,CapTab,hNew,MatNum,
     !            ParD,Con,Cap,ConSat,Ah,AK,ATh,hSat,hTemp,TheTab,ThNew,
     !            thr,ths,lWDep,TempN,Iter,ConO,Kappa,AThS,ThRR,ConR,
     !            AKS,AhW,AThW,AKW,iModel)
      if(Iter.eq.1.and.iHyst.gt.0)
     !call Hyster(NumNP,NMat,hOld,MatNum,ParD,ParW,ThNew,ThOld,Kappa,
     !            AThS,ThRR,ConO,ConR,AKS,KappaO,Ah,AK,iHyst,iModel)
      call Reset (NumNP,rTop,rBot,CosAlf,dt,x,hOld,Con,Cap,WLayer,hNew,
     !            Sink,P,R,S,PB,RB,SB,PT,RT,ST,FreeD,qGWLF,Aqh,Bqh,
     !            GWL0L,ThNew,ThOld,vTop)
      call Shift (KodTop,rTop,rBot,hTop,hBot,hCritA,hCritS,CosAlf,
     !            WLayer,Con(NumNP),Con(NumNP-1),hNew(NumNP),
     !            hNew(NumNP-1),x(NumNP),x(NumNP-1),TopInf,KodBot,
     !            Con(1),Con(2),hNew(1),hNew(2),x(1),x(2),SeepF,
     !            ThNew(NumNP),ThOld(NumNP),Sink(NumNP),ThNew(1),
     !            ThOld(1),Sink(1),dt)
      do 13 i=1,NumNP
        hTemp(i)=hNew(i)
13    continue
      call Gauss (NumNP,KodTop,KodBot,hTop,hBot,hNew,P,R,S,PB,RB,SB,PT,
     !            RT,ST,rMin,rMax)
      Iter =Iter+1
      ItCum=ItCum+1

*     Test for convergence
      ItCrit=.true.
      do 14 i=1,NumNP
        m=MatNum(i)
        EpsTh=0.
        EpsH=0.
        if(hTemp(i).lt.hSat(m).and.hNew(i).lt.hSat(m)) then
          Th=ThNew(i)+Cap(i)*(hNew(i)-hTemp(i))/(ths(m)-thr(m))/ATh(i)
          EpsTh=abs(ThNew(i)-Th)
        else
          EpsH=abs(hNew(i)-hTemp(i))
        end if
        if(EpsTh.gt.TolTh.or.EpsH.gt.TolH) then
          ItCrit=.false.
          goto 15
        end if
14    continue
15    continue
      if(ItCrit) then
        do 26 i=1,NumNP
          ThNew(i)=ThNew(i)+Cap(i)*(hNew(i)-hTemp(i))
26      continue
      end if
      if(.not.ItCrit) then
        if(Iter.lt.MaxIt) then
          goto 12
        else if(dt.le.dtMin) then
          ConvgF=.false.
          write(*,*) ' The numerical solution has not converged ! '
	    write(*,*) dt,dtmin
	    write(*,*) Iter
          return
        else
          do 16 i=1,NumNP
            if(iHyst.gt.0) Kappa(i)=KappaO(i)
            hNew(i) =hOld(i)
            hTemp(i)=hOld(i)
16        continue
          KodTop=KTOld
          KodBot=KBOld
          dt=amax1(dt/3,dtMin)
          dtOpt=dt
          t=tOld+dt
          goto 11
        end if
      end if
      return
      end

************************************************************************

      subroutine Reset(N,rTop,rBot,CosAlf,dt,x,hOld,Con,Cap,WLayer,hNew,
     !                 Sink,P,R,S,PB,RB,SB,PT,RT,ST,FreeD,qGWLF,Aqh,Bqh,
     !                 GWL0L,ThNew,ThOld,vTop)

      logical WLayer,FreeD,qGWLF
      double precision P,R,S,PB,RB,SB,PT,RT,ST,A2,A3,B,F2
      dimension x(N),hOld(N),hNew(N),P(N),R(N),S(N),Con(N),Cap(N),
     !          Sink(N),ThNew(N),ThOld(N)

*     Finite differences

      dxB=x(2)-x(1)
      dx=dxB/2.
      ConB=(Con(1)+Con(2))/2.
      S(1)=-ConB/dxB
      if(FreeD) then
        RB=ConB/dxB
        SB=-RB
        PB=0.d0
      else
        RB=ConB/dxB+Cap(1)*dx/dt
        SB=-ConB/dxB
        if(qGWLF) rBot=Fqh(hNew(1)-GWL0L,Aqh,Bqh)
        PB=ConB*CosAlf-Sink(1)*dx+Cap(1)*dx*hNew(1)/dt-
     !     (ThNew(1)-ThOld(1))*dx/dt+rBot
      end if
      do 11 i=2,N-1
        dxA=x(i)-x(i-1)
        dxB=x(i+1)-x(i)
        dx=(dxA+dxB)/2.
        ConA=(Con(i)+Con(i-1))/2.
        ConB=(Con(i)+Con(i+1))/2.
        A2=ConA/dxA+ConB/dxB
        A3=-ConB/dxB
        B =(ConA-ConB)*CosAlf
        F2=Cap(i)*dx
        R(i)=A2+F2/dt
        P(i)=F2*hNew(i)/dt-(ThNew(i)-ThOld(i))*dx/dt-B-Sink(i)*dx
        S(i)=A3
11    continue
      dxA=x(N)-x(N-1)
      dx=dxA/2.
      ConA=(Con(N)+Con(N-1))/2.
      RT=ConA/dxA+Cap(N)*dx/dt
      ST=-ConA/dxA
      PT=Cap(N)*dx*hNew(N)/dt-(ThNew(N)-ThOld(N))*dx/dt-Sink(N)*dx-
     !   ConA*CosAlf
      vTop=-sngl(ST)*hNew(N-1)-sngl(RT)*hNew(N)+sngl(PT)
      PT=PT-rTop
      if(WLayer) then
        if(hNew(N).gt.0.) then
          RT=RT+1./dt
          PT=PT+amax1(hOld(N),0.)/dt
        else
          PT=PT+amax1(hOld(N),0.)/dt
        end if
      end if
      return
      end

***********************************************************************

      subroutine Gauss(N,KodTop,KodBot,hTop,hBot,hNew,P,R,S,PB,RB,SB,PT,
     !                 RT,ST,rMin,rMax)

      double precision P,R,S,PB,RB,SB,PT,RT,ST,rMin
      dimension hNew(N),P(N),R(N),S(N)

*     Forward
      if(KodBot.ge.0) then
        P(2)=P(2)-S(1)*hBot
      else
        if(dabs(RB).lt.rMin) RB=rMin
        P(2)=P(2)-PB*S(1)/RB
        R(2)=R(2)-SB*S(1)/RB
      end if
      do 11 i=3,N-1
        if(dabs(R(i-1)).lt.rMin) R(i-1)=rMin
        P(i)=P(i)-P(i-1)*S(i-1)/R(i-1)
        R(i)=R(i)-S(i-1)*S(i-1)/R(i-1)
11    continue
      if(KodTop.ge.0) then
        P(N-1)=P(N-1)-S(N-1)*hTop
      else
        if(dabs(R(N-1)).lt.rMin) R(N-1)=rMin
        P(N)=PT-P(N-1)*ST/R(N-1)
        R(N)=RT-S(N-1)*ST/R(N-1)
      end if

*     Back
      if(dabs(R(N-1)).lt.rMin) R(N-1)=rMin
      if(KodTop.ge.0) then
        hNew(N)=hTop
        hNew(N-1)=sngl(P(N-1)/R(N-1))
      else
        hNew(N)=sngl(P(N)/R(N))
        hNew(N-1)=sngl((P(N-1)-S(N-1)*hNew(N))/R(N-1))
      end if
      do 12 i=N-2,2,-1
        if(dabs(R(i)).lt.rMin) R(i)=rMin
        hNew(i)=sngl((P(i)-S(i)*hNew(i+1))/R(i))
12    continue
      if(KodBot.ge.0) then
        hNew(1)=hBot
      else
        if(dabs(RB).lt.rMin) RB=rMin
        hNew(1)=sngl((PB-SB*hNew(2))/RB)
      end if
      do 13 i=1,N
        if(abs(hNew(i)).gt.rMax) hNew(i)=sign(rMax,hNew(i))
13    continue        
      return
      end

***********************************************************************

      subroutine Shift(KodTop,rTop,rBot,hTop,hBot,hCritA,hCritS,CosAlf,
     !                 WLayer,ConTop,ConBlw,hNT,hNB,xTop,xBlw,TopInf,
     !                 KodBot,ConBot,ConAbv,hN1,hN2,xBot,xAbv,SeepF,
     !                 ThNewT,ThOldT,SinkTop,ThNewB,ThOldB,SinkBot,dt)

      logical WLayer,TopInf,SeepF

*     Seepage face at the bottom
      if(SeepF) then
        if(KodBot.ge.0) then
          vBot=-(ConBot+ConAbv)/2.*((hN2-hN1)/(xAbv-xBot)+CosAlf)-
     !          (xAbv-xBot)/2.*((ThNewB-ThOldB)/dt+SinkBot)
          if(vBot.gt.0.) then
            KodBot=-2
            rBot=0.
          end if
        else
          if(hN1.ge.0.) then
            KodBot=2
            hBot=0.
          end if
        end if
      end if

*     Atmospheric boundary condition
      if(TopInf) then
        if(KodTop.ge.0) then
          vTop=-(ConTop+ConBlw)/2.*((hNT-hNB)/(xTop-xBlw)+CosAlf)-
     !         (ThNewT-ThOldT)*(xTop-xBlw)/2./dt-SinkTop*(xTop-xBlw)/2.
          if(abs(vTop).gt.abs(rTop).or.vTop*rTop.le.0) KodTop=-4
        else
          if(.not.WLayer) then
            if(hNT.gt.hCritS) then
              KodTop=4
              hTop=hCritS
            end if
          end if
          if(hNT.le.hCritA) then
            KodTop=4
            hTop=hCritA
          end if
        end if
      end if
      return
      end

************************************************************************

      subroutine SetMat(NumNP,NTab,NTabD,NMat,hTab,ConTab,CapTab,hNew,
     !                  MatNum,ParD,Con,Cap,ConSat,Ah,AK,ATh,hSat,hTemp,
     !                  TheTab,theta,thr,ths,lWDep,TempN,Iter,ConO,
     !                  Kappa,AThS,ThRR,ConR,AKS,AhW,AThW,AKW,iModel)

      logical lWDep
      dimension hTab(NTab),ConTab(NTabD,NMat),CapTab(NTabD,NMat),
     !          hNew(NumNP),MatNum(NumNP),ParD(9,NMat),Con(NumNP),
     !          Cap(NumNP),ConSat(NMat),Ah(NumNP),AK(NumNP),
     !          ATh(NumNP),hSat(NMat),hTemp(NumNP),TheTab(NTabD,NMat),
     !          theta(NumNP),thr(NMat),ths(NMat),TempN(NumNP),
     !          ConO(NumNP),Kappa(NumNP),AThS(NumNP),ThRR(NUmNP),
     !          ConR(NumNP),AKS(NumNP),AhW(NMat),AThW(NMat),AKW(NMat)

      alh1=alog10(-hTab(1))
      dlh =(alog10(-hTab(NTab))-alh1)/(NTab-1)
      do 11 i=1,NumNP
        if(lWDep.and.Iter.eq.0) then
          Temp=20.
          Ah(i)=(75.6-0.1425*TempN(i)-2.38e-4*TempN(i)**2)/
     !           (75.6-0.1425*Temp    -2.38e-4*Temp    **2)
          AK(i)= (1.787-0.007*Temp    )/(1.+0.03225*Temp    )/
     !           ((1.787-0.007*TempN(i))/(1.+0.03225*TempN(i)))*
     !           (1.-7.37e-6*(TempN(i)-4.)**2+3.79e-8*(TempN(i)-4.)**3)/
     !           (1.-7.37e-6*(Temp    -4.)**2+3.79e-8*(Temp    -4.)**3)
        end if
        M=MatNum(i)
        if(Kappa(i).eq.-1) then
          hi1=amin1(hSat(M),hTemp(i)/Ah(i))
          hi2=amin1(hSat(M), hNew(i)/Ah(i))
        else if(Kappa(i).eq.+1) then
          hi1=amin1(hSat(M),hTemp(i)/Ah(i)/AhW(M))
          hi2=amin1(hSat(M), hNew(i)/Ah(i)/AhW(M))
        end if
        hiM=0.1*hi1+0.9*hi2
        if(hi1.ge.hSat(M).and.hi2.ge.hSat(M)) then
          Coni=ConSat(M)
        else if(hiM.gt.hTab(NTab).and.hiM.le.hTab(1)) then
          iT=int((alog10(-hiM)-alh1)/dlh)+1
          dh=(hiM-hTab(iT))/(hTab(iT+1)-hTab(iT))
          Coni=ConTab(iT,M)+(ConTab(iT+1,M)-ConTab(iT,M))*dh
        else
          Coni=FK(iModel,hiM,ParD(1,M))
        end if
        if(hiM.ge.hSat(M)) then
          Capi=0.
          Thei=ths(M)
        else if(hiM.ge.hTab(NTab).and.hiM.le.hTab(1)) then
          iT=int((alog10(-hiM)-alh1)/dlh)+1
          dh=(hiM-hTab(iT))/(hTab(iT+1)-hTab(iT))
          Capi=CapTab(iT,M)+(CapTab(iT+1,M)-CapTab(iT,M))*dh
          Thei=TheTab(iT,M)+(TheTab(iT+1,M)-TheTab(iT,M))*dh
        else
          Capi=FC(iModel,hiM,ParD(1,M))
          Thei=FQ(iModel,hiM,ParD(1,M))
        end if
        if(Kappa(i).eq.-1) then
          Con(i)=Coni*AK(i)*AKS(i)
          Cap(i)=Capi*ATh(i)*AThS(i)/Ah(i)
        else
          Con(i)=ConR(i)+Coni*AK(i)*AKS(i)*AKW(M)
          Cap(i)=Capi*ATh(i)*AThS(i)*AThW(M)/Ah(i)/AhW(M)
        end if
        if(Iter.eq.0) then
          ConO(i)=Con(i)
        end if
        if(Kappa(i).eq.-1) theta(i)=thr(M)+(Thei-thr(M))*ATh(i)*AThS(i)
        if(Kappa(i).eq.+1) theta(i)=ThRR(i)+AThW(M)*ATh(i)*AThS(i)*
     !                              (Thei-thr(M))
11    continue
      return
      end

***********************************************************************

      real function Fqh(GWL,Aqh,Bqh)
      Fqh=Aqh*exp(Bqh*abs(GWL))
      return
      end

************************************************************************

*     To calculate the velocities

      subroutine Veloc(N,hNew,Con,x,CosAlf,v,ThNewT,ThOldT,SinkTop,
     !                 ThNewB,ThOldB,SinkBot,dt)

      dimension hNew(N),x(N),Con(N),v(N)

      v(N)=-(Con(N)+Con(N-1))/2.*((hNew(N)-hNew(N-1))/(x(N)-x(N-1))+
     !     CosAlf)-(x(N)-x(N-1))/2.*((ThNewT-ThOldT)/dt+SinkTop)
      do 11 i=2,N-1
        dxA=x(i+1)-x(i)
        dxB=x(i)-x(i-1)
        vA=-(Con(i)+Con(i+1))/2.*((hNew(i+1)-hNew(i))/dxA+CosAlf)
        vB=-(Con(i)+Con(i-1))/2.*((hNew(i)-hNew(i-1))/dxB+CosAlf)
        v(i)=(vA*dxB+vB*dxA)/(dxA+dxB)
11    continue
      v(1)=-(Con(1)+Con(2))/2.*((hNew(2)-hNew(1))/(x(2)-x(1))+CosAlf)-
     !      (x(1)-x(2))/2.*((ThNewB-ThOldB)/dt+SinkBot)
      return
      end

************************************************************************

      subroutine Hyster(NumNP,NMat,hOld,MatNum,ParD,ParW,ThNew,ThOld,
     !                  Kappa,AThS,ThRR,ConO,ConR,AKS,KappaO,Ah,AK,
     !                  iHyst,iModel)

      real KsD,KsW,Ks,KW

      dimension MatNum(NumNP),ThOld(NumNP),hOld(NumNP),ParD(9,NMat),
     !          ThNew(NumNP),Kappa(NumNP),AThS(NumNP),ThRR(NUmNP),
     !          ConO(NumNP),ConR(NumNP),AKS(NUmNP),KappaO(NumNP),
     !          Ah(NumNP),AK(NumNP),ParW(9,NMat)

      do 11 i=1,NumNP

*     Check for reversal
      KappaO(i)=Kappa(i)
      if((ThNew(i)-ThOld(i))*Kappa(i).ge.-.00001) goto 11
      Kappa(i)=-Kappa(i)
      m=MatNum(i)
      Thr=ParD(1,m)
      ThsD=ParD(2,m)
      ThsW=ParW(2,m)
      KsD=ParD(5,m)
      KsW=ParW(5,m)

*     Update Ths and Ks for wetting scanning curve
      if(Kappa(i).eq.1) then
        if(ThsW.ge.0.999*ThsD) then
          Ths=ThsD
        else
          RR=1./(ThsD-ThsW)-1./(ThsD-Thr)
          Ths=ThsD-(ThsD-ThOld(i))/(1.+RR*(ThsD-ThOld(i)))
        end if
        if(KsW.ge.0.999*KsD) then
          Ks=KsD
        else
          RR=1./(KsD-KsW)-1./KsD
          Ks=KsD-(KsD-ConO(i))/(1.+RR*(KsD-ConO(i)))
        end if
      end if

*     Update parameters for scanning curve
      if(Kappa(i).eq.1) then ! Wetting
        AThS(i)=1.
        SeW=FS(iModel,hOld(i)/Ah(i),ParW(1,m))
        if(SeW.lt.0.999) AThS(i)=(ThOld(i)-Ths)/(1.-SeW)/(Thr-ThsW)
        ThRR(I)=Ths-AThS(i)*(ThsW-Thr)
        AKS(i)=1.
        ConR(i)=0.
        if(iHyst.eq.2) then
          KW=AK(i)*FK(iModel,hOld(i)/Ah(i),ParW(1,m))
          if(KW.lt.0.999*KsW) AKS(i)=(ConO(i)-Ks)/(KW-KsW)
          ConR(i)=Ks-AKS(i)*KsW
        end if
      else ! Drying
        AThS(i)=(ThOld(i)-Thr)/FS(iModel,hOld(i)/Ah(i),ParD(1,m))/
     !          (ThsD-Thr)
        ThRR(i)=Thr
        AKS(i)=1.
        ConR(i)=0.
        if(iHyst.eq.2)AKS(i)=ConO(i)/FK(iModel,hOld(i)/Ah(i),ParD(1,m))/
     !                       AK(i)
      end if
11    continue
      return
      end

* ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||