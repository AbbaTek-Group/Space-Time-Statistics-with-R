* Source file OUTPUT.FOR |||||||||||||||||||||||||||||||||||||||||||||||

      subroutine TLInf(Con1,Con2,ConN,ConB,x1,x2,xN,xB,CosAlf,t,dt,
     !                 IterW,IterC,TLevel,ShortO,TPrint,rTop,rRoot,
     !                 vRoot,hNewN,hNewB,hRoot,hNew1,hNew2,CumQ,ItCum,
     !                 KodTop,KodBot,ConvgF,lWat,lChem,cRoot,NS,NSD,
     !                 NumNP,Conc,cvTop,cvBot,cvCh0,cvCh1,Peclet,
     !                 Courant,wCumT,wCumA,cCumT,cCumA,CumCh,ThNewT,
     !                 ThOldT,SinkTop,ThNewB,ThOldB,SinkBot,lScreen,
     !                 ierr)
      integer TLevel
      logical ShortO,ConvgF,lWat,lChem,lScreen
      dimension CumQ(5),cRoot(NS),Conc(NSD,NumNP),cvTop(NS),cvBot(NS),
     !          cvCh0(NS),cvCh1(NS),cCumT(NS),cCumA(NS),CumCh(4,NS)

      vTop=-(ConN+ConB)/2*((hNewN-hNewB)/(xN-xB)+CosAlf)-
     !     (ThNewT-ThOldT)*(xN-xB)/2./dt-SinkTop*(xN-xB)/2.
      vBot=-(Con1+Con2)/2*((hNew2-hNew1)/(x2-x1)+CosAlf)-
     !     (ThNewB-ThOldB)*(x1-x2)/2./dt-SinkBot*(x1-x2)/2.
      CumQ(1)=CumQ(1)+rTop *dt
      CumQ(2)=CumQ(2)+rRoot*dt
      CumQ(3)=CumQ(3)+vTop *dt
      CumQ(4)=CumQ(4)+vRoot*dt
      CumQ(5)=CumQ(5)+vBot *dt
      wCumT=wCumT+(vBot-vTop-vRoot)*dt
      wCumA=wCumA+(abs(vBot)+abs(vTop)+abs(vRoot))*dt
      if(lChem) then
        do 11 jj=1,NS
          CumCh(1,jj)=CumCh(1,jj)-cvTop(jj)*dt
          CumCh(2,jj)=CumCh(2,jj)+cvBot(jj)*dt
          CumCh(3,jj)=CumCh(3,jj)+cvCh0(jj)*dt
          CumCh(4,jj)=CumCh(4,jj)+cvCh1(jj)*dt
          cCumT(jj)=cCumT(jj)+(cvTop(jj)-cvBot(jj)-cvCh0(jj)-cvCh1(jj))
     !              *dt
          cCumA(jj)=cCumA(jj)+(abs(cvBot(jj))+abs(cvTop(jj))+
     !                         abs(cvCh0(jj))+abs(cvCh1(jj)))*dt
11      continue
      end if
CYAP      if(lScreen) then
cYAP       if(float((TLevel+19)/20 ).eq.(TLevel+19)/20.)
cYAP++
CYAP     !  write(*,110)
CYAP        write(*,120) t,IterW,IterC,ItCum,vTop,CumQ(3),CumQ(4),CumQ(5),
CYAP     !               hNewN,hRoot,hNew1
cYAP--
cYAP      end if
      if(TLevel.eq.1) then
        write(71,130,err=901)
        if(lChem) then
          write(70,150,err=902)
          do 12 jj=1,NS
            write(80+jj,160,err=903)
12        continue
        else
          write(70,140,err=902)
        end if
      end if
      if(.not.ShortO.or.abs(TPrint-t).lt.0.001*dt) then
        if(lWat.or.TLevel.eq.1)
     !    write(71,170,err=901) t,rTop,rRoot,vTop,vRoot,vBot,
     !                         (CumQ(i),i=1,5),hNewN,hRoot,hNew1,TLevel
        if(lChem) then
          write(70,180,err=902) TLevel,t,dt,IterW,IterC,ItCum,KodTop,
     !                          KodBot,ConvgF,Peclet,Courant
          do 13 jj=1,NS
            write(80+jj,190,err=903) t,-cvTop(jj),cvBot(jj),CumCh(1,jj),
     !                       CumCh(2,jj),CumCh(3,jj),CumCh(4,jj),
     !                       Conc(jj,NumNP),cRoot(jj),Conc(jj,1),TLevel
13        continue
        else
          write(70,200,err=901) TLevel,t,dt,IterW,ItCum,KodTop,KodBot,
     !                          ConvgF
        end if
      end if
      return

*     Error when writing into an output file 
901   ierr=1
      return
902   ierr=2
      return
903   ierr=3
      return

110   format(/
     !'       Time ItW ItC ItCum  vTop    SvTop    SvRoot   SvBot    ',
     !'hTop hRoot hBot'/)
120   format(f12.3,2i3,i6,4e9.2,3f6.0)
130   format(/
     !'       Time        rTop      rRoot      vTop       vRoot      v',
     !'Bot     sum(rTop) sum(rRoot)  sum(vTop) sum(vRoot)',
     !'  sum(vBot)    hTop       hRoot      hBot      t-level'/
     !'        [T]       [L/T]      [L/T]      [L/T]      [L/T]   ',
     !'   [L/T]       [L]        [L]        [L]        [L]   ',
     !'     [L]        [L]        [L]        [L]'/)
140   format(//' TLevel   Time         dt      Iter  ItCum  KodT',
     !'  KodB  Convergency'/)
150   format(//' TLevel   Time         dt     ItrW ItrC  ItCum  KodT',
     !'  KodB Converg Peclet   Courant'/)
160   format(' All solute fluxes and cumulative solute fluxes are positi
     !ve out of the region'//
     !'       Time       cvTop      cvBot    Sum(cvTop) Sum(cvBo'
     !,'t)   cvCh0      cvCh1       cTop      cRoot       cBot    TLevel
     !'/'        [T]      [M/L2/T]   [M/L2/T]     [M]        [M] ',
     !'       [M]        [M]       [M/L3]    [M/L3]      [M/L3]')
170   format(f14.5,13e11.3,i8)
180   format(i5,e12.3,e12.3,2i5,3i6,l6,2f10.3)
190   format(f14.5,9e11.3,i8)
200   format(i5,e12.5,e12.3,i5,3i6,l6,2f10.3)
      end

************************************************************************

      subroutine ALInf(t,CumQ,hNewN,hRoot,hNew1,ALevel,ierr)

      integer ALevel
      dimension CumQ(5)

      if (ALevel.eq.1) write(72,110,err=901)
      write(72,120,err=901) t,(CumQ(i),i=1,5),hNewN,hRoot,hNew1,ALevel
      return

*     Error when writing into an output file 
901   ierr=1
      return

110   format(//
     !'   Time      sum(rTop)  sum(rRoot) sum(vTop)  sum(vRoot)',
     !' sum(vBot)    hTop       hRoot      hBot      A-level'/
     !'    [T]        [L]        [L]        [L]        [L]    ',
     !'    [L]        [L]         [L]       [L] '/)
120   format(f12.5,8e11.3,i8)
      end

************************************************************************

      subroutine SubReg(N,NMat,NLay,hNew,ThN,ThO,x,MatNum,LayNum,t,dt,
     !                  CosAlf,Con1,Con2,ConN,ConM,lChem,Conc,ChPar,
     !                  PLevel,ths,wCumT,wCumA,cCumT,cCumA,wVolI,cVolI,
     !                  WatIn,SolIn,lWat,lTemp,Temp,TPar,TDep,NS,NSD,
     !                  Sorb,lLinear,lEquil,ierr)

      logical lWat,lChem,lTemp,lLinear(NS),lEquil
      integer PLevel
      dimension hNew(N),ThN(N),ThO(N),x(N),MatNum(N),LayNum(N),
     !          Conc(NSD,N),ChPar(NSD*16+3,NMat),ths(NMat),cCumA(NS),
     !          cCumT(NS),cVolI(NS),WatIn(N),SolIn(N),Temp(N),
     !          TPar(10,NMat),TDep(NSD*16+3),Sorb(NSD,N),hMean(10),
     !          cMean(5,10),TMean(10),SubVol(10),SubCha(10),ConVol(5),
     !          ConSub(5,10),cTot(5),SubT(10),Area(10)

      ATot=0.
      Tr=293.15
      R=8.314
      if(lWat.or.PLevel.eq.0) then
        Volume=0.
        Change=0.
        hTot=0.
        DeltW=0.
      end if
      if(lTemp) then
        TTot=0.
        TVol=0.
      end if
      if(lChem) then
        do 11 jj=1,NS
          cTot(jj)=0.
          ConVol(jj)=0.
11      continue
        DeltC=0.
      end if
      do 13 Lay=1,NLay
        Area(Lay)=0.
        if(lWat.or.PLevel.eq.0) then
          SubVol(Lay)=0.
          SubCha(Lay)=0.
          hMean(Lay)=0.
        end if
        if(lTemp) then
          SubT(Lay)=0.
          TMean(Lay)=0.
        end if
        if(lChem) then
          do 12 jj=1,NS
            ConSub(jj,Lay)=0.
            cMean(jj,Lay)=0.
12        continue
        end if
13    continue

      do 15 i=N-1,1,-1
        cEl=0.
        Mi=MatNum(i)
        Mj=MatNum(i+1)
        Lay=LayNum(i)
        dx=x(i+1)-x(i)
        Area(Lay)=Area(Lay)+dx
        ATot=ATot+dx
        TT=(Temp(i)+Temp(i+1))/2.+273.15
        if(lWat.or.PLevel.eq.0) then
          hE=(hNew(i)+hNew(i+1))/2.
          VNewi=dx*(ThN(i)+ThN(i+1))/2.
          VOldi=dx*(ThO(i)+ThO(i+1))/2.
          Volume=Volume+VNewi
          Change=Change+(VNewi-VOldi)/dt
          SubCha(Lay)=SubCha(Lay)+(VNewi-VOldi)/dt
          SubVol(Lay)=SubVol(Lay)+VNewi
          hTot=hTot+hE*dx
          hMean(Lay)=hMean(Lay)+hE*dx
        end if
        if(lTemp) then
          TE=(Temp(i)+Temp(i+1))/2.
          TNewE=dx*((Temp(i)+273.15)*(TPar(2,Mi)*TPar(7,Mi)+
     !                    TPar(3,Mi)*TPar(8,Mi)+TPar(9,Mi)*ThN(i))+
     !              (Temp(i+1)+273.15)*(TPar(2,Mj)*TPar(7,Mj)+
     !                    TPar(3,Mj)*TPar(8,Mj)+TPar(9,Mj)*ThN(i+1)))/2.
          TVol=TVol+TNewE
          SubT(Lay)=SubT(Lay)+TNewE
          TTot=TTot+TE*dx
          TMean(Lay)=TMean(Lay)+TE*dx
        end if
        if(lChem) then
          do 14 jj=1,NS
            jjj=(jj-1)*16
            cE=(Conc(jj,i)+Conc(jj,i+1))/2.
            TTi=(Temp(i)-Tr)/R/TT/Tr
            xKsi=ChPar(jjj+6,Mi)*exp(TDep(jjj+6)*TTi)
            xNui=ChPar(jjj+7,Mi)*exp(TDep(jjj+7)*TTi)
            fExpi=ChPar(jjj+8,Mi)*exp(TDep(jjj+8)*TTi)
            Henryi=ChPar(jjj+9,Mi)*exp(TDep(jjj+9)*TTi)
            TTj=(Temp(i+1)-Tr)/R/TT/Tr
            xKsj=ChPar(jjj+6,Mj)*exp(TDep(jjj+6)*TTj)
            xNuj=ChPar(jjj+7,Mj)*exp(TDep(jjj+7)*TTj)
            fExpj=ChPar(jjj+8,Mj)*exp(TDep(jjj+8)*TTj)
            Henryj=ChPar(jjj+9,Mj)*exp(TDep(jjj+9)*TTj)
            C1=1.
            C2=1.
            if(.not.lLinear(jj)) then
              if(Conc(jj,i).gt.0.) C1=Conc(jj,i)**
     !             (fExpi-1.)/(1.+xNui*Conc(jj,i)**fExpi)
              if(Conc(jj,i+1).gt.0.) C2=Conc(jj,i+1)**
     !             (fExpj-1.)/(1.+xNuj*Conc(jj,i+1)**fExpj)
            end if
            cNewi=dx*
     !          (Conc(jj,i)*  (thN(i)+ChPar(3,Mi)*ChPar(1,Mi)*xKsi*C1+
     !                         (ths(Mi)-ThN(i))*Henryi)+
     !           Conc(jj,i+1)*(thN(i+1)+ChPar(3,Mj)*ChPar(1,Mj)*xKsj*C2+
     !                         (ths(Mj)-ThN(i+1))*Henryj))/2.
              if(.not.lEquil) cNewi=cNewi+
     !           dx*(ChPar(1,Mi)*Sorb(jj,i)+ChPar(1,Mj)*Sorb(jj,i+1))/2.
            ConVol(jj)=ConVol(jj)+cNewi
            ConSub(jj,Lay)=ConSub(jj,Lay)+cNewi
            cTot(jj)=cTot(jj)+cE*dx
            cMean(jj,Lay)=cMean(jj,Lay)+cE*dx
            if(jj.eq.1) cEl=cNewi
14        continue
        end if
        if(PLevel.eq.0) then
          if(lWat) WatIn(i)=vNewi
          if(lChem) SolIn(i)=cEl
        else
          if(lWat) DeltW=DeltW+abs(WatIn(i)-vNewi)
          if(lChem) DeltC=DeltC+abs(SolIn(i)-cEl)
        end if
15    continue
      do 17 Lay=1,NLay
        if(lWat.or.PLevel.eq.0) hMean(Lay)=hMean(Lay)/Area(Lay)
        if(lTemp) TMean(Lay)=TMean(Lay)/Area(Lay)
        do 16 jj=1,NS
          if(lChem) cMean(jj,Lay)=cMean(jj,Lay)/Area(Lay)
16      continue
17    continue
      if(lWat.or.PLevel.eq.0) hTot=hTot/ATot
      if(lTemp) TTot=TTot/ATot
      do 18 jj=1,NS
        if(lChem) cTot(jj)=cTot(jj)/ATot
18    continue
      v1=-(Con1+Con2)/2.*((hNew(2)-hNew(1))/(x(2)-x(1))+CosAlf)
      vN=-(ConN+ConM)/2.*((hNew(N)-hNew(N-1))/(x(N)-x(N-1))+CosAlf)

      write(76,110,err=901) t
      write(76,120,err=901) (i,i=1,NLay)
      write(76,130,err=901)
      write(76,140,err=901)   ATot,  (Area(i),i=1,NLay)
      if(lWat.or.PLevel.eq.0) then
        write(76,150,err=901) Volume,(SubVol(i),i=1,NLay)
        write(76,160,err=901) Change,(SubCha(i),i=1,NLay)
        write(76,170,err=901) hTot,  ( hMean(i),i=1,NLay)
      end if
      if(lTemp) then
        write(76,180,err=901) TVol,  (  SubT(i),i=1,NLay)
        write(76,190,err=901) TTot,  ( TMean(i),i=1,NLay)
      end if
      if(lChem) then
        do 19 jj=1,NS
          write(76,200,err=901) jj,ConVol(jj),(ConSub(jj,i),i=1,NLay)
          write(76,210,err=901) jj,cTot(jj),  ( cMean(jj,i),i=1,NLay)
19      continue
      end if
      if(lWat.or.PLevel.eq.0) write(76,220,err=901) vN,v1

*     Mass balance calculation
      if(PLevel.eq.0) then
        wVolI=Volume
        do 20 jj=1,NS
          cVolI(jj)=ConVol(jj)
20      continue
      else
        if(lWat) then
          wBalT=Volume-wVolI-wCumT
          write(76,230,err=901) wBalT
          ww=amax1(DeltW,wCumA)
          if(ww.gt.1.e-25) then
            wBalR=abs(wBalT)/ww*100.
            write(76,240,err=901) wBalR
          end if
        end if
        if(lChem) then
          do 21 jj=1,NS
            cBalT=ConVol(jj)-cVolI(jj)+cCumT(jj)
            write(76,250,err=901) jj,cBalT
            cc=amax1(DeltC,cCumA(jj))
            if(cc.gt.1.e-25) then
              cBalR=abs(cBalT)/cc*100.
              write(76,260,err=901) jj,cBalR
            end if
21        continue
        end if
      end if
      write(76,130,err=901)
      return

*     Error when writing into an output file 
901   ierr=1
      return

110   format(/'-----------------------------------------------------'/
     !        ' Time       [T]',f12.4/
     !        '-----------------------------------------------------')
120   format( ' Sub-region num.             ',9(I7,4x))
130   format( '-----------------------------------------------------')
140   format( ' Area     [L]      ',e11.3,9e11.3)
150   format( ' W-volume [L]      ',e11.3,9e11.3)
160   format( ' In-flow  [L/T]    ',e11.3,9e11.3)
170   format( ' h Mean   [L]      ',e11.3,9e11.3)
180   format( ' HeatVol  [M/T2]   ',e11.3,10e11.3)
190   format( ' tMean    [K]      ',f11.3,10f11.3)
200   format( ' ConcVol  [M/L2] ',i1,1x,e11.3,10e11.3)
210   format( ' cMean    [M/L3] ',i1,1x,e11.3,10e11.3)
220   format( ' Top Flux [L/T]    ',e11.3/
     !        ' Bot Flux [L/T]    ',e11.3)
230   format( ' WatBalT  [L]      ',e11.3)
240   format( ' WatBalR  [%]      ',f11.3)
250   format( ' CncBalT  [M]    ',i1,1x,e11.3)
260   format( ' CncBalR  [%]    ',i1,1x,f11.3)
      end

***********************************************************************

      subroutine NodOut(N,NMat,hNew,thN,Con,x,xSurf,CosAlf,TPrint,
     !                  MatNum,Cap,Bxz,Sink,ConS,NS,NSD,Conc,Temp,Sorb,
     !                  ierr)

      dimension hNew(N),thN(N),Con(N),x(N),MatNum(N),Cap(N),Bxz(N),
     !          Sink(N),ConS(NMat),Conc(NSD,N),Temp(N),Sorb(NSD,N)

      M=MatNum(N)
      ConSN=ConS(M)*Bxz(N)
      vN=-(Con(N)+Con(N-1))/2*((hNew(N)-hNew(N-1))/(x(N)-x(N-1))+CosAlf)
      write(75,110,err=901) TPrint
      write(75,120,err=901) 1,x(N)-xSurf,hNew(N),thN(N),Con(N),Cap(N),
     !                      vN,Sink(N),Con(N)/ConSN,vN/ConSN,Temp(N),
     !                      (Conc(jj,N),jj=1,NS),(Sorb(jj,N),jj=1,NS)
      do 11 i=N-1,2,-1
        dxA=x(i+1)-x(i)
        dxB=x(i)-x(i-1)
        vA=-(Con(i)+Con(i+1))/2.*((hNew(i+1)-hNew(i))/dxA+CosAlf)
        vB=-(Con(i)+Con(i-1))/2.*((hNew(i)-hNew(i-1))/dxB+CosAlf)
        vi= (vA*dxA+vB*dxB)/(dxA+dxB)
        write(75,120,err=901) N-i+1,x(i)-xSurf,hNew(i),thN(i),Con(i),
     !                        Cap(i),vi,Sink(i),Con(i)/ConSN,vi/ConSN,
     !                        Temp(i),(Conc(jj,i),jj=1,NS),
     !                        (Sorb(jj,i),jj=1,NS)
11    continue
      M=MatNum(1)
      v1=-(Con(1)+Con(2))/2.*((hNew(2)-hNew(1))/(x(2)-x(1))+CosAlf)
      write(75,120,err=901) N,x(1)-xSurf,hNew(1),thN(1),Con(1),Cap(1),
     !                      v1,Sink(1),Con(1)/ConSN,v1/ConSN,Temp(1),
     !                      (Conc(jj,1),jj=1,NS),(Sorb(jj,1),jj=1,NS)
      write(75,'(''end'')',err=901)
      return

*     Error when writing into an output file 
901   ierr=1
      return

110   format(//' Time:',f13.6//
     !' Node    Depth      Head Moisture       K          C         Fl',
     !'ux        Sink        K/KsTop  v/KsTop   Temp      Conc(1...NS)'/
     !'         [L]       [L]     [-]        [L/T]      [1/L]      [L/',
     !'T]        [1/T]         [-]      [-]      [C]        [M/L*3]'/)
120   format(i4,1x,f8.2,1x,f10.2,1x,f6.3,2x,4e12.3,2f10.3,f8.2,20e12.3)
      end

***********************************************************************

      subroutine ObsNod(t,N,NObs,NS,NSD,Node,Conc,hNew,ThNew,TempN,
     !                  lChem,ierr)

      dimension Node(NObs),Conc(NSD,N),ThNew(N),TempN(N),hNew(N)
      logical lChem
CYAP++
      common /ptf/iptf
CYAP--
      if(.not.lChem) then
        write(77,100,err=901) t,(hNew(Node(i)),ThNew(Node(i)),
     !                        TempN(Node(i)),i=1,NObs)
CYAP++
	  write(600+iptf,'(F12.3,10F12.4)') t,(ThNew(Node(i)),i=1,NObs)
CYAP--
      else
        if(NS.eq.1)write(77,110,err=901)t,(hNew(Node(i)),ThNew(Node(i)),
     !                 TempN(Node(i)),(Conc(j,Node(i)),j=1,NS),i=1,NObs)
        if(NS.eq.2)write(77,120,err=901)t,(hNew(Node(i)),ThNew(Node(i)),
     !                 TempN(Node(i)),(Conc(j,Node(i)),j=1,NS),i=1,NObs)
        if(NS.eq.3)write(77,130,err=901)t,(hNew(Node(i)),ThNew(Node(i)),
     !                 TempN(Node(i)),(Conc(j,Node(i)),j=1,NS),i=1,NObs)
        if(NS.eq.4)write(77,140,err=901)t,(hNew(Node(i)),ThNew(Node(i)),
     !                 TempN(Node(i)),(Conc(j,Node(i)),j=1,NS),i=1,NObs)
        if(NS.eq.5)write(77,150,err=901)t,(hNew(Node(i)),ThNew(Node(i)),
     !                 TempN(Node(i)),(Conc(j,Node(i)),j=1,NS),i=1,NObs)
      end if

100   format(2x,f12.3,10(f12.3,f8.4,f9.3,       2x))
110   format(2x,f12.3,10(f12.3,f8.4,f9.3, e11.3,2x))
120   format(2x,f12.3,10(f12.3,f8.4,f9.3,2e11.3,2x))
130   format(2x,f12.3,10(f12.3,f8.4,f9.3,3e11.3,2x))
140   format(2x,f12.3,10(f12.3,f8.4,f9.3,4e11.3,2x))
150   format(2x,f12.3,10(f12.3,f8.4,f9.3,5e11.3,2x))
      return

*     Error when writing into an output file 
901   ierr=1
      return

      end

* ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||