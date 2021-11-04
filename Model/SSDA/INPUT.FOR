* Source file INPUT.FOR ||||||||||||||||||||||||||||||||||||||||||||||||

      subroutine BasInf(CosAlf,MaxIt,TolTh,TolH,TopInF,BotInF,ShortO,
     !                  lWat,lChem,SinkF,WLayer,qGWLF,FreeD,SeepF,AtmBC,
     !                  KodTop,KodBot,rTop,rRoot,rBot,hCritS,hCritA,
     !                  GWL0L,Aqh,Bqh,kTOld,kBOld,NUnitD,iUnit,NMat,
     !                  NMatD,NLay,lRoot,lTemp,lWDep,lEquil,lScreen,
     !                  ierr)

      character*72 Hed
      character*5  LUnit,TUnit,MUnit
      logical TopInF,BotInF,ShortO,lWat,lChem,lTemp,SinkF,WLayer,qGWLF,
     !        FreeD,SeepF,AtmBC,lRoot,lWDep,lEquil,lScreen
      dimension iUnit(NUnitD)

      read(30,*,err=901)
      read(30,*,err=901)
      read(30,'(a)',err=901) Hed
      read(30,*,err=901)
      read(30,'(a)',err=901) LUnit
      read(30,'(a)',err=901) TUnit
      read(30,'(a)',err=901) MUnit
      read(30,*,err=901)
      read(30,*,err=901) lWat,lChem,lTemp,SinkF,lRoot,ShortO,lWDep,
     !                   lScreen,AtmBC,lEquil
      read(30,*,err=901)
      read(30,*,err=901) NMat,NLay,CosAlf
      if(NMat.gt.NMatD) then
        ierr=4
        return
      end if
      read(30,*,err=902)
      read(30,*,err=902)
      read(30,*,err=902) MaxIt,TolTh,TolH
      read(30,*,err=902)
      read(30,*,err=902) TopInF,WLayer,KodTop
      read(30,*,err=902)
      read(30,*,err=902) BotInF,qGWLF,FreeD,SeepF,KodBot
      if(.not.TopInF.and..not.BotInF.and.
     !                              (KodTop.eq.-1.or.KodBot.eq.-1)) then
        read(30,*,err=902)
        read(30,*,err=902) rTop,rBot,rRoot  !,hCritS,hCritA
      end if
      if(qGWLF) then
        read(30,*,err=902)
        read(30,*,err=902) GWL0L,Aqh,Bqh
      end if

*     Input modifications
      rRoot=abs(rRoot)
      hCritA=-abs(hCritA)
      if(TopInF) KodTop=isign(3,KodTop)
      if(BotInF) KodBot=isign(3,KodBot)
      if(AtmBC) then
        hCritS=0
        KodTop=-4
      end if
      if(WLayer) KodTop=-iabs(KodTop)
      if(qGWLF)  KodBot=-7
      if(FreeD)  KodBot=-5
      if(SeepF)  KodBot=-2
      kTOld=KodTop
      kBOld=KodBot

      if(lScreen) then
c        write(*,*)'----------------------------------------------------'
c        write(*,*)'|                                                  |'
c        write(*,*)'|                    HYDRUS                        |'
c        write(*,*)'|                                                  |'
c        write(*,*)'|   Code for simulating one-dimensional variably   |'
c        write(*,*)'|    saturated water flow, heat transport, and     |'
c        write(*,*)'|   transport of solutes involved in sequential    |'
c        write(*,*)'|         first-order decay reactions              |'
c        write(*,*)'|                                                  |'
c        write(*,*)'|                  version 6.0                     |'
c        write(*,*)'|                                                  |'
c        write(*,*)'|          Last modified: September, 1996          |'
c        write(*,*)'|                                                  |'
c        write(*,*)'----------------------------------------------------'
c        write(*,*)
        write(*,*) Hed
        write(*,*)
      end if
      do 11 i=1,NUnitD
        write(iUnit(i),*,err=903)'******* Program HYDRUS'
        write(iUnit(i),*,err=903)'******* ',Hed
C GM Uncomment next two lines
        call getdat(ii,imonth,iday)
        call gettim(ihours,mins,isecs,ii)
        write(iUnit(i),100,err=903) iday,imonth,ihours,mins,isecs
        write(iUnit(i),*,err=903)'Units: L = ',LUnit,', T = ',TUnit,
     !                           ', M = ',MUnit
11    continue

      write(50,*,err=903)
      write(50,*,err=903) 'CosAlf,MaxIt,TolTh,  TolH'
      write(50,110,err=903) CosAlf,MaxIt,TolTh,TolH
      write(50,*,err=903)
      write(50,*,err=903) 'TopInF,BotInF,AtmBC,SinkF,WLayer,qGWLF,FreeD,
     !SeepF,lWat,lChem,lTemp,lRoot,lWDep'
      write(50,120,err=903) TopInF,BotInF,AtmBC,SinkF,WLayer,qGWLF,
     !                      FreeD,SeepF,lWat,lChem,lTemp,lRoot,lWDep
      if(.not.(TopInF.or.BotInF)) write(72,130,err=903)
      return

*     Error when reading from an input file
901   close(30)
      ierr=1
      return
902   close(30)
      ierr=2
      return
*     Error when writing into an output file
903   ierr=3
      return

100   format(' Date: ',i3,'.',i2,'.','    Time: ',i3,':',i2,':',i2)
110   format(f6.3,i5,f8.3,f8.5)
120   format(13l6)
130   format(////)
      end

************************************************************************

      subroutine NodInf(NumNPD,NumNP,NObsD,NObs,hTop,hBot,x,hNew,hOld,
     !                  MatNum,hTemp,LayNum,Beta,Ah,AK,ATh,Conc,Sorb,
     !                  TempN,TempO,Node,NSD,NS,xSurf,lChem,lTemp,
     !                  lEquil,lScreen,ierr)

      character*30 Text1,Text2,Text3
      dimension x(NumNPD),hNew(NumNPD),hOld(NumNPD),MatNum(NumNPD),
     !          hTemp(NumNPD),LayNum(NumNPD),Beta(NumNPD),Ah(NumNPD),
     !          AK(NumNPD),ATh(NumNPD),Conc(NSD,NumNPD),TempN(NumNPD),
     !          TempO(NumNPD),Node(NObsD),Sorb(NSD,NumNPD),SConc(5),
     !          SSorb(5),C(5),S(5)
      logical lChem,lTemp,lEquil,lScreen

      if(lScreen) write(*,*)'reading nodal information'
      read(32,*,err=901) n
      do 11 i=1,n
        read(32,*,err=901)
11    continue
      read(32,*,err=901) NumNP,NS
      if(NumNP.gt.NumNPD) then
        ierr=3
        return
      end if
      if(NS.gt.NSD) then
        ierr=5
        return
      end if

*     Read nodal point information
      j=NumNP+1
12    continue
      j=j-1
      if(.not.lChem.and..not.lTemp) then
        read(32,*,err=901) n,x1,h,M,L,B,Ax,Bx,Dx
      else if(.not.lChem) then
        read(32,*,err=901) n,x1,h,M,L,B,Ax,Bx,Dx,Te
      else if(lEquil) then
        read(32,*,err=901) n,x1,h,M,L,B,Ax,Bx,Dx,Te,(C(ii),ii=1,NS)
      else
        read(32,*,err=901) n,x1,h,M,L,B,Ax,Bx,Dx,Te,(C(ii),ii=1,NS),
     !             (S(ii),ii=1,NS)
      end if
      n=NumNP-n+1
      x(n)=x1
      hOld(n)=h
      MatNum(n)=M
      LayNum(n)=L
      Beta(n)=B
      Ah(n)=Ax
      AK(n)=Bx
      ATh(n)=Dx
      TempO(n)=Te
      do 1 ii=1,NS
        Conc(ii,n)=C(ii)
        Sorb(ii,n)=S(ii)
1     continue

      if(j-n) 13,18,14
13    write(*,*)'ERROR in NodInf at node =', n
      stop
14    continue
      dx=x(nOld)-x(n)
      ShOld=(hOld(nOld)-hOld(n))/dx
      SBeta=(Beta(nOld)-Beta(n))/dx
      SAh=(Ah(nOld)-Ah(n))/dx
      SAK=(AK(nOld)-AK(n))/dx
      SATh=(ATh(nOld)-ATh(n))/dx
      if(lTemp) STemp=(TempO(nOld)-TempO(n))/dx
      if(lChem) then
        do 15 ii=1,NS
          SConc(ii)=(Conc(ii,nOld)-Conc(ii,n))/dx
          SSorb(ii)=(Sorb(ii,nOld)-Sorb(ii,n))/dx
15      continue
      end if
      do 17 i=nOld-1,n+1,-1
        dx=x(nOld)-x(i)
        hOld(i)=hOld(nOld)-ShOld*dx
        Beta(i)=Beta(nOld)-SBeta*dx
        Ah(i)=Ah(nOld)-SAh*dx
        AK(i)=AK(nOld)-SAK*dx
        ATh(i)=ATh(nOld)-SATh*dx
        if(lTemp) TempO(i)=TempO(nOld)-STemp*dx
        if(lChem) then
          do 16 ii=1,NS
            Conc(ii,i)=Conc(ii,nOld)-SConc(ii)*dx
            Sorb(ii,i)=Sorb(ii,nOld)-SSorb(ii)*dx
16        continue
        end if
        MatNum(i)=MatNum(i+1)
        LayNum(i)=LayNum(i+1)
17    continue
      j=n
18    continue
      nOld=n
      if(j.gt.1) goto 12

      SBeta=Beta(NumNP)*(x(NumNP)-x(NumNP-1))/2.
      do 19 i=2,NumNP-1
        SBeta=SBeta+Beta(i)*(x(i+1)-x(i-1))/2.
19    continue
      do 20 i=2,NumNP
        if(SBeta.gt.0.) then
          Beta(i)=Beta(i)/SBeta
        else
          Beta(i)=0.
        end if
20    continue
      xSurf=x(NumNP)

*     Print nodal information
      write(50,110,err=902)
      do 21 n=NumNP,1,-1
        if(.not.lChem.and..not.lTemp) then
          write(50,120,err=902) NumNP-n+1,x(n),hOld(n),MatNum(n),
     !                          LayNum(n),Beta(n),Ah(n),AK(n),ATh(n)
        else if(.not.lChem) then
          write(50,120,err=902) NumNP-n+1,x(n),hOld(n),MatNum(n),
     !                          LayNum(n),Beta(n),Ah(n),AK(n),ATh(n),
     !                          TempO(n)
        else if(lEquil) then
          write(50,120,err=902) NumNP-n+1,x(n),hOld(n),MatNum(n),
     !                          LayNum(n),Beta(n),Ah(n),AK(n),ATh(n),
     !                          TempO(n),(Conc(ii,n),ii=1,NS)
        else
          write(50,120,err=902) NumNP-n+1,x(n),hOld(n),MatNum(n),
     !                          LayNum(n),Beta(n),Ah(n),AK(n),ATh(n),
     !                          TempO(n),(Conc(ii,n),ii=1,NS),
     !                          (Sorb(ii,n),ii=1,NS)
        end if
        hNew(n) =hOld(n)
        hTemp(n)=hOld(n)
        TempN(n)=TempO(n)
21    continue
      write(50,'(''end'')',err=902)
      hBot=hNew(1)
      hTop=hNew(NumNP)
      write(50,130,err=902) NS

      read(32,*,err=901) NObs
      if(NObs.gt.NObsD) then
        ierr=4
        return
      end if
      if(Nobs.gt.0) then
        read(32,*,err=901) (Node(i),i=1,NObs)
        do 22 i=1,NObs
          Node(i)=NumNP-Node(i)+1
22      continue
        Text1='    h        theta    Temp   '
        Text2='  Conc     '
        Text3='Node('
        if (.not.lChem) then
          write(77,240,err=902) (Text3,NumNP-Node(j)+1,j=1,NObs)
          write(77,250,err=902) (Text1,i=1,NObs)
        else
          if(NS.eq.1) write(77,140,err=902)
     !      (Text3,NumNP-Node(j)+1,j=1,NObs)
          if(NS.eq.2) write(77,150,err=902)
     !      (Text3,NumNP-Node(j)+1,j=1,NObs)
          if(NS.eq.3) write(77,160,err=902)
     !      (Text3,NumNP-Node(j)+1,j=1,NObs)
          if(NS.eq.4) write(77,170,err=902)
     !      (Text3,NumNP-Node(j)+1,j=1,NObs)
          if(NS.eq.5) write(77,180,err=902)
     !      (Text3,NumNP-Node(j)+1,j=1,NObs)
         if(NS.eq.1)write(77,190,err=902)(Text1,(Text2,j=1,NS),i=1,NObs)
         if(NS.eq.2)write(77,200,err=902)(Text1,(Text2,j=1,NS),i=1,NObs)
         if(NS.eq.3)write(77,210,err=902)(Text1,(Text2,j=1,NS),i=1,NObs)
         if(NS.eq.4)write(77,220,err=902)(Text1,(Text2,j=1,NS),i=1,NObs)
         if(NS.eq.5)write(77,230,err=902)(Text1,(Text2,j=1,NS),i=1,NObs)
        end if
      end if
      return

*     Error when reading from an input file
901   ierr=1
      return
*     Error when writing into an output file
902   ierr=2
      return

110   format (/'Nodal point information'//
     !'node      x         hOld    MatN LayN  Beta      Ah       AK ',
     !'     ATh     Temp    Conc(1...NS)         Sorb(1...NS)'/)
120   format (i4,2f11.3,2i5,f8.3,3f9.3,f8.2,10e12.4,10e12.4)
130   format (/' Number of species in the chain : ',i3)
240   format (///14x,10(15x,a5,i3,')', 9x))
140   format (///14x,10(15x,a5,i3,')',18x))
150   format (///14x,10(15x,a5,i3,')',29x))
160   format (///14x,10(15x,a5,i3,')',40x))
170   format (///14x,10(15x,a5,i3,')',51x))
180   format (///14x,10(15x,a5,i3,')',62x))
250   format (/'       time     ',10(a29,     2x))
190   format (/'       time     ',10(a29, a11,2x))
200   format (/'       time     ',10(a29,2a11,2x))
210   format (/'       time     ',10(a29,3a11,2x))
220   format (/'       time     ',10(a29,4a11,2x))
230   format (/'       time     ',10(a29,5a11,2x))
      end

************************************************************************

      subroutine MatIn(NMat,ParD,ParW,hTab1,hTabN,lScreen,ierr,NumNP,Ah,
     !                 iHyst,AhW,AThW,AKW,MatNum,hNew,Kappa,AThS,ThRR,
     !                 ConR,AKS,KappaO,iModel)

      dimension ParD(9,NMat),ParW(9,NMat),Ah(NumNP),AhW(NMat),AKW(NMat),
     !          AThW(NMat),MatNum(NumNP),hNew(NumNP),Kappa(NumNP),
     !          AThS(NumNP),ThRR(NumNP),ConR(NumNP),AKS(NumNP),
     !          KappaO(NumNP)
      logical lScreen

      if(lScreen) write(*,*)'reading material information'
      read(30,*,err=901)
      read(30,*,err=901) hTab1,hTabN
      hTab1=-amin1(abs(hTab1),abs(hTabN))
      hTabN=-amax1(abs(hTab1),abs(hTabN))
      read(30,*,err=901)
      read(30,*,err=901) iModel,iHyst
      if(iHyst.gt.0) then
        read(30,*,err=901)
        read(30,*,err=901) IKappa
      else
        IKappa=-1
      endif
      do 11 i=1,NumNP
        Kappa(i)=IKappa
        KappaO(i)=IKappa
11    continue
      if(iModel.eq.2.or.(iModel.eq.0.and.iHyst.eq.0)) then
        write(50,110,err=902)
      else if(iModel.eq.1) then
        write(50,111,err=902)
      else
        write(50,112,err=902)
      end if
      read(30,*,err=901)
      do 12 M=1,NMat
        if(iHyst.eq.0) then
          if(iModel.eq.1) then
            NPar=9
          else
            NPar=5
          end if
          read(30,*,err=901) (ParD(i,M),i=1,NPar)
         if(iModel.eq.1) then
            ParD(6,M)=amax1(ParD(6,M),ParD(2,M))
            ParD(7,M)=amin1(ParD(7,M),ParD(1,M))
          else if(iModel.eq.0) then
            ParD(6,M)=ParD(2,M)
            ParD(7,M)=ParD(1,M)
            ParD(8,M)=ParD(2,M)
            ParD(9,M)=ParD(5,M)
          end if
          write(50,120,err=902) M,(ParD(i,M),i=1,NPar)
        else
          read(30,*) (ParD(i,M),i=1,6),ParW(2,M),ParW(3,M),ParW(5,M)
          ParD(6,M)=amax1(ParD(6,M),ParD(2,M))
          write(50,120,err=902) M,(ParD(i,M),i=1,6),ParW(2,M),
     !                          ParW(3,M),ParW(5,M)
          ParW(1,M)=ParD(1,M)
          ParW(4,M)=ParD(4,M)
          AhW(M)=ParD(3,M)/ParW(3,M)
          AThW(M)=(ParW(2,M)-ParW(1,M))/(ParD(2,M)-ParD(1,M))
          AKW(M)=1.0
          if(iHyst.eq.2) AKW(M)=ParW(5,M)/ParD(5,M)
          ParW(6,M)=ParW(1,M)+AThW(M)*(ParD(6,M)-ParD(1,M))
          ParD(7,M)=ParD(1,M)
          ParD(8,M)=ParD(2,M)
          ParD(9,M)=ParD(5,M)
          ParW(7,M)=ParW(1,M)
          ParW(8,M)=ParD(2,M)
          ParW(9,M)=ParW(5,M)
        end if
12    continue

*     Hysteresis Update for Initial Pressure Head Distributions

      do 13 i=1,NumNP
        M=MatNum(i)
        hNew(i)=amax1(hNew(i),Ah(i)*FH(iModel,0.00000001,ParD(1,M)))
        AThS(i)=1.
        AKS(i)=1.
        ThRR(i)=ParD(1,M)
        ConR(i)=0.
13    continue
      return

*     Error when reading from an input file
901   ierr=1
      return
*     Error when writing into an output file
902   ierr=2
      return

110   format(//'MatNum, Param. array:'//'   Mat     Qr     Qs        ',
     !'Alfa         n          Ks'/)
111   format(//'MatNum, Param. array:'//'   Mat     Qr     Qs        ',
     !'Alfa         n          Ks       Qm     Qa     Qk       Kk'/)
112   format(//'MatNum, Param. array:'//'   Mat     Qr     Qs        ',
     !'Alfa         n          Ks       Qm     QsW  AlfaW      KsW'/)
120   format(i5,2x,2f7.3,3e12.3,3f7.3,e12.3)
      end

************************************************************************

      subroutine GenMat(NTab,NTabD,NMat,thr,ths,hSat,Par,hTab,ConTab,
     !                  CapTab,ConSat,TheTab,iModel,lScreen,ierr)

      dimension hTab(NTab),ConTab(NTabD,NMat),CapTab(NTabD,NMat),
     !          Par(9,NMat),ConSat(NMat),thr(NMat),hSat(NMat),
     !          ths(NMat),TheTab(NTabD,NMat)
      logical lScreen

      if(lScreen) write(*,*)'generating materials'
      write(50,110,err=901)
      write(50,120,err=901)
      hTab1=hTab(1)
      hTabN=hTab(NTab)
      dlh=(alog10(-hTabN)-alog10(-hTab1))/(NTab-1)
      do 11 i=1,NTab
        alh=alog10(-hTab1)+(i-1)*dlh
        hTab(i)=-10**alh
11    continue
      do 13 M=1,NMat
        hSat(M)  =FH(iModel,1.0,Par(1,M))
        ConSat(M)=Par(5,M)
        thr(M)   =Par(1,M)
        ths(M)   =Par(2,M)
        do 12 i=1,NTab
          ConTab(i,M)=FK(iModel,hTab(i),Par(1,M))
          CapTab(i,M)=FC(iModel,hTab(i),Par(1,M))
          TheTab(i,M)=FQ(iModel,hTab(i),Par(1,M))
          Qe         =FS(iModel,hTab(i),Par(1,M))
          a10h=alog10(max(-hTab(i),1e-30))
          a10K=alog10(ConTab(i,M))
          write(50,130,err=901) TheTab(i,M),hTab(i),a10h,CapTab(i,M),
     !                          ConTab(i,M),a10K,Qe
12      continue
        write(50,140,err=901)
13    continue
      return

*     Error when writing into an output file
901   ierr=1
      return

110   format(/7x,'Table of Hydraulic Properties which are interpolated i
     !n simulation'/7x,65('=')/)
120   format('  theta         h        log h        C             K',
     !'        log K          S')
130   format(f8.4,e12.3,e12.4,e12.4,e12.4,e12.4,f10.4)
140   format('end')
      end

************************************************************************

      subroutine TmIn(tInit,tMax,tAtm,tOld,dt,dtMax,DMul,DMul2,dtMin,
     !                TPrint,t,dtOpt,TopInF,BotInF,lScreen,ItMin,
     !                ItMax,MaxAL,hCritS,ierr)

      logical TopInF,BotInF,lScreen
      dimension TPrint(50)

      if(lScreen) write(*,*) 'reading time information'
      read(30,*,err=901)
      read(30,*,err=901)
      read(30,*,err=901) dt,dtMin,dtMax,DMul,DMul2,ItMin,ItMax,MPL
      if(MPL.gt.50) then
        ierr=2
        return
      end if
      read(30,*,err=901)
      read(30,*,err=901) tInit,tMax
      read(30,*,err=901)
      read(30,*,err=901) (TPrint(i),i=1,MPL)
      if(TopInF.or.BotInF) then
        read(31,*,err=901) 
        read(31,*,err=901)
        read(31,*,err=901) MaxAL
        read(31,*,err=901)
        read(31,*,err=901) hCritS
        read(31,*,err=901)
      else
        tAtm=tMax
      end if
      TPrint(MPL+1)=tMax
      tOld=tInit
      t=tInit+dt
      return

*     Error when reading from an input file
901   ierr=1
      return
      end

************************************************************************

      subroutine SinkIn(NMat,lChem,lMoSink,lSolRed,lSolAdd,P0,POptm,
     !                  P2H,P2L,P3,r2H,r2L,aOsm,c50,P3c,NS,lScreen,ierr)

      dimension POptm(NMat),aOsm(NS)
      logical lChem,lMoSink,lSolRed,lSolAdd,lScreen

      if(lScreen) write(*,*) 'reading sink information'
      read(30,*,err=901)
      read(30,*,err=901)
      read(30,*,err=901) iMoSink
      if(iMoSink.eq.0) then
        lMoSink=.true.
      else
        lMoSink=.false.
      end if
      read(30,*,err=901)
      if(lMoSink) then
        read(30,*,err=901) P0,P2H,P2L,P3,r2H,r2L
        read(30,*,err=901)
        read(30,*,err=901) (POptm(i),i=1,NMat)
        P0 =-abs(P0)
        P2L=-abs(P2L)
        P2H=-abs(P2H)
        P3 =-abs(P3)
      else
        read(30,*,err=901) P0,P3
      end if
      if(lChem) then
        read(30,*,err=901)
        read(30,*,err=901) lSolRed
       if(lSolRed) then
          read(30,*,err=901)
          read(30,*,err=901) lSolAdd
          read(30,*,err=901)
          if(lSolAdd) then
            read(30,*,err=901) (aOsm(i),i=1,NS)
         else
            read(30,*,err=901) c50,P3c,(aOsm(i),i=1,NS)
          end if
        end if
      end if
      
      return

*     Error when reading from an input file
901   ierr=1
      return
      end

************************************************************************

      subroutine RootIn(tRMin,tRHarv,xRMin,xRMax,RGR,lScreen,ierr)

      logical lScreen

      if(lScreen) write(*,*) 'reading of root growth information'
      read(30,*,err=901)
      read(30,*,err=901)
      read(30,*,err=901) iRFak,tRMin,tRMed,tRHarv,xRMin,xRMed,xRMax
      if(iRFak.eq.1) then
        tRMed=(tRHarv+tRMin)/2.
        xRMed=(xRMax+xRMin)/2.
      end if
      rtm=tRMed-tRMin
      RGR=-(1./rtm)*alog(amax1(.0001,(xRMin*(xRMax-xRMed)))/
     !                  (xRMed*(xRMax-xRMin)))
      write(50,110,err=902) tRMin,tRHarv,xRMin,xRMax,RGR
      return

*     Error when reading from an input file
901   ierr=1
      return
*     Error when writing into an output file
902   ierr=2
      return

110   format(/' Root growth information'/1x,23('=')/' tRMin = ',f10.3,
     !' tRHarv = ',f10.3/' xRMin = ',f10.3,' xRMax = ',f10.3/
     !' Root growth rate = ',e11.3)
      end

************************************************************************

      subroutine TempIn(NMat,TPar,Ampl,tPeriod,kTopT,tTop,kBotT,tBot,
     !                  TopInf,BotInf,lScreen,ierr)

      logical TopInf,BotInf,lScreen
      dimension TPar(10,NMat)

      if(lScreen) write(*,*) 'reading heat transport information'
      read(30,*,err=901)
      read(30,*,err=901)
      do 11 i=1,NMat
        read(30,*,err=901)    (TPar(j,i),j=1,9)
11    continue
      read(30,*,err=901)
      read(30,*,err=901) Ampl,tPeriod
      write(50,110,err=902) Ampl
      do 12 i=1,NMat
        write(50,120,err=902) (TPar(j,i),j=1,9)
12    continue
      read(30,*,err=901)
      read(30,*,err=901) kTopT,tT,kBotT,tB
      if(.not.TopInf) tTop=tT
      if(.not.BotInf) tBot=tB
      return

*     Error when reading from an input file
901   ierr=1
      return
*     Error when writing into an output file
902   ierr=2
      return

110   format(//' Heat transport information'/1x,26('=')//' ample = ',
     ! f10.3//'   Beta    Qn     Qo        B1         B2         B3
     !    Cn         Co         Cw')
120   format(3f7.3,6e11.3)
      end

************************************************************************

      subroutine Profil(N,NMat,x,MatNum,xSurf,Beta,Ah,AK,ATh,thr,thS,
     !                  ConS,hS,lScreen,ierr)

      dimension x(N),MatNum(N),Beta(N),Ah(N),AK(N),ATh(N),
     !          thr(NMat),thS(NMat),ConS(NMat),hS(NMat)
      logical lScreen

      if(lScreen) write(*,*)'printing profile information'
      write(78,110,err=901)
      ConSN=ConS(MatNum(N))*AK(N)
      do 11 i=N,1,-1
        M=MatNum(i)
        write(78,120,err=901) N-i+1,xSurf-x(i),thr(M),thr(M)+ATh(i)*
     !                        (thS(M)-thr(M)),hS(M)*Ah(i),ConS(M)*AK(i),
     !                        ConS(M)*AK(i)/ConSN,Beta(i),Ah(i),
     !                        AK(i),ATh(i)
11    continue
      write(78,'(''end'')',err=901)
      return

*     Error when writing into an output file
901   ierr=1
      return

110   format(//'    n      depth     THr       THs       hs       Ks',
     !'        Ks/KsTop     Beta      Ah        AK        ATh'/)
120   format(i5,f10.2,2f10.3,f10.1,e12.3,5f10.3)
      end

************************************************************************

*     Read information about solute transport

      subroutine ChemIn(lUpW,lTDep,NMat,NS,NSD,MaxItC,ChPar,TDep,kTopCh,
     !                  cTop,kBotCh,cBot,epsi,tPulse,CumCh,cTolA,cTolR,
     !                  lLinear,lEquil,lArtD,PeCr,lScreen,dSurf,cAtm,
     !                  ierr)

      logical lUpW,lTDep,lLinear(NS),lEquil,lArtD,lScreen
      dimension ChPar(NSD*16+3,NMat),TDep(NS*16+3),cTop(NS),cBot(NS),
     !          CumCh(4,NS)

      if(lScreen) write(*,*) 'reading solute transport information'
      write(50,110,err=902)
      read(30,*,err=901)
      read(30,*,err=901)
      read(30,*,err=901) epsi,lUpW,lArtD,lTDep,cTolA,cTolR,MaxItC,PeCr,
     !                   NS
      PeCr=amax1(PeCr,0.1)
      if(lUpW) then
        write(50,120,err=902)
      else
        write(50,130,err=902)
        if(lArtD) write(50,140,err=902) PeCr
      end if
      write(50,150,err=902) lTDep,cTolA,cTolR,MaxItC
      read(30,*,err=901)
      lEquil=.true.
      do 11 M=1,NMat
        read(30,*,err=901) (ChPar(j,M),j=1,3)
        write(50,160,err=902) M,(ChPar(j,M),j=1,3)
        if(ChPar(3,M).lt.1.) lEquil=.false.
11    continue
      do 13 jj=1,NS
        jjj=(jj-1)*16
        write(50,170,err=902) jj
        read(30,*,err=901)
        read(30,*,err=901) (ChPar(jjj+j,1),j=4,5)
        write(50,180,err=902) (ChPar(jjj+j,1),j=4,5)
        read(30,*,err=901)
        lLinear(jj)=.true.
        do 12 M=1,NMat
          ChPar(jjj+4,M)=ChPar(jjj+4,1)
          ChPar(jjj+5,M)=ChPar(jjj+5,1)
          read(30,*,err=901) (ChPar(jjj+j,M),j=6,19)
          write(50,190,err=902) M,(ChPar(jjj+j,M),j=6,19)
          if(abs(ChPar(jjj+7,M)-0.0).gt.0.001) then
            write(50,200,err=902) M
          else if(abs(ChPar(jjj+8,M)-1.).gt.0.001) then
            write(50,210,err=902) M
          else
            write(50,220,err=902) M
          end if
          if(abs(ChPar(jjj+7,M)-0.0).gt.0.001.or.
     !       abs(ChPar(jjj+8,M)-1.0).gt.0.001) lLinear(jj)=.false.
12      continue
13    continue
      do 14 jj=1,NS*16+3
       TDep(jj)=0.
14    continue
      do 16 jj=1,NS
        do 15 i=1,4
          CumCh(i,jj)=0.
15      continue
        if(lTDep) then
          jjj=(jj-1)*16
          if(jj.eq.1) then
		   read(30,*,err=901)
	    endif
          read(30,*,err=901)
          read(30,*,err=901) (TDep(jjj+j),j=4,5)
          read(30,*,err=901)
          read(30,*,err=901) (TDep(jjj+j),j=6,19)
        end if
16    continue
      read(30,*,err=901)
      read(30,*,err=901) kTopCh,(cTop(jj),jj=1,NS),kBotCh,
     !                   (cBot(jj),jj=1,NS)
      if(kTopCh.eq.-2) then 
        read(30,*,err=901)
        read(30,*,err=901) dSurf,cAtm
      end if
      write(50,230,err=902) kTopCh,(cTop(jj),jj=1,NS)
      write(50,240,err=902) kBotCh,(cBot(jj),jj=1,NS)
      read(30,*,err=901)
      read(30,*,err=901) tPulse
      write(50,250,err=902) tPulse
      return

*     Error when reading from an input file
901   ierr=1
      return
*     Error when writing into an output file
902   ierr=2
      return

110   format(//' Solute transport information'/1X,28('='))
120   format(/' Upstream weighting finite-element method')
130   format(/' Galerkin finite-element method')
140   format (/' Artificial dispersion is added when Peclet number is',
     !         ' higher than',f10.3)
150   format(//' lTDep     cTolA     cTolR   MaxItC'/l3,e13.3,f10.4,i7/
     !        //' Mat.     Bulk.D.    DispL    Fraction')
160   format(i3,f13.4,3f10.4)
170   format(/'    Dif.w.      Dif.g.   ',50('-'),' (',i2,'.solute)')
180   format(2e12.4/' Mat.     KS         Nu         Beta      Henry
     !  SinkL1     SinkS1     SinkG1     SinkL1`    SinkS1`    SinkG1`
     !  SinkL0     SinkS0     SinkG0      Alfa')
190   format(i4,14e11.4)
200   format(/' Langmuir nonlinear adsorption isotherm for material ',
     !       i2)
210   format(/' Freundlich nonlinear adsorption isotherm for material ',
     !       i2)
220   format(/' No adsorption or linear adsorp. isotherm for material ',
     !       i2)
230   format(/' kTopCh      cTop(1...NS)'/i4,7x,6e10.3)
240   format(/' kBotCh      cBot(1...NS)'/i4,7x,6e10.3)
250   format(/' tPulse =   ',f15.3)
      end

* ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||