* Source file SOLUTE.FOR |||||||||||||||||||||||||||||||||||||||||||||||

*     To assemble and solve the solute transport equation
*     Mass-lumping finite elements

      subroutine Solute(N,NMat,NS,NSD,x,dt,t,tPulse,ChPar,MatNum,thO,
     !                  thN,vO,vN,Disp,epsi,kTopCh,cTop,kBotCh,cBot,
     !                  Conc,B,D,E,F,g0,g1,Retard,cvTop,cvBot,cvCh0,
     !                  cvCh1,lUpW,wc,Peclet,Courant,dtMaxC,TempO,TempN,
     !                  cNew,cPrevO,cTemp,TDep,thSat,cTolA,cTolR,IterC,
     !                  MaxItC,vCorr,Sorb,SorbN,lLinear,lEquil,lArtD,
     !                  PeCr,q0,q1,dSurf,cAtm,ierr)

      logical lUpW,lConv,lLinear(NSD),lEquil,lArtD
      double precision B,D,E,F
      dimension x(N),ChPar(NSD*16+3,NMat),MatNum(N),thO(N),thN(N),vO(N),
     !          vN(N),Disp(N),cTop(NS),cBot(NS),Conc(NSD,N),B(N),D(N),
     !          E(N),F(N),g0(N),g1(N),Retard(N),cvTop(NS),cvBot(NS),
     !          cvCh0(NS),cvCh1(NS),wc(N),TempO(N),TempN(N),cNew(N),
     !          cPrevO(N),cTemp(N),TDep(NSD*16+3),thSat(NMat),
     !          vCorr(N),Sorb(NSD,N),SorbN(N),q0(N),q1(N)

      alf=1.-epsi
      IterC=1.
      NLevel=2
      Peclet=0.
      Courant=0.
      dtMaxC=1.e+30

*     Loop on species in the chain

      do 19 jj=1,NS
        Iter=0
        jjj=(jj-1)*16
        cvTop(jj)=0.
        cvBot(jj)=0.
        cvCh0(jj)=0.
        cvCh1(jj)=0.
        if(t.gt.tPulse) cTop(jj)=0.
        if(kBotCh.le.0) cvBot(jj)=alf*Conc(jj,1)*vO(1)
        if(kTopCh.lt.0.and.vO(N).lt.0.) cvTop(jj)=alf*cTop(jj)*vO(N)
        if(kTopCh.eq.-2) then
          M=MatNum(N)
          Tr=293.15
          R=8.314
          jjj=(jj-1)*16
          TT=(TempO(N)+273.15-Tr)/R/(TempO(N)+273.15)/Tr
          Dg=ChPar(5,M)*exp(TDep(5)*TT)
          Henry=ChPar(jjj+9,M)*exp(TDep(jjj+9)*TT)
          cvTop(jj)=cvTop(jj)+alf*Dg/dSurf*Henry*Conc(jj,N)-
     !              Dg/dSurf*cAtm
        end if
        if(.not.lLinear(jj)) then
          do 11 i=1,N
            cNew(i)=Conc(jj,i)
            if(.not.lEquil) SorbN(i)=Sorb(jj,i)
11        continue
        end if

*       Iterative loop for a nonlinear adsorption isotherm

12      Iter=Iter+1
        if(.not.lLinear(jj)) then
          do 13 i=1,N
            cTemp(i)=cNew(i)
13        continue
        end if

*     To construct the matrix equation
        do 15 Level=1,NLevel

*     Calculate the dispersion coefficients, retardation factors, source/
*     decay coefficients, Peclet and Courant numbers, upstream weighting
*     factors
          call Coeff(jj,Level,NLevel,N,NMat,NSD,x,Disp,vO,vN,thO,
     !               thN,thSat,ChPar,MatNum,TempN,TempO,TDep,g0,g1,
     !               Retard,Conc,cNew,cPrevO,dt,Pecl,Cour,dtMxC,
     !               lLinear,lEquil,lUpW,lArtD,Iter,wc,vCorr,SorbN,
     !               Sorb,epsi,PeCr,q0,q1)
          Peclet=amax1(Peclet,Pecl)
          Courant=amax1(Courant,Cour)
          dtMaxC=amin1(dtMaxC,dtMxC)

*     Set up the matrix equation
          call MatSet(jj,N,NS,NSD,Level,epsi,alf,dt,kBotCh,kTopCh,cBot,
     !                cTop,x,thO,thN,vO,vN,Conc,Disp,Retard,wc,g0,g1,B,
     !                D,E,F,E1,D1,F1,BN,DN,FN,NMat,ChPar,TempO,TempN,
     !                TDep,dSurf,cAtm,MatNum)

          do 14 i=1,N
            if(Level.eq.1) vO(i)=vO(i)-vCorr(i)
            if(Level.eq.2) vN(i)=vN(i)-vCorr(i)
            if(Level.eq.1.and.Iter.eq.1.and.i.ne.N) then
              dx=x(i+1)-x(i)
              cvCh0(jj)=cvCh0(jj)+alf*dx*(q0(i)+q0(i+1))/2.
              cvCh1(jj)=cvCh1(jj)+alf*dx*(q1(i)+q1(i+1))/2.
            end if
14        continue

15      continue

*     Solve matrix equation
        call BanSol(N,B,D,E,F)

*       Test for convergence for nonlinear problem
        lConv=.true.
        do 16 i=1,N
          if(NS.gt.1.and.Iter.eq.1) cPrevO(i)=Conc(jj,i)
          if(lLinear(jj)) then
            Conc(jj,i)=amax1(sngl(F(i)),0.)
            if(Conc(jj,i).lt.1.e-30.and.Conc(jj,i).gt.0.)
     !                                                 Conc(jj,i)=1.e-30
          else
            cNew(i)=sngl(F(i))
            if(cNew(i).lt.1.0E-30) cNew(i)=0.
            if(abs(cNew(i)-cTemp(i)).gt.cTolA+cTolR*Conc(jj,i))
     !        lConv=.false.
          end if
16      continue
        if(.not.lLinear(jj)) then
          if(.not.lConv.and.iter.lt.MaxItC) then
            goto 12
          else if(.not.lConv) then
            ierr=1
          end if
          do 17 i=1,N
            Conc(jj,i)=cNew(i)
            if(.not.lEquil) Sorb(jj,i)=SorbN(i)
17        continue
        else if(.not.lEquil) then
*         Calculate sorbed concentration for linear noneq. adsorption
          Tr=293.15
          R=8.314
          jjj=(jj-1)*16
          do 22 i=1,N
            M=MatNum(i)
            TT=(TempN(i)+273.15-Tr)/R/(TempN(i)+273.15)/Tr
            Frac=ChPar(3,M)*exp(TDep(3)*TT)
            xKs=ChPar(jjj+6,M)*exp(TDep(jjj+6)*TT)
            xNu=ChPar(jjj+7,M)*exp(TDep(jjj+7)*TT)
            fExp=ChPar(jjj+8,M)*exp(TDep(jjj+8)*TT)
            GamS=ChPar(jjj+11,M)*exp(TDep(jjj+11)*TT)
            GamS1=ChPar(jjj+14,M)*exp(TDep(jjj+14)*TT)
            Omega=ChPar(jjj+19,M)*exp(TDep(jjj+19)*TT)
            Sorb(jj,i)=Sorb(jj,i)+dt*Omega*(1.-Frac)*xKs*Conc(jj,i)/
     !            (2.+dt*(Omega+GamS+GamS1))
22        continue          
        end if
        do 18 i=1,N-1
          dx=x(i+1)-x(i)
          cvCh0(jj)=cvCh0(jj)+epsi*dx*(q0(i)+q0(i+1))/2.
          cvCh1(jj)=cvCh1(jj)+epsi*dx*(q1(i)+q1(i+1))/2.
18      continue

*     Set up mass fluxes
        if(kTopCh.lt.0) then
          if(vN(N).lt.0.) cvTop(jj)=cvTop(jj)+epsi*vN(N)*cTop(jj)
        else 
          cvTop(jj)=FN-BN*Conc(jj,N-1)-DN*Conc(jj,N)
        end if
        if(kTopCh.eq.-2) then
          M=MatNum(N)
          Tr=293.15
          R=8.314
          jjj=(jj-1)*16
          TT=(TempN(N)+273.15-Tr)/R/(TempN(N)+273.15)/Tr
          Dg=ChPar(5,M)*exp(TDep(5)*TT)
          Henry=ChPar(jjj+9,M)*exp(TDep(jjj+9)*TT)
          cvTop(jj)=cvTop(jj)+epsi*Dg/dSurf*Henry*Conc(jj,N)-
     !              Dg/dSurf*cAtm
        end if
        if(kBotCh.le.0) then
          cvBot(jj)=cvBot(jj)+epsi*Conc(jj,1)*vN(1)
        else
          cvBot(jj)=D1*Conc(jj,1)+E1*Conc(jj,2)-F1
        end if
        IterC=max0(IterC,Iter)
19    continue
      return
      end

************************************************************************    

*     Calculate the dispersion coefficients, retardation factors, source/
*     decay coefficients, Peclet and Courant numbers, upstream weighting
*     factors

      subroutine Coeff(jj,Level,NLevel,NumNP,NMat,NSD,x,Disp,vO,vN,thO,
     !                 thN,thSat,ChPar,MatNum,TempN,TempO,TDep,g0,g1,
     !                 Retard,Conc,cNew,cPrevO,dt,Peclet,Courant,dtMaxC,
     !                 lLinear,lEquil,lUpW,lArtD,Iter,wc,vCorr,SorbN,
     !                 Sorb,epsi,PeCr,q0,q1)

      logical lUpW,lLinear(NSD),lEquil,lArtD
      dimension x(NumNP),Disp(NumNP),vO(NumNP),vN(NumNP),thO(NumNP),
     !          thN(NumNP),thSat(NMat),ChPar(NSD*16+3,NMat),g0(NumNP),
     !          g1(NumNP),MatNum(NumNP),Conc(NSD,NumNP),TempO(NumNP),
     !          TempN(NumNP),TDep(NSD*16+3),cNew(NumNP),cPrevO(NumNP),
     !          Retard(NumNP),wc(NumNP),Sorb(NSD,NumNP),SorbN(NumNP),
     !          vCorr(NumNP),q0(NumNP),q1(NumNP)

      TanH(z)=(exp(z)-exp(-z))/(exp(z)+exp(-z))

      jjj=(jj-1)*16
      if(jj.gt.1) jj1=jjj-16
      Peclet=0.
      Courant=0.
      CourMax=1.
      dtMaxC=1.e+30
      Tr=293.15
      R=8.314

      do 12 i=NumNP,1,-1
        M=MatNum(i)

*     Temperature dependence
        if(Level.eq.NLevel) then
          ThW=ThN(i)
          ThG=amax1(0.,thSat(M)-ThW)
          v=vN(i)
          if(i.ne.NumNP) vj=vN(i+1)
          if(i.ne.NumNP) Thj=ThN(i+1)
          TT=(TempN(i)+273.15-Tr)/R/(TempN(i)+273.15)/Tr
          if(jj.gt.1) cPrev=Conc(jj-1,i)
        else
          ThW=ThO(i)
          ThG=amax1(0.,thSat(M)-ThW)
          v=vO(i)
          if(i.ne.NumNP) vj=vO(i+1)
          if(i.ne.NumNP) Thj=ThO(i+1)
          TT=(TempO(i)+273.15-Tr)/R/(TempO(i)+273.15)/Tr
          if(jj.gt.1) cPrev=cPrevO(i)
        end if
        ro=ChPar(1,M)*exp(TDep(1)*TT)
        Frac=ChPar(3,M)*exp(TDep(3)*TT)
        Dw=ChPar(4,M)*exp(TDep(4)*TT)
        Dg=ChPar(5,M)*exp(TDep(5)*TT)
        xKs=ChPar(jjj+6,M)*exp(TDep(jjj+6)*TT)
        xNu=ChPar(jjj+7,M)*exp(TDep(jjj+7)*TT)
        fExp=ChPar(jjj+8,M)*exp(TDep(jjj+8)*TT)
        Henry=ChPar(jjj+9,M)*exp(TDep(jjj+9)*TT)
        GamL=ChPar(jjj+10,M)*exp(TDep(jjj+10)*TT)
        GamS=ChPar(jjj+11,M)*exp(TDep(jjj+11)*TT)
        GamG=ChPar(jjj+12,M)*exp(TDep(jjj+12)*TT)
        GamL1=ChPar(jjj+13,M)*exp(TDep(jjj+13)*TT)
        GamS1=ChPar(jjj+14,M)*exp(TDep(jjj+14)*TT)
        GamG1=ChPar(jjj+15,M)*exp(TDep(jjj+15)*TT)
        xMuL=ChPar(jjj+16,M)*exp(TDep(jjj+16)*TT)
        xMuS=ChPar(jjj+17,M)*exp(TDep(jjj+17)*TT)
        xMuG=ChPar(jjj+18,M)*exp(TDep(jjj+18)*TT)
        Omega=ChPar(jjj+19,M)*exp(TDep(jjj+19)*TT)
        if(jj.gt.1) then
          xKsP=ChPar(jj1+6,M)*exp(TDep(jj1+6)*TT)
          xNuP=ChPar(jj1+7,M)*exp(TDep(jj1+7)*TT)
          fExpP=ChPar(jj1+8,M)*exp(TDep(jj1+8)*TT)
          HenryP=ChPar(jj1+9,M)*exp(TDep(jj1+9)*TT)
          GamL1P=ChPar(jj1+13,M)*exp(TDep(jj1+13)*TT)
          GamS1P=ChPar(jj1+14,M)*exp(TDep(jj1+14)*TT)
          GamG1P=ChPar(jj1+15,M)*exp(TDep(jj1+15)*TT)
        end if
        if(Level.eq.NLevel) then
          TTO=(TempO(i)+273.15-Tr)/R/(TempO(i)+273.15)/Tr
          xKsO=ChPar(jjj+6,M)*exp(TDep(jjj+6)*TTO)
          xNuO=ChPar(jjj+7,M)*exp(TDep(jjj+7)*TTO)
          fExpO=ChPar(jjj+8,M)*exp(TDep(jjj+8)*TTO)
          HenryO=ChPar(jjj+9,M)*exp(TDep(jjj+9)*TTO)
          GamSO=ChPar(jjj+11,M)*exp(TDep(jjj+11)*TTO)
          GamS1O=ChPar(jjj+14,M)*exp(TDep(jjj+14)*TTO)
          xMuSO=ChPar(jjj+17,M)*exp(TDep(jjj+17)*TTO)
          OmegaO=ChPar(jjj+19,M)*exp(TDep(jjj+19)*TTO)
          dKs=(xKs-xKsO)/dt
          dNu=(xNu-xNuO)/dt
          ddExp=(fExp-fExpO)/dt
          dHenry=(Henry-HenryO)/dt
          if(i.ne.1) TTi=(TempN(i-1)+273.15-Tr)/R/(TempN(i-1)+273.15)/Tr
          if(i.ne.NumNP) 
     !      TTj=(TempN(i+1)+273.15-Tr)/R/(TempN(i+1)+273.15)/Tr
        else
          TTN=(TempN(i)+273.15-Tr)/R/(TempN(i)+273.15)/Tr
          xKsN=ChPar(jjj+6,M)*exp(TDep(jjj+6)*TTN)
          xNuN=ChPar(jjj+7,M)*exp(TDep(jjj+7)*TTN)
          fExpN=ChPar(jjj+8,M)*exp(TDep(jjj+8)*TTN)
          HenryN=ChPar(jjj+9,M)*exp(TDep(jjj+9)*TTN)
          dKs=(xKsN-xKs)/dt
          dNu=(xNuN-xNu)/dt
          ddExp=(fExpN-fExp)/dt
          dHenry=(HenryN-Henry)/dt
          if(i.ne.1) TTi=(TempO(i-1)+273.15-Tr)/R/(TempO(i-1)+273.15)/Tr
          if(i.ne.NumNP) 
     !      TTj=(TempO(i+1)+273.15-Tr)/R/(TempO(i+1)+273.15)/Tr
        end if
        if(i.ne.1) Henryi=ChPar(jjj+9,MatNum(i-1))*exp(TDep(jjj+9)*TTi)
        if(i.ne.NumNP)
     !    Henryj=ChPar(jjj+9,MatNum(i+1))*exp(TDep(jjj+9)*TTj)

        dSConc=1.
        dConc=1.
        SConcP=1.
        SConc=1.
        if(.not.lEquil) then
          SEqui=0.
          SEquiO=0.
        end if

*       Effects of nonlinear adsorption
        if(.not.lLinear(jj)) then
          cc=Conc(jj,i)
          cMid=(Conc(jj,i)+cNew(i))/2.
          if(Level.eq.NLevel) cc=cNew(i)
          if(cc.gt.0.) then
            dSConc=fExp*cc**(fExp-1.)/(1+xNu*cc**fExp)**2
            SConc=cc**(fExp-1.)/(1.+xNu*cc**fExp)
            if(.not.lEquil) SEqui=xKs*cc**(fExp-1)/(1.+xNu*cc**fExp)
          end if
          if(cMid.gt.0.)
     !      dConc=fExp*cMid**(fExp-1.)/(1.+xNu*cMid**fExp)**2
          if(Level.eq.NLevel.and..not.lEquil.and.Conc(jj,i).gt.0.)
     !      SEquiO=xKsO*Conc(jj,i)**fExpO/(1.+xNuO*Conc(jj,i)**fExpO)
          if(jj.gt.1.and.cPrev.gt.0.)
     !      SConcP=cPrev**(fExpP-1.)/(1.+xNuP*cPrev**fExpP)
        end if

*     Calculate the retardation factors
        Retard(i)=(ro*Frac*xKs*dConc+ThG*Henry)/ThW+1.

*     Calculate the dispersion coefficients
        TauW=ThW**(7./3.)/thSat(M)**2
        TauG=ThG**(7./3.)/thSat(M)**2
        Disp(i)=ChPar(2,M)*abs(v)/ThW+Dw*TauW+ThG/ThW*Dg*Henry*TauG
        if(.not.lArtD.and..not.lUpW) then
          fi=0.
          if(cMid.gt.0.)
     !      fi=6.*ThW*ro*xKs*cMid**(fExp-1.)*
     !         (fExp/(1.+xNu*cMid**fExp)**2-1./(1.+xNu*cMid**fExp))
          DPom=dt/(6.*ThW*(ThW+ro*Frac*xKs*dSConc+ThG*Henry)+fi)
          if(Level.ne.NLevel) then
            Disp(i)=Disp(i)+v*v*DPom
          else
            Disp(i)=Disp(i)-v*v*DPom
          end if
        end if
        if(lArtD) then
          DD=0.
          if(PeCr.ne.0.and.abs(v).gt.1.e-15) DD=v*v*dt/thW/thW/
     !       Retard(i)/PeCr
          if(DD.gt.Disp(i)) Disp(i)=DD
        end if

*     Calculate the adsorbed concentration on kinetic sites
        if(.not.lEquil) then
          SSorb=Sorb(jj,i)
          if(Level.eq.NLevel) then
            if(.not.lLinear(jj)) then
              SorbN(i)=Sorb(jj,i)+dt*
     !             (epsi*    (Omega*((1.-Frac)*SEqui*cc-SorbN(i))-
     !                       (GamS+GamS1)*SorbN(i)+(1.-Frac)*xMuS)+
     !             (1.-epsi)*(OmegaO*((1.-Frac)*SEquiO-SSorb)-
     !                       (GamSO+GamS1O)*SSorb+(1.-Frac)*xMuSO))
              SSorb=SorbN(i)
            else
              Sorb(jj,i)=((2.-(OmegaO+GamSO+GamS1O)*dt)*Sorb(jj,i)+
     !                    dt*(1.-Frac)*OmegaO*xKsO*Conc(jj,i)+
     !                    dt*(1.-Frac)*(xMuSO+xMuS))/
     !                    (2.+dt*(Omega+GamS+GamS1))
              SSorb=Sorb(jj,i)
            end if
          end if
        end if

*       Calculate zero-order coefficient g0
        g0(i)=xMuL*ThW+Frac*ro*xMuS+ThG*xMuG
        q0(i)=xMuL*ThW+ro*xMuS+ThG*xMuG
c        if(.not.lEquil) g0(i)=g0(i)-Omega*ro*((1.-Frac)*SEqui*cc-SSorb)
        if(.not.lEquil) g0(i)=g0(i)+Omega*ro*SSorb
        if(jj.gt.1) then
          cG=cPrev*(GamL1P*ThW+ro*Frac*xKsP*GamS1P*SConcP+
     !              ThG*HenryP*GamG1P)
          if(.not.lEquil) cG=cG+GamS1P*ro*Sorb(jj-1,i)
          g0(i)=g0(i)+cG
          q0(i)=q0(i)+cG
        end if
        if(cMid.gt.0.) g0(i)=g0(i)-ro*Frac*
     !          (cMid**fExp/(1.+xNu*cMid**fExp)*dKs-
     !           xKs*cMid**(2.*fExp)/(1.+xNu*cMid**fExp)**2*dNu+
     !           xKs*alog(cMid)*cMid**fExp/(1.+xNu*cMid**fExp)**2*ddExp)
     
*       Calculate first-order coefficient g1
        g1(i)=-(GamL+GamL1)*ThW-(GamS+GamS1)*ro*Frac*xKs*SConc-
     !         (GamG+GamG1)*ThG*Henry
        if(.not.lEquil) then
          if(lLinear(jj)) then
            g1(i)=g1(i)-Omega*ro*(1.-Frac)*xKs
            if(Level.eq.NLevel) g1(i)=g1(i)+Omega*ro*
     !               (dt*Omega*(1.-Frac)*xKs/(2.+dt*(Omega+GamS+GamS1)))
          else
            g1(i)=g1(i)-Omega*ro*(1.-Frac)*SEqui
          end if
        end if
        q1(i)=(-(GamL+GamL1)*ThW-(GamS+GamS1)*ro*Frac*xKs*SConc-
     !         (GamG+GamG1)*ThG*Henry)*Conc(jj,i)
      if(.not.lEquil) q1(i)=q1(i)-(GamS+GamS1)*ro*Sorb(jj,i)

*     Velocity corrections
       if(i.eq.1) then
         dx=x(2)-x(1)
         derK=(Henryj-Henry)/dx
       else if(i.eq.NumNP) then
         dx=x(NumNP)-x(NumNP-1)
         derK=(Henry-Henryi)/dx
       else
         dx=(x(i+1)-x(i-1))/2.
         derK=(Henryj-Henryi)/dx
       end if
       vCorr(i)=ThG*Dg*TauG*derK
       if(Level.eq.1) vO(i)=vO(i)-vCorr(i)
       if(Level.eq.NLevel) vN(i)=vN(i)-vCorr(i)

*     Calculate the maximum local Peclet and Courant numbers
        if(i.ne.NumNP) then
          dx=x(i+1)-x(i)
          vv=(abs(v)/ThW+abs(vj)/Thj)/2.
          DD=(Disp(i)+Disp(i+1))/2.
          if(Level.eq.NLevel) then
            Pec=99999.
            dtMax=1.e+30
            vMax=amax1(abs(v)/ThW,abs(vj)/Thj)
            RMin=amin1(Retard(i),Retard(i+1))
            if(DD.gt.0.) Pec=abs(vv)*dx/DD
            Cour=vMax*dt/dx/RMin
            Peclet=amax1(Peclet,Pec)
            Courant=amax1(Courant,Cour)
            Cour1=CourMax
            if(.not.lUpW.and..not.lArtD) then
              if(Pec.ne.99999.) Cour1=amin1(1.,PeCr/amax1(0.5,Pec))
            end if
            if(vMax.gt.1.e-20) dtMax=Cour1*dx*RMin/vMax
            dtMaxC=amin1(dtMaxC,dtMax)

*     To calculate upstream weighting factors
          else if(lUpW.and.Iter.eq.1) then
            Pe2=11.
            if(DD.gt.0.) Pe2=dx*vv/DD/2.
            if(abs(vv).lt.1.e-30) then
              wc(i)=0.
            else if(abs(Pe2).gt.10.) then
              if(vv.gt.0.) wc(i)=1.
              if(vv.lt.0.) wc(i)=-1
            else
              wc(i)=1./TanH(Pe2)-1./Pe2
            end if
          end if
        end if
12    continue
      return
      end

************************************************************************

      subroutine MatSet(jj,N,NS,NSD,Level,epsi,alf,dt,kBotCh,kTopCh,
     !                  cBot,cTop,x,thO,thN,vO,vN,Conc,Disp,Retard,wc,
     !                  g0,g1,B,D,E,F,E1,D1,F1,BN,DN,FN,NMat,ChPar,
     !                  TempO,TempN,TDep,dSurf,cAtm,MatNum)

      double precision B,D,E,F
      dimension cBot(NS),cTop(NS),x(N),thO(N),thN(N),vO(N),vN(N),
     !          Conc(NSD,N),Disp(N),Retard(N),wc(N),g0(N),g1(N),B(N),
     !          D(N),E(N),F(N),ChPar(NSD*16+3,NMat),TempO(N),TempN(N),
     !          TDep(NSD*16+3),MatNum(N)

*     Lower boundary condition
      b1=x(2)-x(1)
      if(Level.eq.1) then
        F1=           Conc(jj,1)*
     !     (b1/2./dt*thO(1)*Retard(1)+
     !      alf*(-(thO(1)*Disp(1)+thO(2)*Disp(2))/b1/2.-
     !           ((2.+3.*wc(1))*vO(1)+vO(2))/6.+
     !           b1/12.*(3.*g1(1)+g1(2))))+
     !                Conc(jj,2)*
     !     alf*((thO(1)*Disp(1)+thO(2)*Disp(2))/b1/2.-
     !          (vO(1)+(2.-3.*wc(1))*vO(2))/6.+b1/12.*(g1(1)+g1(2)))+
     !                alf*b1/6.*(2.*g0(1)+g0(2))

*     3. type  BC
        if(kBotCh.eq.-1) F(1)=F1+alf*cBot(jj)*vO(1)
      else
        E1=epsi*(-(thN(1)*Disp(1)+thN(2)*Disp(2))/b1/2.+
     !           (vN(1)+(2.-3.*wc(1))*vN(2))/6.-b1/12.*(g1(1)+g1(2)))
        D1=b1/2./dt*thN(1)*Retard(1)+
     !     epsi*((thN(1)*Disp(1)+thN(2)*Disp(2))/b1/2.+
     !           ((2.+3.*wc(1))*vN(1)+vN(2))/6.-b1/12.*(3.*g1(1)+g1(2)))
        F2=epsi*b1/6.*(2.*g0(1)+g0(2))
        F1=F1+F2

*     1.type BC
        if(kBotCh.eq.1) then
          D(1)=1.
          E(1)=0.
          F(1)=cBot(jj)
        end if

*     3. type  BC
        if(kBotCh.eq.-1) then
          E(1)=E1
          D(1)=D1
          F(1)=F(1)+F2+epsi*cBot(jj)*vN(1)
        end if

*     Free drainage
        if(kBotCh.eq.0) then
          D(1)=-1.
          E(1)=1.
          F(1)=0.
        end if
      end if

      do 11 i=2,N-1
        a1=b1
        b1=x(i+1)-x(i)
        dx=(x(i+1)-x(i-1))/2.
        if(Level.eq.1) then
          F(i)=       Conc(jj,i-1)*
     !       alf*((thO(i-1)*Disp(i-1)+thO(i)*Disp(i))/a1/2.+
     !            ((2.+3.*wc(i-1))*vO(i-1)+vO(i))/6.+
     !            a1/12.*(g1(i-1)+g1(i)))+
     !                Conc(jj,i)*
     !      (dx/dt*thO(i)*Retard(i)+
     !      alf*(-(thO(i-1)*Disp(i-1)+thO(i)*Disp(i))/a1/2.-
     !           (thO(i+1)*Disp(i+1)+thO(i)*Disp(i))/b1/2.-
     !           (vO(i+1)+3.*(wc(i-1)+wc(i))*vO(i)-vO(i-1))/6.+
     !           (a1*(g1(i-1)+3.*g1(i))+b1*(3.*g1(i)+g1(i+1)))/12.))+
     !                Conc(jj,i+1)*
     !      alf*((thO(i+1)*Disp(i+1)+thO(i)*Disp(i))/b1/2.-
     !           (vO(i)+(2.-3.*wc(i))*vO(i+1))/6.+
     !           b1/12.*(g1(i)+g1(i+1)))+
     !              alf*(a1*(g0(i-1)+2.*g0(i))+b1*(2.*g0(i)+g0(i+1)))/6.
        else
          B(i)=epsi*(-(thN(i-1)*Disp(i-1)+thN(i)*Disp(i))/a1/2.-
     !               ((2.+3.*wc(i-1))*vN(i-1)+vN(i))/6.-
     !               a1/12.*(g1(i-1)+g1(i)))
          D(i)=dx/dt*thN(i)*Retard(i)+
     !         epsi*((thN(i-1)*Disp(i-1)+thN(i)*Disp(i))/a1/2.+
     !               (thN(i+1)*Disp(i+1)+thN(i)*Disp(i))/b1/2.+
     !               (vN(i+1)+3.*(wc(i-1)+wc(i))*vN(i)-vN(i-1))/6.-
     !               (a1*(g1(i-1)+3.*g1(i))+b1*(3.*g1(i)+g1(i+1)))/12.)
          E(i)=epsi*(-(thN(i+1)*Disp(i+1)+thN(i)*Disp(i))/b1/2.+
     !               (vN(i)+(2.-3.*wc(i))*vN(i+1))/6.-
     !               b1/12.*(g1(i)+g1(i+1)))
          F(i)=F(i)+epsi*(a1*(g0(i-1)+2.*g0(i))+
     !                    b1*(2.*g0(i)+g0(i+1)))/6.
        end if
11    continue

*     Upper boundary condition
      if(Level.eq.1) then
        FN=           Conc(jj,N-1)*
     !    alf*((thO(N-1)*Disp(N-1)+thO(N)*Disp(N))/b1/2.+
     !       ((2.+3.*wc(N-1))*vO(N-1)+vO(N))/6.+b1/12.*(g1(N-1)+g1(N)))+
     !                Conc(jj,N)*
     !     (b1/2./dt*thO(N)*Retard(N)+
     !      alf*(-(thO(N-1)*Disp(N-1)+thO(N)*Disp(N))/b1/2.+
     !           (vO(N-1)+(2.-3.*wc(N-1))*vO(N))/6.+
     !           b1/12.*(g1(N-1)+3*g1(N))))+
     !                alf*b1/6.*(g0(N-1)+2.*g0(N))

*     3. type BC
        if(kTopCh.le.0) then
          F(N)=FN
          if(vO(N).lt.0.) F(N)=F(N)-alf*vO(N)*cTop(jj)
          if(kTopCh.eq.-2) then
            M=MatNum(N)
            Tr=293.15
            R=8.314
            jjj=(jj-1)*16
            TT=(TempO(N)+273.15-Tr)/R/(TempO(N)+273.15)/Tr
            Dg=ChPar(5,M)*exp(TDep(5)*TT)
            Henry=ChPar(jjj+9,M)*exp(TDep(jjj+9)*TT)
            F(N)=F(N)-alf*Dg/dSurf*Henry*Conc(jj,N)+Dg/dSurf*cAtm
          end if
        end if
      else
        BN=epsi*(-(thN(N-1)*Disp(N-1)+thN(N)*Disp(N))/b1/2.-
     !        ((2.+3.*wc(N-1))*vN(N-1)+vN(N))/6.-b1/12.*(g1(N-1)+g1(N)))
        DN=b1/2./dt*thN(N)*Retard(N)+
     !     epsi*((thN(N-1)*Disp(N-1)+thN(N)*Disp(N))/b1/2.-
     !           (vN(N-1)+(2.-3.*wc(N-1))*vN(N))/6.-
     !           b1/12.*(g1(N-1)+3.*g1(N)))
        FE=epsi*b1/6.*(g0(N-1)+2.*g0(N))
        FN=FN+FE

*     1. type BC
        if(kTopCh.gt.0) then
          B(N)=0.
          D(N)=1.
          F(N)=cTop(jj)

*     3. type BC
        else
          B(N)=BN
          D(N)=DN
          F(N)=F(N)+FE
          if(vN(N).lt.0.) F(N)=F(N)-epsi*vN(N)*cTop(jj)
          if(kTopCh.eq.-2) then
            M=MatNum(N)
            Tr=293.15
            R=8.314
            jjj=(jj-1)*16
            TT=(TempN(N)+273.15-Tr)/R/(TempN(N)+273.15)/Tr
            Dg=ChPar(5,M)*exp(TDep(5)*TT)
            Henry=ChPar(jjj+9,M)*exp(TDep(jjj+9)*TT)
c            F(N)=F(N)-epsi*Dg/dSurf*Henry*Conc(jj,N)
            D(N)=D(N)+epsi*Dg/dSurf*Henry
          end if
        end if
      end if
      return
      end

*************************************************************************

*     Solve matrix equation

      subroutine BanSol(N,A,B,C,F)
      double precision A,B,C,F
      dimension A(N),B(N),C(N),F(N)

      do 11 i=2,N
        B(i)=B(i)-A(i)*C(i-1)/B(i-1)
        F(i)=F(i)-A(i)*F(i-1)/B(i-1)
11    continue
      F(N)=F(N)/B(N)
      do 12 i=2,N
        j=N-i+1
        F(j)=(F(j)-C(j)*F(j+1))/B(j)
12    continue
      return
      end

* ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||