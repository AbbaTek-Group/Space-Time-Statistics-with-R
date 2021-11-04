* Source file TIME.FOR |||||||||||||||||||||||||||||||||||||||||||||||||

      subroutine TmCont(dt,dtMaxW,dtOpt,dMul,dMul2,dtMin,Iter,tPrint,
     !                  tAtm,t,tMax,dtMaxC,ItMin,ItMax,lMinStep,dtInit)
      logical lMinStep

      if(lMinStep) then
        dtMax=amin1(dtMaxW,dtMaxC,dtInit)
        dtOpt=dtMax
        lMinStep=.false.
      else
        dtMax=amin1(dtMaxW,dtMaxC)
      end if
      tFix=amin1(tPrint,tAtm,tMax)
      if(Iter.le.ItMin.and.(tFix-t).ge.dMul*dtOpt) 
     !  dtOpt=amin1(dtMax,dMul*dtOpt)
      if(Iter.ge.ItMax)
     !  dtOpt=amax1(dtMin,dMul2*dtOpt)
      dt=amin1(dtOpt,tFix-t)
      dt=amin1((tFix-t)/anint((tFix-t)/dt),dtMax)
      if(tFix-t.ne.dt.and.dt.gt.(tFix-t)/2.) dt=(tFix-t)/2.
      return
      end

************************************************************************

      subroutine SetBC(tMax,tAtm,rTop,rRoot,rBot,hCritA,hBot,hTop,GWL0L,
     !                 TopInF,BotInF,cTop,cBot,NS,tTop,tBot,Ampl,lTemp,
     !                 lChem,ierr,lMinStep)

      logical TopInF,BotInF,lTemp,lChem,lMinStep
      dimension cTop(NS),cBot(NS),cT(5),cB(5)
      character*3 TheEnd

      read(31,101) TheEnd
      if(TheEnd.eq.'end') then
        tMax=tAtm
        return
      else
        backspace 31
        if(.not.lChem.and..not.lTemp) then
          read(31,*,err=901) tAtm,Prec,rSoil,rR,hCA,rB,hB,hT
        else if(lTemp.and..not.lChem) then
          read(31,*,err=901) tAtm,Prec,rSoil,rR,hCA,rB,hB,hT,tTop,tBot,
     !                       Ampl
        else
          read(31,*,err=901) tAtm,Prec,rSoil,rR,hCA,rB,hB,hT,tTop,tBot,
     !                       Ampl,(cT(jj),cB(jj),jj=1,NS)
        end if
      end if

*     Top of the profile
      if(TopInF) then
        rTopOld=rTop
        rTop=abs(rSoil)-abs(Prec)
        if(abs(rTopOld-rTop).gt.1.e-8) lMinStep=.true.
        if(abs(hTop-hT).gt.1.e-8) lMinStep=.true.
        hTop=hT
        rRoot=abs(rR)
        hCritA=-abs(hCA)
        do 11 jj=1,NS
          cTop(jj)=cT(jj)
          cBot(jj)=cB(jj)
11      continue
      end if

*     Bottom of the profile
      if(BotInF) then
        if(abs(rBot-rB).gt.1.e-8) lMinStep=.true.
        rBot=rB
        if(abs(hBot-hB-GWL0L).gt.1.e-8) lMinStep=.true.
        hBot=hB+GWL0L
      end if
      return

*     Error when reading from an input file 
901   ierr=1
      return

101   format(a3)
      end

* ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||