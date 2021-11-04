* Source file MATERIAL.FOR |||||||||||||||||||||||||||||||||||||||||||||

      real function FK(iModel,h,Par)

      implicit double precision(A-H,O-Z)
      double precision n,m,Ks,Kr,Kk,Lambda
      real h,Par(9)
      integer PPar

      Qr=Par(1)
      Qs=Par(2)
      Alfa=Par(3)
      n=Par(4)
      Ks=Par(5)
      if(iModel.le.1) then
        BPar=.5d0
        PPar=2
        Qm=Par(6)
        Qa=Par(7)
        Qk=Par(8)
        Kk=Par(9)
        m=1.d0-1.d0/n
        HMin=-1.d300**(1.d0/n)/max(Alfa,1.d0)
        HH=max(dble(h),HMin)
        Qees=dmin1((Qs-Qa)/(Qm-Qa),.999999999999999d0)
        Qeek=dmin1((Qk-Qa)/(Qm-Qa),Qees)
        Hs=-1.d0/Alfa*(Qees**(-1.d0/m)-1.d0)**(1.d0/n)
        Hk=-1.d0/Alfa*(Qeek**(-1.d0/m)-1.d0)**(1.d0/n)
        if(dble(h).lt.Hk) then
          Qee=(1.d0+(-Alfa*HH)**n)**(-m)
          Qe =(Qm-Qa)/(Qs-Qa)*Qee
          Qek=(Qm-Qa)/(Qs-Qa)*Qeek
          FFQ =1.d0-(1.d0-Qee **(1.d0/m))**m
          FFQk=1.d0-(1.d0-Qeek**(1.d0/m))**m
          if(FFQ.le.0.d0) FFQ=m*Qee**(1.d0/m)
          Kr=(Qe/Qek)**Bpar*(FFQ/FFQk)**PPar*Kk/Ks
          FK=sngl(max(Ks*Kr,1.d-37))
          return
        end if
        if(dble(h).ge.Hk.and.dble(h).lt.Hs) then
          Kr=(1.d0-Kk/Ks)/(Hs-Hk)*(dble(h)-Hs)+1.d0
          FK=sngl(Ks*Kr)
        end if
        if(dble(h).ge.Hs) FK=sngl(Ks)
      else
        BPar=1.d0
        Lambda=2.d0  !  !=2 for Mualem Model, =1.5 for Burdine model
        Hs=-1.d0/Alfa
        if(h.lt.Hs) then
          Kr=1.d0/(-Alfa*h)**(n*(BPar+Lambda)+2.d0)
          FK=sngl(max(Ks*Kr,1d-37))
        else
          FK=sngl(Ks)
        end if
      end if
      return
      end

************************************************************************

      real function FC(iModel,h,Par)

      implicit double precision(A-H,O-Z)
      double precision n,m
      real h,Par(9)

      Qr=Par(1)
      Qs=Par(2)
      Alfa=Par(3)
      n=Par(4)
      if(iModel.le.1) then
        Qm=Par(6)
        Qa=Par(7)
        m=1.d0-1.d0/n
        HMin=-1.d300**(1.d0/n)/max(Alfa,1.d0)
        HH=max(dble(h),HMin)
        Qees=dmin1((Qs-Qa)/(Qm-Qa),.999999999999999d0)
        Hs=-1.d0/Alfa*(Qees**(-1.d0/m)-1.d0)**(1.d0/n)
        if(dble(h).lt.Hs) then
          C1=(1.d0+(-Alfa*HH)**n)**(-m-1.d0)
          C2=(Qm-Qa)*m*n*(Alfa**n)*(-HH)**(n-1.d0)*C1
          FC=sngl(max(C2,1.d-37))
          return
        else
          FC=0.0
        end if
      else
        Hs=-1.d0/Alfa
        if(h.lt.Hs) then
          C2=(Qs-Qr)*n*Alfa**(-n)*(-h)**(-n-1.d0)
          FC=sngl(max(C2,1d-37))
        else
          FC=0.0
        end if
      end if
      return
      end

************************************************************************

      real function FQ(iModel,h,Par)

      implicit double precision(A-H,O-Z)
      double precision n,m
      real h,Par(9)

      Qr=Par(1)
      Qs=Par(2)
      Alfa=Par(3)
      n=Par(4)
      if(iModel.le.1) then
        Qm=Par(6)
        Qa=Par(7)
        m=1.d0-1.d0/n
        HMin=-1.d300**(1.d0/n)/max(Alfa,1.d0)
        HH=max(dble(h),HMin)
        Qees=dmin1((Qs-Qa)/(Qm-Qa),.999999999999999d0)
        Hs=-1.d0/Alfa*(Qees**(-1.d0/m)-1.d0)**(1.d0/n)
        if(dble(h).lt.Hs) then
          Qee=(1.d0+(-Alfa*HH)**n)**(-m)
          FQ=sngl(max(Qa+(Qm-Qa)*Qee,1.d-37))
          return
        else
          FQ=sngl(Qs)
        end if
      else
        Hs=-1.d0/Alfa
        if(h.lt.Hs) then
          Qee=(-Alfa*h)**(-n)
          FQ=sngl(max(Qr+(Qs-Qr)*Qee,1d-37))
        else
          FQ=sngl(Qs)
        end if
      end if
      return
      end

************************************************************************

      real function FH(iModel,Qe,Par)

      implicit double precision(A-H,O-Z)
      double precision n,m
      real Qe,Par(9)

      Qr=Par(1)
      Qs=Par(2)
      Alfa=Par(3)
      n=Par(4)
      if(iModel.le.1) then
        Qm=Par(6)
        Qa=Par(7)
        m=1.d0-1.d0/n
        HMin=-1.d300**(1.d0/n)/max(Alfa,1.d0)
        QeeM=(1.d0+(-Alfa*HMin)**n)**(-m)
        Qee=dmin1(dmax1(Qe*(Qs-Qa)/(Qm-Qa),QeeM),.999999999999999d0)
        FH=sngl(max(-1.d0/Alfa*(Qee**(-1.d0/m)-1.d0)**(1.d0/n),-1.d37))
      else
        FH=sngl(max(-1.d0/Alfa*Qe**(-1.d0/n),-1.d37))
      end if
      return
      end

************************************************************************

      real function FS(iModel,h,Par)

      implicit double precision(A-H,O-Z)
      double precision n,m
      real h,Par(9)

      Qr=Par(1)
      Qs=Par(2)
      Alfa=Par(3)
      n=Par(4)
      if(iModel.le.1) then
        Qm=Par(6)
        Qa=Par(7)
        m=1.d0-1.d0/n
        Qees=dmin1((Qs-Qa)/(Qm-Qa),.999999999999999d0)
        Hs=-1.d0/Alfa*(Qees**(-1.d0/m)-1.d0)**(1.d0/n)
        if(h.lt.Hs) then
          HMin=-1.d300**(1./n)/max(Alfa,1.d0)
          HH=max(dble(h),HMin)
          Qee=(1.d0+(-Alfa*HH)**n)**(-m)
          Qe=Qee*(Qm-Qa)/(Qs-Qa)
          FS=sngl(max(Qe,1.d-37))
        else
          FS=1.0
        end if
      else
        Hs=-1.d0/Alfa
        if(h.lt.Hs) then
          Qe=(-Alfa*h)**(-n)
          FS=sngl(max(Qe,1.d-37))
        else
          FS=1.
        end if
      end if
      return
      end

* ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||