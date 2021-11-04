* Source file TEMPER.FOR |||||||||||||||||||||||||||||||||||||||||||||||

*     Calculation of heat transport

      subroutine Temper(N,NMat,x,dt,t,MatNum,TempO,TempN,TPar,Ampl,B,D,
     !                  E,F,vOld,vNew,ThOld,ThNew,Cap,Cond,Sink,tPeriod,
     !                  kTopT,tTop,kBotT,tBot)

      double precision B,D,E,F
      dimension x(N),MatNum(N),TempO(N),TempN(N),TPar(10,NMat),B(N),
     !          D(N),E(N),F(N),vOld(N),vNew(N),ThOld(N),ThNew(N),Cap(N),
     !          Cond(N),Sink(N)

      Cw=TPar(9,1)
      PI=3.141592654
      tTopA=tTop
      if(tPeriod.gt.0.) tTopA=tTop+Ampl*sin(2.*PI*t/tPeriod-7.*PI/12.)
      do 11 i=1,N
        TempO(i)=TempN(i)
11    continue

      do 14 Level=1,2
        do 12 i=1,N
          M=MatNum(i)
          if(Level.eq.1) then
            th=ThOld(i)
            v=vOld(i)
          else
            th=ThNew(i)
            v=vNew(i)
          end if
          Cap(i)=TPar(7,M)*TPar(1,M)+TPar(8,M)*TPar(2,M)+TPar(9,M)*th
          Cond(i)=TPar(4,M)+TPar(5,M)*th+TPar(6,M)*sqrt(th)+
     !            TPar(9,M)*TPar(3,M)*abs(v)
12      continue
        dx=x(2)-x(1)
        if(kBotT.gt.0) then
          D(1)=1.
          E(1)=0.
          F(1)=tBot
        else if(kBotT.lt.0) then
          if(Level.eq.2) then
            D(1)=dx/12./dt*(3.*Cap(1)+Cap(2))+(Cond(1)+Cond(2))/dx/4.+
     !           Cw*(2.*vNew(1)+vNew(2))/12.+
     !           dx/24.*Cw*(3.*Sink(1)+Sink(2))
            E(1)=-(Cond(1)+Cond(2))/4./dx+dx/12./dt*(Cap(1)+Cap(2))+
     !           Cw*(2.*vNew(2)+vNew(1))/12.+dx/24.*Cw*(Sink(1)+Sink(2))
          else
            F(1)=TempO(1)*(dx/12./dt*(3.*Cap(1)+Cap(2))-
     !                     (Cond(1)+Cond(2))/dx/4.-
     !                     Cw*(2.*vOld(1)+vOld(2))/12.-
     !                     dx/24.*Cw*(3.*Sink(1)+Sink(2)))+
     !           TempO(2)*((Cond(1)+Cond(2))/4./dx+
     !                     dx/12./dt*(Cap(1)+Cap(2))-
     !                     Cw*(2.*vOld(2)+vOld(1))/12.-
     !                     dx/24.*Cw*(Sink(1)+Sink(2)))+
     !                     tBot*Cw*(vNew(1)+vOld(1))/2.
          end if
        else
          D(1)=-1.
          E(1)=1.
          F(1)=0.
        end if
        do 13 i=2,N-1
          dxA=x(i)-x(i-1)
          dxB=x(i+1)-x(i)
          dx=(x(i+1)-x(i-1))/2.
          if(Level.eq.2) then
            B(i)=-(Cond(i)+Cond(i-1))/4./dxA-Cw*(vNew(i)+
     !           2.*vNew(i-1))/12.+dxA/12./dt*(Cap(i-1)+Cap(i))+
     !           dxA/24.*Cw*(Sink(i-1)+Sink(i))
            D(i)=(Cond(i-1)+Cond(i))/4./dxA+(Cond(i)+Cond(i+1))/4./dxB+
     !           (dxA*(Cap(i-1)+3.*Cap(i))+
     !            dxB*(Cap(i+1)+3.*Cap(i)))/12./dt+
     !           Cw*(vNew(i+1)-vNew(i-1))/12.+
     !           dxA/24.*Cw*(Sink(i-1)+3.*Sink(i))+
     !           dxB/24.*Cw*(3.*Sink(i)+Sink(i+1))
            E(i)=-(Cond(i)+Cond(i+1))/4./dxB+
     !           dxB/12./dt*(Cap(i+1)+Cap(i))+
     !           Cw*(2.*vNew(i+1)+vNew(i))/12.+
     !           dxB/24*Cw*(Sink(i+1)+Sink(i))
          else
            F(i)=TempO(i-1)*(dxA/12./dt*(Cap(i-1)+Cap(i))+
     !                       (Cond(i)+Cond(i-1))/4./dxA+
     !                       Cw*(vOld(i)+2.*vOld(i-1))/12.-
     !                       dxA/24.*Cw*(Sink(i-1)+Sink(i)))+
     !           TempO(i)*(-Cw*(vOld(i+1)-vOld(i-1))/12.+
     !                     (dxA*(Cap(i-1)+3.*Cap(i))+
     !                     dxB*(3.*Cap(i)+Cap(i+1)))/12./dt-
     !                     (Cond(i+1)+Cond(i))/4./dxB-
     !                     (Cond(i)+Cond(i-1))/4./dxA-
     !                     dxA/24.*Cw*(Sink(i-1)+3.*Sink(i))-
     !                     dxB/24.*Cw*(3.*Sink(i)+Sink(i+1)))+
     !           TempO(i+1)*(dxB/12./dt*(Cap(i)+Cap(i+1))+
     !                       (Cond(i+1)+Cond(i))/4./dxB-
     !                       Cw*(2.*vOld(i+1)+vOld(i))/12.-
     !                       dxB/24.*Cw*(Sink(i+1)+Sink(i)))
          end if
13      continue
        if(kTopT.gt.0) then
          B(N)=0.
          D(N)=1.
          F(N)=tTopA
        else if(kTopT.lt.0) then
          dx=x(N)-x(N-1)
          if(Level.eq.2) then
            B(N)=-(Cond(N)+Cond(N-1))/4./dx-
     !           Cw*(vNew(N)+2.*vNew(N-1))/12.+
     !           dx/12./dt*(Cap(N-1)+Cap(N))+
     !           dx/24.*Cw*(Sink(N-1)+Sink(N))
            D(N)=dx/12./dt*(Cap(N-1)+3.*Cap(N))+
     !           (Cond(N-1)+Cond(N))/4./dx-
     !           Cw*(2.*vNew(N)+vNew(N-1))/12.+
     !           dx/24.*Cw*(Sink(N-1)+3.*Sink(N))
          else
            F(N)=TempO(N-1)*((Cond(N)+Cond(N-1))/4./dx+
     !                       Cw*(vOld(N)+2.*vOld(N-1))/12.+
     !                       dx/12./dt*(Cap(N-1)+Cap(N))-
     !                       dx/24.*Cw*(Sink(N-1)+Sink(N)))+
     !           TempO(N)*(dx/12./dt*(Cap(N-1)+3.*Cap(N))-
     !                     (Cond(N-1)+Cond(N))/4./dx+
     !                     Cw*(2.*vOld(N)+vOld(N-1))/12.-
     !                     dx/24.*Cw*(Sink(N-1)+3.*Sink(N)))-
     !           tTopA*Cw*(vNew(N)+vOld(N))/2.
          end if
        end if
14    continue

*     Solve matrix equation
      call BanSol(N,B,D,E,F)
      do 15 i=1,N
        TempN(i)=sngl(F(i))
15    continue
      return
      end

* ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||