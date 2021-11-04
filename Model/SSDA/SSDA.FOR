c================================================================
c SSDA (Soil Sensor Data Assimilation)
c Program to assimilate soil water observations at several depths 
c in HYDRUS 1D simulations of water flow in soil profiles
c================================================================
      parameter (nntobs=1000, nndepth=100, nnsensor=30,nntPTF=10, 
     @           nstsim=10000)
      parameter (nnDAtimes=1000, nnnode=1001)
c    nntobs:maximum total number of observation times  
c    nndepth: maximum total number of reference depths    
c    nnsensor: maximum total number of sensors at a given depth
c    nntPTF: maximum total number of pedotransfer functions
c    nnDAtimes: maximum total number of data assimilation times
c    nnnode: maximum number of nodes in in soil water simulations with HYDRUS
      Dimension imodtype(nntPTF)
c    imodtype: 0 - van Genuchten, 2 Brooks and Corey
 	Dimension thr(nndepth,nntPTF),ths(nndepth,nntPTF),
     &  alpha(nndepth,nntPTF),vgn(nndepth,nntPTF), 
     &  aks(nndepth,nntPTF),al(nndepth,nntPTF)
	Dimension thr1(nndepth),ths1(nndepth),alpha1(nndepth),
     C vgn1(nndepth)
c    thr: residual water content, cm^3 cm^-3
c    ths: volumetric water content at saturation, cm^3 cm^-3
c    alpha: van Genuchten's alpha or inverse of air entry pressure in
c          Brook Corey equation, cm^ 1
c    vgn: van Genuchten's 'n' or Brooks Corey lambda
c    aks: saturated hydraulic conductivity cm d^ 1
c    al: tortuosity parameter, [] 
c    thr1, ths1, alpha1, vgn1: auxuliary
	Dimension ipick(nndepth)
c    ipick: shows if the reference depth has measurements: 1 - yes, 0 no
	Dimension DATIMES(nnDAtimes)
c    DATIMES times when the data assimilation is performed
 	Dimension x(nnnode),h(nnnode),xref(nndepth),
     &  winit1(nndepth), winit(nndepth,nntPTF),
     & ThNew(nnnode),xhydrus(nnnode),Thetas(nnnode,nntPTF),
     & Thetasav(nnnode)
c    x: elevation of nodes  relative soil surface, cm: x(1)=0,x(i)<0,i>1
c    h: nodal pressure head, cm
c    xref: reference depths, cm; they include measurement depths and 
c          additional interpolation depths
c    winit1: soil water content distribution in soil profile, cm^3 cm^-3;
c            auxiliary variable
c    ThNew: water content simulated with HYDRUS for the time of 
c           the assimilation update
c    xhydrus: nodal coordinates used in HYDRUS1D; they are measured 
c             upwards from the bottom, and nodal count is also up from bottom
c    Thetas: nodal soil water contents simulated with the PTF ensemble
c    for the assimilation update time
c    Thetasav: average simulated nodal soil water content for the
c    assimilation update time
	Dimension tobs(nntobs)
c    tobs: observation times, d
      Dimension wnew(nntobs,nndepth,nnsensor)
CYAP++
      Dimension wnewf(nntobs,nndepth,nnsensor)

	Dimension v(nnsensor,nndepth),aver(nndepth),
     &  cm(nndepth,nndepth),nob(nndepth),corm(nndepth,nndepth)
c    wnew: array to store observed water contents, cm^3 cm^-3
c    v: bias corrected deviation of the observed soil water content 
c       for a specific sensor from the average at a given depth, cm^3 cm^-3
c    aver: average observed water content at a given depth, cm^3 cm^-3 
c    cm: covariance matrix of water content observations, cm^6 cm^-6
c    cm: correlation matrix of water content observations, []
c    nob: number of available observations at a given depth (<= number of sensors) 
	Dimension bb(nnsensor,nndepth),bias(nnsensor,nndepth)
c    bias: sensor biases 
      Dimension winitmn(nndepth,nntPTF), xpoutmn(nndepth,nntPTF)
c    winitmn: soil water content in reference points right before update, cm^3 cm^-3 
c    xpoutmn: soil water content in reference points after update, cm^3 cm^-3
	Dimension tout_obs(nntobs),wout_obs(nntobs,nndepth),
     &  tout_sim(nntobs),wout_sim(nntobs,nndepth)
	Dimension PTFid(nntPTF)
c    PTFidnumber: PTF as numbered in Guber Pachepsky, 2011
	dimension mat(nnnode)
c    PTF error array and average of the assimilated wc if errors
	Dimension errors(nntPTF), winitav(nndepth)
c    nodal soil material id number
	Dimension wout_sim1(nstsim,nndepth,nntPTF),m(nndepth),
     .	wout_sim2(nstsim,nndepth,nntPTF)
c    wout_sim1: simulated soil water contents by depths, PTFs and simulation times, cm^3 cm^-3
c    m: material number for each reference depth
	character cFileName*80,cDataPath*80,LINE*200,A1*1,A2*1,A3*1
      character FileName*80, ProfileFile*80,SelectorFile*80
      character AtmosphFile*80
	integer*2 ptferror, nerrors !GM
c==========================================================
c This common field is used in the HYDRUS subroutine OBSNOD 
c to write down simulation results for each time step
c==========================================================
      common /ptf/iptf
C=====================
      open(10,file='datadir.txt')
	read(10,'(A)') cDataPath
	close(10)
      iLengthPath = Len_Trim(cDataPath)
      cFileName = cDataPath(1:iLengthPath)//'\Input.dat'
	open(20,file=cFileName, status='old')
c----read the number of observation times, reference depths, sensors and materials
        read(20,*)
        read(20,*)
        read(20,*) ntobs, ndepth, nsensor,nmat
c----read the observation depth
        read(20,*)
        read(20,*) (xref(i),i=1,ndepth)
c----read the DA selection for the depths        
	read(20,*)
        read(20,*) (ipick(i),i=1,ndepth)
c----number of measurement depths used in assimilation
	NMD=0
	Do k=1,ndepth
	  if(ipick(k).GT.0) NMD=NMD+1
	Enddo
cc================================
c Read observations
c================================
      Read(20,*)
	mm=1
      Do k=1,ndepth
	  If(ipick(k).EQ.1) then
           read(20,'(A)') filename
           cFileName = cDataPath(1:iLengthPath)//"\"//filename
	     open(23,file=cFileName, status='old')
	     Do i=1,ntobs
	       Read(23,*) tobs(i),(wnew(i,mm,j),j=1,nsensor)
           Enddo
	     Close(23)
	     mm=mm+1
	  Endif
      Enddo
c=======================
c read the sensor biases
c=======================
      Read(20,*)
      Read(20,'(A)') filename
 	cFileName = cDataPath(1:iLengthPath)//"\"//filename
	open(23,file=cFileName,status='old')
	mm=1
	Do k=1,ndepth
	  If(ipick(k).EQ.1) then
          Read(23,*) (bias(j,mm),j=1,nsensor)
	    mm=mm+1
        Endif
	Enddo
      Close(23)
c================================================================
c-----read the HYDRUS template files 
       read(20,*)
       read(20,'(A)')ProfileFile
       read(20,'(A)')AtmosphFile
       read(20,'(A)')SelectorFile
c==================================
c Read data assimilation times
c================================== 
	Read(20,*)
	Read(20,*) NDATIMES
	Read(20,*) LINE
	Read(20,*) (DATIMES(idat),idat=1,Ndatimes+1)
c=========================================================       
c File of initial water contents at all reference depths
c=========================================================
      read(20,*)
      read(20,'(A)') filename
      cFileName = cDataPath(1:iLengthPath)//"\"//filename
	open(23,file=cFileName, status='old')
	Do k=1,ndepth
	   Read(23,*) winit1(k)
      Enddo
	Close(23)
c======================================
c Read hydraulic properties of each PTF
c======================================
      read(20,*)
  	read(20,*) NPTF
  	read(20,*) 
	Do iptf=1,NPTF
        read(20,*) PTFid(iptf),imodtype(iptf)
c-----read hydraulic properties of each material and PTF	
	  Do k=1,nmat 
	    read(20,*) thr(k,iptf),ths(k,iptf),alpha(k,iptf),vgn(k,iptf),
     &       aks(k,iptf),al(k,iptf)
	  Enddo
	Enddo
c-----Done with the 'input.dat' file
      close(20)
c=======================================================================
c Read nodal depths to build pressure head profiles by converting w to h
c and soil materials by nodes
c Reads the default profile.dat file
c=======================================================================
      cFileName = cDataPath(1:iLengthPath)//"\"//ProfileFile
      open(32,file=cFileName, status='old')
	Read(32,*)
	Read(32,*)
	Read(32,*)
	Read(32,*) NNODES
	Do i=1,NNODES
	Read(32,*) idum, x(i),h(i),mat(i)
	Enddo
	close(32)
c========================================================
c Open files for simulation results in reference nodes
c========================================================
      Do iptf=1,NPTF
	  L0=600+iptf
	  L1=iptf/10
	  L2=iptf-10*(L1)
	  A1=char(L1+48)
	  A2=char(L2+48)
        cFileName = cDataPath(1:iLengthPath)//'\sim_'//A1//A2//'.txt'
	  open(L0,file=cFileName)
	Enddo
c================================================================
c Initial soil water content distribution for each ensemble model
c================================================================
	Do k=1,ndepth
	Do iptf=1,NPTF
	 winit(k,iPTF)=winit1(k)
	Enddo
	Enddo
c====================================
c Nodal number for observation depths
c====================================
      Do j=1,ndepth
        difmax=1000.
        Do i=1,NNODES 
        delt=abs(x(i)-xref(j))
        If(delt.LT.difmax) then
          difmax=delt
          m(j)=i
        Endif
        Enddo
      Enddo
c========================================
c check PTF ths vs initial water contents
c========================================
 	Do k=1,ndepth
 	Do iptf=1,NPTF
	if(winit(k,iPTF).LT.thr(mat(m(k)),iptf).OR.
     &	winit(k,iPTF).GT.ths(mat(m(k)),iptf)) then 
	   print *, k,iptf, winit(k,iPTF),thr(k,iptf),ths(k,iptf)
	   stop
	 endif
	Enddo
	Enddo
c=======================
c Data assimilation loop
c=======================
      DO idatimes=1,NDATIMES
      nerrors=0
      t1=DATIMES(idatimes)
 	t2=DATIMES(idatimes+1)
c===========================================================
c  Build atmospheric file for the period from t1 through t2
c  Atmospheric data are on daily basis
c===========================================================
      cFileName = cDataPath(1:iLengthPath)//"\"//AtmosphFile
      open(25,file=cFileName, status='old')
      cFileName = cDataPath(1:iLengthPath)//'\Atmosph.in'
      open(26,file=cFileName, status='unknown')
c
      Do k=1,6
	Read(25,'(a)') LINE
	write(26,'(a)') LINE
	Enddo
c
900	read(25,*) T
	If(T.GT.t1+1-0.1) then
	Backspace(25)
  	Read(25,'(a)') LINE
 	write(26,'(a)') LINE
	endif
	if (T.LT.t2-0.1) goto 900
      write(26,'(a)') 'end'
	close(25)
	close(26)
c=============================================================
c PTF loop
c Build input files and do simulations for each ensemble model
c=============================================================
	DO iptf=1,NPTF
c==============================
c Build the 'selector.in' file
c==============================

      cFileName = cDataPath(1:iLengthPath)//"\"//SelectorFile
	open(21,file=cFileName, status='old')
      cFileName = cDataPath(1:iLengthPath)//'\Selector.in'
	open(22,file=cFileName, status='unknown')
c
      do k =1,21 
	Read(21,'(a)') 	LINE
	write(22,'(a)') LINE
	Enddo
c
      read(21,*) imodel, ihyster 
	write(22,'(I7,4x,I7)') imodtype(iptf), ihyster
	read(21,'(a)') LINE 
	write(22,'(a)') LINE
c
	Do k=1,nmat 
	  Read(21,'(a)') 	LINE
	write(22,'(6F10.4)') thr(k,iptf),ths(k,iptf),
     &    alpha(k,iptf),vgn(k,iptf),aks(k,iptf)
	Enddo
c
      Do k=1,4
	  Read(21,'(a)')  LINE
	  write(22,'(a)') LINE
	Enddo
c
	read(21,*) LINE
	write(22,*)	t1,t2
	Read(21,'(a)') 	LINE
	write(22,'(a)') LINE
 	read(21,*) LINE
      write(22,*) t2     
	Read(21,'(a)') 	LINE
	write(22,'(a)') LINE
	close(21)
	close(22)
c================================================
c  Build the Profile.dat file for the current PTF
c================================================
c-----The auxuliary arrays are filled to convert W to h
      Do k=1,nmat
	thr1(k)=thr(k,iptf)
	ths1(k)=ths(k,iptf)
	alpha1(k)=alpha(k,iptf)
	vgn1(k)=vgn(k,iptf)
	imodel=imodtype(iptf)
	Enddo
      Do k=1,ndepth
		winit1(k)=winit(k,iptf)
	if (winit1(k).gt.ths1(mat(m(k)))) winit1(k)=ths1(mat(m(k)))-0.003 
	Enddo
c
      if(IDATIMES.EQ.1) then 
      Call W_TO_H(NNODES,x,h,mat,ndepth,xref,winit1,thr1,ths1,alpha1,
     &     vgn1,imodel)
	Else
	Do i=1,NNODES
	  Thetasav(i)=Thetas(i,iptf)
	Enddo
      Call W_TO_H_B(NNODES,x,h,mat,ndepth,xref,winit1,thr1,ths1,
     &alpha1,vgn1,Thetasav,imodel)
      Endif

      cFileName = cDataPath(1:iLengthPath)//"\"//ProfileFile
      open(24,file=cFileName, status='old')
      cFileName = cDataPath(1:iLengthPath)//'\Profile.dat'
      open(32,file=cFileName, status='unknown') 
c
	Do k=1,4
	  read(24,'(a)') LINE
	  write(32,'(a)') LINE
	Enddo
c
	Do i=1,NNODES
	  read(24,*) np,xp,hp,Mp,Lp,Bp,Axp,Bxp,Dxp
	  write(32,'(I4,F10.3,F26.3, 2I4,4F10.3)') 
     &      np,xp,h(i),Mp,Lp,Bp,Axp,Bxp,Dxp
	Enddo
CYAP
      write(32,*) ndepth
	write(32,'(100I5)') (m(j),j=1,ndepth)
	close(24)
      close(32)	
	ptferror=0
c===========
c run HYDRUS
c===========
	print *,'Ensemble model ',iptf, ';   Start day = ',t1
      CALL HYDRUS(xhydrus,ThNew,ptferror)
c==========================================
cGM	If errors have been found during the Hydrus run
	if (ptferror.GT.0) then
		errors(iptf)=1
		nerrors=nerrors+1
	endif
c==========================================
c==========================================================
c save nodal water contents and water contents at reference
c depths to use them as winits on DA updates
c==========================================================
	Do k=1,ndepth
	  Do i=1,NNODES
	   Thetas(i,iptf)=ThNew(NNODES+1-i)
	   if(abs(xref(k)-xhydrus(NNODES+1-i)).LT.0.1) 
     &               winit(k,iptf)=ThNew(NNODES+1-i)
 	  Enddo
 	Enddo

	ENDDO ! PTF; done with the ensemple modeling before DA update
c==========================================
c===============
c DA starts here
c===============
c==========================================
c Observation data on the assimilation time
c==========================================
	i2=1
c
	Do i=1,ntobs
	 if(abs(tobs(i)-t2).LT.abs(tobs(i2)-t2)) i2=i
	Enddo
c	
c=================================================================
c Build the average and the covariance matrix (cm) of observations
c=================================================================
	Do mm=1,NMD
	Do j=1,nsensor
        v(j,mm)=wnew(i2,mm,j) 
 	Enddo
	Enddo
c
	Do mm=1,NMD
	  nob(mm)=0
	  aver(mm)=0.0
	  Do j=1,nsensor
	    if(v(j,mm).GT.0.0) then
		   aver(mm)=aver(mm)+v(j,mm)
	       nob(mm)=nob(mm)+1
	    endif
	  Enddo
	  aver(mm)=aver(mm)/float(nob(mm))
	Enddo
c
	Do mm=1,NMD
	  Do j=1,nsensor
          v(j,mm)=v(j,mm)-aver(mm)- bias(j,mm)
	  Enddo
	Enddo
c
	Do k1=1,NMD
	  Do k2=1,NMD
	    cm(k1,k2)=0.
	    s0=0.
          Do j=1,nsensor	  
            if(wnew(i2,k1,j).GT.0.AND.wnew(i2,k2,j).GT.0) then
		    cm(k1,k2)=cm(k1,k2)+ v(j,k1)*v(j,k2)
	        s0=s0+1.
	      endif
	    Enddo
	    cm(k1,k2)=cm(k1,k2)/(s0-1.)
	  Enddo
	Enddo


	Do k=1,ndepth	  
	  i=0
	  Do j=1,NPTF
		if (errors(j).EQ.0) then
		  i=i+1
	      winitmn(k,i)=winit(k,j)		  
		endif
	  Enddo
	Enddo
	NPTF2=NPTF-nerrors
c
      CALL DA(NMD,ndepth,NPTF2,winitmn,aver,cm,ipick,XPOUTmn,
     &        nndepth,nntPTF)
c=======================================================
c Saving DA results to be used as new initial conditions
c=======================================================
c
      Do k=1,ndepth
	  i=0
	  winitav(k)=0
	  Do j=1,NPTF
		if (errors(j).EQ.0) then
		  i=i+1
	      winit(k,j)=XPOUTmn(k,i)
		  winitav(k)=winit(k,j)+winitav(k)
		endif
	  Enddo
		winitav(k)=winitav(k)/NPTF2
	Enddo
c
cGM update of winit in PTF where errors were found with the average of the others
      Do k=1,ndepth
	  i=0
	  Do j=1,NPTF
		if (errors(j).EQ.1) then		  
	      winit(k,j)=winitav(k)
		endif
	  Enddo
	Enddo
cGM================================================================================== End
c
	ENDDO ! DA TIMES

	Do iptf=1,NPTF
	  close(600+iptf)
	Enddo
c========
c Outputs 
c========
c====================================================================================
c Printouts are for the beginning of day. It is relatively easy 
c to change this to the arbitrary time of day or to average over some predefined time
c====================================================================================

c=========================================
c Generate average observed water contents
c=========================================
      iout=1
      i=1
15    if(i.GT.ntobs) goto 16
	if(tobs(i).GE.DATIMES(1).AND.tobs(i).LE.DATIMES(NDATIMES+1))then 
	   tout_obs(iout)=tobs(i)
	   Do j=1,NMD
	     s0=0.
	     s1=0.
	     Do k=1,nsensor
 	      if(wnew(i,j,k).GT.0.0) then      
			 s0=s0+1.0
               s1=s1+wnew(i,j,k)
	       Endif
	     Enddo
	     wout_obs(iout,j)=s1/s0
	   Enddo
	   iout=iout+1
	Endif
	i=i+1
	goto 15
16    continue
	jj=0
      Do j=1,ndepth
	if(ipick(j).EQ.1) then
	jj=jj+1 
	L0=nint(abs(xref(j)))
	L1=L0/100
	L2=(L0-100*L1)/10
	L3=L0-L1*100-L2*10
	  A1=char(L1+48)
	  A2=char(L2+48)
	  A3=char(L3+48)
	  cFileName = cDataPath(1:iLengthPath)//'\wout_obs'
     %       //A1//A2//A3//'.txt'
	  open(L0,file=cFileName)
 	  Do i=1,iout-1
	    write(L0,'(F7.1,F10.4)') tout_obs(i),wout_obs(i,jj) 
	  Enddo
	  close(L0)
	endif
	Enddo

c===================
c Simulation results
c===================

	Do iout=1,NINT(DATIMES(NDATIMES+1)-DATIMES(1)+1)
		tout_sim(iout)=DATIMES(1)+float(iout)-1.0
          Do iptf=1,NPTF
          Do k=1,ndepth
             wout_sim2(iout,k,iptf) = 0.
          Enddo
          Enddo
	Enddo
c =================
c  Initial time
c =================
	  Do k=1,ndepth
	    wout_sim(1,k)=wout_obs(1,k)
        Enddo
c ==
c Loop through times
c
	Do iout=2,NINT(DATIMES(NDATIMES+1)-DATIMES(1)+1)        
c ====
c Loop through PTFs
c ====
      Do iptf=1,nptf       
c =================
c  Input file name 
c =================
	  L0=600+iptf
	  L1=iptf/10
	  L2=iptf-10*(L1) 
	  A1=char(L1+48)
	  A2=char(L2+48)
        cFileName = cDataPath(1:iLengthPath)//'\sim_'//A1//A2//'.txt'
	  open(L0,file=cFileName)
c ==================================================================
c Loop through simulation output to find the outpuit on tout (iout)
c ==================================================================
181	  read(L0,*,ERR=183) t,(wout_sim2(iout,k,iptf),k=1,ndepth)
	  if(abs(t-tout_sim(iout)).LT.1.E-04) goto 183
        goto 181
183     close(L0)
        Enddo ! iPTF
c =====================================
c Average PTF resuts for the given iout
c =====================================
	   Do k=1,ndepth
            wout_sim(iout,k)=0.
            n=0
		  Do iptf=1,nptf
			if (wout_sim2(iout,k,iptf).GT.0) then
			  wout_sim(iout,k)=wout_sim(iout,k)+wout_sim2(iout,k,iptf)
			  n=n+1
			endif
		  Enddo
		  wout_sim(iout,k)=wout_sim(iout,k)/float(n)
         Enddo

      Enddo ! iout
c ========================
c write simulation results
c ========================
c
      Do k=1,ndepth
	  L0=nint(abs(xref(k)))
	  L1=L0/100
	  L2=(L0-100*L1)/10
	  L3=L0-L1*100-L2*10
	  A1=char(L1+48)
	  A2=char(L2+48)
	  A3=char(L3+48)
	  cFileName = cDataPath(1:iLengthPath)//'\wout_sim'//A1//A2//A3//
     &'.txt'
	  open(L0,file=cFileName)
 	  Do iout=1,NINT(DATIMES(NDATIMES+1)-DATIMES(1)+1)
	    write(L0,*) tout_sim(iout),wout_sim(iout,k)
	  Enddo
	  close(L0)
	Enddo
c========================================================================= End GM
	STOP
	END
c=================================
c HYDRUS converted to a subroutine
c=================================
*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
*                                                                      *
*     HYDRUS   - Numerical model of one-dimensional variably saturated *
*                water flow, heat transport, and transport of solutes  *
*                involved in sequential first-order decay reactions    *
*                version 6.0                                           *
*                                                                      *
*                Version coupled with user-friendly interface for      *
*                MS Windows 3.1 environment                            *
*                                                                      *
*     Designed by J.Simunek, K. Huang and M. Th. van Genuchten (1996)  *
*                                                                      *
*     Based on model SWMS_1D - J. Simunek, 1992,                       *
*                                                                      *
*                                       Last modified: February, 1998  *
*                                                                      *
*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*

      SUBROUTINE HYDRUS(x,ThNew,ptferr)

      parameter (NumNPD=1001,
     !           NMatD =20,
     !           NTabD =100,
     !           NObsD =10,
     !           NSD   =5,
     !           NUnitD=8)

      integer PLevel,Alevel,TLevel,err
      logical SinkF,WLayer,qGWLF,TopInF,ShortO,lWat,lChem,lTemp,ConvgF,
     !        FreeD,SeepF,BotInF,AtmBC,lRoot,lUpW,lWDep,lTDep,lEquil,
     !        lLinear(NSD),lArtD,lMoSink,lSolRed,lSolAdd,lScreen,
     !        lMinStep
      double precision P,R,S,Q,RTime1,RTime2
      character cFileName*80,cDataPath*80,cName*12,ch1*1
      integer*2 stat,ptferr

      dimension x(NumNPD),hNew(NumNPD),hOld(NumNPD),hTemp(NumNPD),
     !  MatNum(NumNPD),ParD(9,NMatD),TPrint(50),Sink(NumNPD),thr(NMatD),
     !  POptm(NMatD),Beta(NumNPD),LayNum(NumNPD),CumQ(5),hTab(NTabD),
     !  ConTab(NTabD,NMatD),CapTab(NTabD,NMatD),TheTab(NTabD,NMatD),
     !  Con(NumNPD),Cap(NumNPD),ConSat(NMatD),ths(NMatD),P(NumNPD),
     !  R(NumNPD),S(NumNPD),hSat(NMatD),ParW(9,NMatD),
     !  ChPar(NSD*16+3,NMatD),Conc(NSD,NumNPD),vNew(NumNPD),
     !  vOld(NumNPD),Disp(NumNPD),Retard(NumNPD),g0(NumNPD),g1(NumNPD),
     !  Q(NumNPD),wc(NumNPD),ThNew(NumNPD),ThOld(NumNPD),WatIn(NumNPD),
     !  SolIn(NumNPD),Node(NObsD),iUnit(NUnitD),cvCh0(NSD),cvCh1(NSD),
     !  cTop(NSD),cBot(NSD),cRoot(NSD),CumCh(4,NSD),cVOlI(NSD),
     !  cCumA(NSD),cCumT(NSD),cvTop(NSD),cvBot(NSD),cNew(NumNPD),
     !  cTemp(NumNPD),cPrevO(NumNPD),TDep(NSD*16+3),TPar(10,NMatD),
     !  TempO(NumNPD),TempN(NumNPD),SorbN(NumNPD),Sorb(NSD,NumNPD),
     !  q0(NumNPD),q1(NumNPD),aOsm(NSD),
     !  Ah(NumNPD),AK(NumNPD),ATh(NumNPD),AhW(NMatD),AThW(NMatD),
     !  AKW(NMatD),Kappa(NumNPD),AThS(NumNPD),ThRR(NumNPD),
     !  ConR(NumNPD),ConO(NumNPD),AKS(NumNPD),KappaO(NumNPD)

      data CosAlf,NTab,ItCum,TLevel,ALevel,PLevel,hRoot,vRoot,cRoot
     !    /  1.  ,100 , 0   ,  1   ,  1   ,   1  ,  0. ,  0. ,NSD*0. /
     !
     !     CumQ,Sink     ,wc       ,iUnit                  ,IterW,IterC
     !    /5*0.,NumNPD*0.,NumNPD*0.,50,70,71,72,75,76,77,78, 0   , 0  /
     !
     !     dtMaxC,wCumT,wCumA,cCumT ,cCumA ,ierr,lMinStep
     !    /1.e+30,  0. ,  0. ,NSD*0.,NSD*0.,  0 ,.false./
	Plevel=1
	Tlevel=1
	Alevel=1
	ItCum=0
	IterW=0
	IterC=0
      iCount = NARGS()

      open(10,file='datadir.txt')
	read(10,'(A)') cDataPath
	close(10)

      iLengthPath = Len_Trim(cDataPath)
      cFileName = cDataPath(1:iLengthPath)//'\Selector.in'
      open(30,file=cFileName, status='old',err=901)
      cFileName = cDataPath(1:iLengthPath)//'\Profile.dat'
      open(32,file=cFileName, status='old',err=901)
      cFileName = cDataPath(1:iLengthPath)//'\I_CHECK.OUT'
      open(50,file=cFileName, status='unknown',err=902)
      cFileName = cDataPath(1:iLengthPath)//'\RUN_INF.OUT'
      open(70,file=cFileName, status='unknown',err=902)
      cFileName = cDataPath(1:iLengthPath)//'\T_LEVEL.OUT'
      open(71,file=cFileName, status='unknown',err=902)
      cFileName = cDataPath(1:iLengthPath)//'\NOD_INF.OUT'
      open(75,file=cFileName, status='unknown',err=902)
      cFileName = cDataPath(1:iLengthPath)//'\BALANCE.OUT'
      open(76,file=cFileName, status='unknown',err=902)
      cFileName = cDataPath(1:iLengthPath)//'\OBS_NODE.OUT'
      open(77,file=cFileName, status='unknown',err=902)
      cFileName = cDataPath(1:iLengthPath)//'\A_LEVEL.OUT'
      open(72,file=cFileName, status='unknown',err=902)
      cFileName = cDataPath(1:iLengthPath)//'\PROFILE.OUT'
      open(78,file=cFileName, status='unknown',err=902)

*     Read input data

      call BasInf(CosAlf,MaxIt,TolTh,TolH,TopInF,BotInF,ShortO,lWat,
     !            lChem,SinkF,WLayer,qGWLF,FreeD,SeepF,AtmBC,KodTop,
     !            KodBot,rTop,rRoot,rBot,hCritS,hCritA,GWL0L,Aqh,Bqh,
     !            kTOld,kBOld,NUnitD,iUnit,NMat,NMatD,NLay,lRoot,lTemp,
     !            lWDep,lEquil,lScreen,err)
      if(err.ne.0) goto (905,906,903,916) err
      if(TopInF.or.BotInF) then
        cFileName = cDataPath(1:iLengthPath)//'\ATMOSPH.IN'
        open(31,file=cFileName, status='old',err=901)
      end if
      call NodInf(NumNPD,NumNP,NObsD,NObs,hTop,hBot,x,hNew,hOld,
     !            MatNum,hTemp,LayNum,Beta,Ah,AK,ATh,Conc,Sorb,
     !            TempN,TempO,Node,NSD,NS,xSurf,lChem,lTemp,lEquil,
     !            lScreen,err)
      if(err.ne.0) goto (912,903,914,915,929) err
      call MatIn (NMat,ParD,ParW,hTab(1),hTab(NTab),lScreen,err,NumNP,
     !            Ah,iHyst,AhW,AThW,AKW,MatNum,hNew,Kappa,AThS,ThRR,
     !            ConR,AKS,KappaO,iModel)
      if(err.ne.0) goto (906,903) err
      call GenMat(NTab,NTabD,NMat,thr,ths,hSat,ParD,hTab,ConTab,CapTab,
     !            ConSat,TheTab,iModel,lScreen,err)
      if(err.eq.1) goto 903
      call SetMat(NumNP,NTab,NTabD,NMat,hTab,ConTab,CapTab,hNew,MatNum,
     !            ParD,Con,Cap,ConSat,Ah,AK,ATh,hSat,hTemp,TheTab,ThOld,
     !            thr,ths,lWDep,TempN,IterW,ConO,Kappa,AThS,ThRR,ConR,
     !            AKS,AhW,AThW,AKW,iModel)
      call TmIn  (tInit,tMax,tAtm,tOld,dt,dtMax,DMul,DMul2,dtMin,TPrint,
     !            t,dtOpt,TopInF,BotInF,lScreen,ItMin,ItMax,MaxAL,
     !            hCritS,err)
      if(err.ne.0) goto (907,928) err
      dtInit=dt
      if(lRoot) then
        call RootIn(tRMin,tRHarv,xRMin,xRMax,RGR,lScreen,err)
        if(err.ne.0) goto (908,903) err
        call SetRG(NumNP,x,Beta,t,tRMin,tRHarv,xRMin,xRMax,RGR)
      end if
      if(lTemp) then
        call TempIn(NMat,TPar,Ampl,tPeriod,kTopT,tTop,kBotT,tBot,TopInF,
     !              BotInF,lScreen,err)
        if(err.ne.0) goto (910,903) err
      end if

      if(lChem) then
        call ChemIn(lUpW,lTDep,NMat,NS,NSD,MaxItC,ChPar,TDep,kTopCh,
     !              cTop,kBotCh,cBot,epsi,tPulse,CumCh,cTolA,cTolR,
     !              lLinear,lEquil,lArtD,PeCr,lScreen,dSurf,cAtm,err)
        if(err.ne.0) goto (911,903) err
        do 8 i=1,NS
          write(ch1,'(i1)') i
          cName = '\solutex.out'
          cName(8:8) = ch1
          cFileName = cDataPath(1:iLengthPath)//cName
          open(80+i,file=cFileName, status='unknown',err=902)
8       continue
      end if
      if(TopInF.or.BotInF) then
        call SetBC(tMax,tAtm,rTop,rRoot,rBot,hCritA,hBot,hTop,GWL0L,
     !             TopInF,BotInF,cTop,cBot,NS,tTop,tBot,Ampl,lTemp,
     !             lChem,err,lMinStep)
        if(err.eq.1) goto 913
      end if
      if(SinkF) then
        call SinkIn(NMat,lChem,lMoSink,lSolRed,lSolAdd,P0,POptm,P2H,
     !              P2L,P3,r2H,r2L,aOsm,c50,P3c,NS,lScreen,err)
        if(err.ne.0) goto (909,903) err
        call SetSnk(NumNP,NMat,MatNum,x,hRoot,vRoot,Sink,rRoot,hNew,
     !              lMoSink,lSolRed,lSolAdd,P0,POptm,P2H,P2L,P3,r2H,r2L,
     !              aOsm,c50,P3c,Beta,lChem,NS,NSD,Conc,cRoot)
      end if
      close(30)
      close(32)
      close(50)

      call Profil (NumNP,NMat,x,MatNum,xSurf,Beta,Ah,AK,ATh,thr,ths,
     !             ConSat,hSat,lScreen,err)
      if(err.eq.1) goto 923
      call NodOut (NumNP,NMat,hNew,ThOld,Con,x,xSurf,CosAlf,tInit,
     !             MatNum,Cap,AK,Sink,ConSat,NS,NSD,Conc,TempO,Sorb,
     !             err)
      if(err.eq.1) goto 920
      call SubReg (NumNP,NMat,NLay,hNew,ThOld,ThOld,x,MatNum,LayNum,
     !             t-dt,dt,CosAlf,Con(1),Con(2),Con(NumNP),Con(NumNP-1),
     !             lChem,Conc,ChPar,0,ths,wCumT,wCumA,cCumT,cCumA,wVolI,
     !             cVolI,WatIn,SolIn,lWat,lTemp,TempN,TPar,TDep,NS,NSD,
     !             Sorb,lLinear,lEquil,err)
      if(err.eq.1) goto 921
      if(lChem.or.lTemp)
     !  call Veloc(NumNP,hOld,Con,x,CosAlf,vOld,ThOld(NumNP),
     !             ThOld(NumNP),Sink(NumNP),ThOld(1),ThOld(1),Sink(1),
     !             dt)
      if(.not.lWat) then
        do 11 i=1,NumNP
          vNew(i)=vOld(i)
          ThNew(i)=ThOld(i)
11      continue
      end if

      if(lScreen) write(*,*) 'beginning of numerical solution'

      call getdat (i,i,iday)

      call gettim (ihours,mins,isecs,i)
      Rtime1=iday*24.*60.*60.+ihours*60.*60.+mins*60.+isecs

*     Solve water movement
12    if(lWat) then
        call WatFlow(NumNP,NTab,NTabD,NMat,hTab,ConTab,CapTab,hNew,hOld,
     !             MatNum,ParD,ParW,Con,Cap,ConSat,Ah,AK,ATh,hSat,hTemp,
     !               KodTop,KodBot,rTop,rBot,CosAlf,t,dt,x,Sink,P,R,S,
     !               FreeD,SeepF,qGWLF,Aqh,Bqh,GWL0L,hTop,hBot,hCritA,
     !               hCritS,WLayer,IterW,ItCum,TopInf,kTOld,kBOld,TolTh,
     !               TolH,MaxIt,dtMin,tOld,dtOpt,ConvgF,TheTab,ThNew,
     !               ThOld,thr,ths,lWDep,TempN,Kappa,KappaO,AThS,ThRR,
     !               ConO,ConR,AKS,AhW,AThW,AKW,iHyst,iModel)
cyap++
	 if(.not.ConvgF) goto 2002
cyap--
      else
        iterW=1
        ItCum=ItCum+1
      end if

*     To calculate the velocities
      if(lWat.and.(lTemp.or.lChem))
     !  call Veloc(NumNP,hNew,Con,x,CosAlf,vNew,ThNew(NumNP),
     !             ThOld(NumNP),Sink(NumNP),ThOld(1),ThNew(1),Sink(1),
     !             dt)

*     Root zone calculations
      if(lRoot)
     !  call SetRG(NumNP,x,Beta,t,tRMin,tRHarv,xRMin,xRMax,RGR)
      if(SinkF)
     !  call SetSnk(NumNP,NMat,MatNum,x,hRoot,vRoot,Sink,rRoot,hNew,
     !              lMoSink,lSolRed,lSolAdd,P0,POptm,P2H,P2L,P3,r2H,
     !              r2L,aOsm,c50,P3c,Beta,lChem,NS,NSD,Conc,cRoot)

*     Calculation of heat transport
      if(lTemp)
     !  call Temper(NumNP,NMat,x,dt,t,MatNum,TempO,TempN,TPar,Ampl,P,R,
     !              S,Q,vOld,vNew,ThOld,ThNew,Retard,Disp,Sink,tPeriod,
     !              kTopT,tTop,kBotT,tBot)

*     Calculations of the solute transport
      if(lChem) then
        call Solute(NumNP,NMat,NS,NSD,x,dt,t,tPulse,ChPar,MatNum,ThOld,
     !              ThNew,vOld,vNew,Disp,epsi,kTopCh,cTop,kBotCh,cBot,
     !              Conc,P,R,S,Q,g0,g1,Retard,cvTop,cvBot,cvCh0,cvCh1,
     !              lUpW,wc,Peclet,Courant,dtMaxC,TempO,TempN,cNew,
     !              cPrevO,cTemp,TDep,ths,cTolA,cTolR,IterC,MaxItC,
     !              hTemp,Sorb,SorbN,lLinear,lEquil,lArtD,PeCr,q0,q1,
     !              dSurf,cAtm,err)
        if(err.ne.0) goto 927
      end if

*     T-level information
      call TLInf(Con(1),Con(2),Con(NumNP),Con(NumNP-1),x(1),x(2),
     !           x(NumNP),x(NumNP-1),CosAlf,t,dt,IterW,IterC,TLevel,
     !           ShortO,TPrint(PLevel),rTop,rRoot,vRoot,hNew(NumNP),
     !           hNew(NumNP-1),hRoot,hNew(1),hNew(2),CumQ,ItCum,KodTop,
     !           KodBot,ConvgF,lWat,lChem,cRoot,NS,NSD,NumNP,Conc,cvTop,
     !           cvBot,cvCh0,cvCh1,Peclet,Courant,wCumT,wCumA,cCumT,
     !           cCumA,CumCh,ThNew(NumNP),ThOld(NumNP),Sink(NumNP),
     !           ThNew(1),ThOld(1),Sink(1),lScreen,err)
      if(err.ne.0) goto (919,918,926) ierr
      if(NObs.gt.0) then
        call ObsNod(t,NumNP,NObs,NS,NSD,Node,Conc,hNew,ThNew,TempN,
     !              lChem,err)
        if(err.eq.1) goto 922
      end if

*     P-level information
      if(abs(TPrint(PLevel)-t).lt.0.001*dt) then
        call NodOut(NumNP,NMat,hNew,ThNew,Con,x,xSurf,CosAlf,
     !              TPrint(PLevel),MatNum,Cap,AK,Sink,ConSat,NS,NSD,
     !              Conc,TempN,Sorb,err)
        if(err.eq.1) goto 920
        call SubReg(NumNP,NMat,NLay,hNew,ThNew,ThOld,x,MatNum,LayNum,t,
     !              dt,CosAlf,Con(1),Con(2),Con(NumNP),Con(NumNP-1),
     !              lChem,Conc,ChPar,PLevel,ths,wCumT,wCumA,cCumT,cCumA,
     !              wVolI,cVolI,WatIn,SolIn,lWat,lTemp,TempN,TPar,TDep,
     !              NS,NSD,Sorb,lLinear,lEquil,err)
        if(err.eq.1) goto 921
        PLevel=PLevel+1
      end if

*     A-level information
      if(abs(t-tAtm).le.0.001*dt.and.(TopInF.or.BotInF)) then
        call ALInf(t,CumQ,hNew(NumNP),hRoot,hNew(1),ALevel,err)
        if(err.ne.0) goto (925) ierr
        call SetBC(tMax,tAtm,rTop,rRoot,rBot,hCritA,hBot,hTop,GWL0L,
     !             TopInF,BotInF,cTop,cBot,NS,tTop,tBot,Ampl,lTemp,
     !             lChem,err,lMinStep)
        if(err.eq.1) goto 913
        ALevel=ALevel+1
      end if

*     Time governing
      if(abs(t-tMax).le.0.001*dt) then
        write(70,'(''end'')')
        Rtime2=iday*24.*60.*60.+ihours*60.*60.+mins*60.+isecs
        write(70,*)
        write(70,*) 'Real time [sec]',Rtime2-RTime1
        write(71,'(''end'')')
        write(72,'(''end'')')
        write(77,'(''end'')')
        if(lChem) then
          do 13 jj=1,NS
            write(80+jj,'(''end'')')
13        continue
        end if
	  close(31)
	  goto 2001
C        write(*,*) 'Press Enter to continue'
C        read(*,*)
C        stop
      else
        tOld=t
        dtOld=dt
        kTOld=KodTop
        kBOld=KodBot
        Iter=max0(IterW,IterC)
        call TmCont(dt,dtMax,dtOpt,DMul,DMul2,dtMin,Iter,TPrint(PLevel),
     !              tAtm,t,tMax,dtMaxC,ItMin,ItMax,lMinStep,dtInit)
        t=t+dt
        TLevel=TLevel+1
      end if

*     New pressure heads
      if(lWat) then
        do 14 i=1,NumNP
          hTemp(i)=hNew(i)+(hNew(i)-hOld(i))*dt/dtOld
          hOld(i) =hNew(i)
          hNew(i) =hTemp(i)
          ThOld(i)=ThNew(i)
          if(lTemp.or.lChem) vOld(i)=vNew(i)
14      continue
      end if

      goto 12

* --- End of time loop -------------------------------------------------

*     Error messages
901   ierr=1
      goto 1000
902   ierr=2
      goto 1000
903   ierr=3
      goto 1000
904   ierr=4
      goto 1000
905   ierr=5
      goto 1000
906   ierr=6
      goto 1000
907   ierr=7
      goto 1000
908   ierr=8
      goto 1000
909   ierr=9
      goto 1000
910   ierr=10
      goto 1000
911   ierr=11
      goto 1000
912   ierr=12
      goto 1000
913   ierr=13
      goto 1000
914   ierr=14
      goto 1000
915   ierr=15
      goto 1000
916   ierr=16
      goto 1000
c917   ierr=17
c      goto 1000
918   ierr=18
      goto 1000
919   ierr=19
      goto 1000
920   ierr=20
      goto 1000
921   ierr=21
      goto 1000
922   ierr=22
      goto 1000
923   ierr=23
      goto 1000
c924   ierr=24
c      goto 1000
925   ierr=25
      goto 1000
926   ierr=26
      goto 1000
927   ierr=27
      goto 1000
928   ierr=28
      goto 1000
929   ierr=29
      goto 1000

1000  call ErrorOut(ierr,cFileName,cDataPath,iLengthPath,lScreen)
	ptferr=1 !GM
      write(*,*) 'Press Enter to continue'
C      read(*,*)
C      stop

101   format(a)
2002	ptferr=1 !GM
2001  CONTINUE
	close(50)
	close(70)
	close(71)
	close(72)
      close(75)
	close(76)
	close(77)
	close(78)
	close(26)
	close(31)
      RETURN
      end

************************************************************************

      subroutine ErrorOut(ierr,cFileName,cDataPath,iLengthPath,lScreen)

      character*90 cErr(30),cFileNameErr
      Character*80 cFileName,cDataPath
      logical lScreen

      cErr( 1)='Open file error in file :'
      cErr( 2)='File already exists or hard disk is full ! Open file err
     !or in output file : '
      cErr( 3)='Error when writing to an output file !'
      cErr( 4)='Error when reading from an input file Level_01.dir data
     !pathway !'
      cErr( 5)='Error when reading from an input file Selector.in Basic
     !Informations !'
      cErr( 6)='Error when reading from an input file Selector.in Water
     !Flow Informations !'
      cErr( 7)='Error when reading from an input file Selector.in Time I
     !nformations !'
      cErr( 8)='Error when reading from an input file Selector.in Root G
     !rowth Informations !'
      cErr( 9)='Error when reading from an input file Selector.in Sink I
     !nformations !'
      cErr(10)='Error when reading from an input file Selector.in Heat T
     !ransport Informations !'
      cErr(11)='Error when reading from an input file Selector.in Solute
     ! Transport Informations !'
      cErr(12)='Error when reading from an input file Profile.dat !'
      cErr(13)='Error when reading from an input file Atmosph.in !'
      cErr(14)='Dimension in NumNPD is exceeded !'
      cErr(15)='Dimension in NObsD is exceeded !'
      cErr(16)='Dimension in NMatD is exceeded !'
      cErr(17)='Error when writing into an output file I_CHECK.OUT !'
      cErr(18)='Error when writing into an output file RUN_INF.OUT !'
      cErr(19)='Error when writing into an output file T_LEVEL.OUT !'
      cErr(20)='Error when writing into an output file NOD_INF.OUT !'
      cErr(21)='Error when writing into an output file BALANCE.OUT !'
      cErr(22)='Error when writing into an output file OBS_NODE.OUT !'
      cErr(23)='Error when writing into an output file PROFILE.OUT !'
      cErr(24)='Error when writing into an output file SOLID.OUT !'
      cErr(25)='Error when writing into an output file A_LEVEL.OUT !'
      cErr(26)='Error when writing into an output file SOLUTE.OUT !'
      cErr(27)='Does not converge in the solute transport module !'
      cErr(28)='Number of Print-Times is exceeded !'
      cErr(29)='Dimension in NSD is exceeded !'

      cFileNameErr = cDataPath(1:iLengthPath)//'\Error.msg'
      open(99,file=cFileNameErr,status='unknown')
      if(ierr.le.2) then
        if(lScreen) write( *,*) cErr(ierr),cFileName
        write(99,*) cErr(ierr),cFileName
      else
        if(lScreen) write( *,*) cErr(ierr)
        write(99,*) cErr(ierr)
      end if
      close(99)

      return
      end

* ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||