      program main
c*******************************************
c
c
c                   DYME-2D 
c    (A 2D Dynamic Meshfree (explicit) Program)
c
c    -------- h-adaptive Version ------
c                (Plane Strain)
c  
c  This is a 2-D FEM/RKPM Fortarn Program 
c  by using explicit time integration.
c  The code is designated as a sovler for
c  Large deformation problems in solid mechanics.
c
c  Currently, it can deal with
c
c  (1) Hyperelastic materials (Rubber);
c  (2) Elastic-plastic materials;
c  (3) Elasto-viscoplastic materials;
c  (4) Thermo-elasto-viscoplastic material;
c
c   Quadrilateral & Triangular Mesh 
c   (Triangular Background Cell has not been implemented yet)
c  
c  Three time integration scheme is implemented:
c
c   1. Predictor-Corrector Integrator via Newmark-beta;
c   2. Central Difference;
c   3. Improved Euler-theta;
c
c
c  A total Lagrangian formulation is implemented.
c 
c
c
c  Author: Shaofan Li
c  Copyright @
c  Northwestern University,  June, 1998 -- Oct, 1999
c
c
c 
ci***************************************************************
c
      implicit none
      include 'parameter.h'
c
      real*8 xs(2,maxNumnp),xsk(2,maxGP),
     &       xm(2,maxNumnp),xmk(2,maxGP),
     &       dxm(2,maxNumnp),
     &       dvm(maxNumnp),dvmk(maxGP)
c  
      integer LmDispbcX(maxDispbcX),LmDispbcY(maxDispbcY),
     &        LmVelbcX(maxVelbcX),LmVelbcY(maxVelbcY),
     &        LmTract(2,maxTract)
c
      real*8  vDispbcX(maxDispbcX),vDispbcY(maxDispbcY),
     &        vVelbcX(maxVelbcX),vVelbcY(maxVelbcY),
     &        vTract(2,2,maxTract)
c
      integer igp_elemind(maxGP),igp_locintInd(maxGP)
      real*8  gp_loc2D(2,maxIntElem)
c
      integer  Lnp(maxNumnp),Lmnp(mnsch,maxNumnp),
     &         Lgp(maxGP),Lmgp(mnsch,maxGP),
     &         lnods(maxNode,maxElem),
     &         lgops(maxIntElem,maxElem),gpelem(maxGP),
     &         Leff_adp(maxNumnp)
c
      real*8   shpk(mnsch,maxGP),shpkdx(mnsch,maxGP),
     &         shpkdy(mnsch,maxGP),shpn(mnsch,maxNumnp),
     &         shBCX(maxEssbcX,maxEssbcX),
     &         shBCY(maxEssbcY,maxEssbcY)
      real*8   fint(2,maxNumnp),fext(2,maxNumnp),press(maxGP)
      integer  nxtime(mxndt),nytime(mxndt)
c
      real*8   disp(2,maxNumnp),vel(2,maxNumnp),acc(2,maxNumnp),
     &         disp_last(2,maxNumnp),vel_last(2,maxNumnp),
     &         acc_last(2,maxNumnp)
      real*8   disp_node(2,maxNumnp),d_disp(2,maxNumnp),
     &         vel_node(2,maxNumnp),
     &         acc_node(2,maxNumnp)
      real*8   amass(maxNumnp)
      real*8   sts_CHY(3,3),	! Cauchy stress for EP material
     &         sts_PK1(3,3),	! PK1
     &         sts_KH(3,3)
c
      character*60 fhead,frname,fwname
      character*60 tmpstr1,tmpstr2,tmpstr3
      character*26 suffix
      character*1  single
c
      real*8  ys_k(maxGP),
     &        stsgp(4,maxGP),effstsgp(maxGP),
     &        effstegp(maxGP),thermo(maxGP),
     &        effRgp(maxGP)
c
      integer iextf_LoadingType,ivel_loadingType
      integer i_ask_fhead  !=1 (ask for fhead anyway )
      			   !=0 (can be input as command parameter)
      integer i_rate       ! For elasto-plastic problem:
                           !    i_rate = 0 : Jaumann ; 
		           !    i_rate = 1 : Green-Naghdi ;
		           ! For elsato-viscoplastic problem:
		           !    i_rate = 0 : Jaumann-Cauchy;
		           !    i_rate = 1 : Jaumann-Kirchhoff; 
c
      integer ifwtimestep,ifwdispt,ifwvelpt,ifwstspt,
     &        ifwecho,
     &        ifwload,ifwvelload,ifwadapt,ifwout,
     &        ifwmesh,ifwgmesh
c
      integer numnp,mgk,nelem,
     &        nndxdt,nndydt,ioutput,istep_debug,itmp
c
      integer i,ie,ik,inode,ip,ipt,istep,ists,
     &        idim,i_dtype,
     &        j,jp,jpt,jLoop,mloop,
     &        k,kp,kq,kr,n_zero,
     &        iz,ist,ike,ik_eID
c
      real*8 eps,teps,elminx,elminy,elengx,elengy,
     &       dt00,dvtemp,
     &       dixxx,diyyy,
     &       dttem,
     &       xle00,
     &       temp,
     &       tom,v_istep,effsts,effep,xkc1,xkc2,
     &       alpha,alpha0,alpha1
c
       real*8 effSnpt(maxNumnp),effEnpt(maxNumnp),
     &        effEopt(maxNumnp),stsnpt(4,maxNumnp),
     &        effTemp(maxNumnp),effRate(maxNumnp),
     &        effMax(maxNumnp)
c
      real*8  ax,ay,afact1,afact2,
     &        dt,time,
     &        rho0,cc1,cc2,rambda,
     &        el_E, el_V, el_G, el_K,
     &        pl_k0, pl_H, pl_EP,
     &        r0,xc1,xc2,zimp,cr1,cr2,
     &        beta,gamma,theta,
     &        Temp_0,Cp
c
      integer nstep,nOutputFrq,nPlotFrq,iplot,iplot_gp,
     &        iMaterType,nnode,nintElem,iLumping,
     &        nnc1,nnc2,
     &        imethRKPM,imethFEM,imeth,iIntMeth,iGL,
     &        iInter,iAdapt,numadp,iad_stage,Istage,
     &        nTract,nDispbcX,nDispbcY,nVelbcX,nVelbcY
c
       integer nTractf,LmTractF(maxTF),
     &         Idamage_mark(maxGP),idamage,
     &         ifail
       real*8  TF,BW,DFL,bweight(maxTF),
     &         dispar,velpar,accpar,
     &         test_m1,test_m2,test_m3,
     &         temp_max,test_min(3)
c
       integer iCrack
       real*8  Cyplane,Cxmin,Cxmax
c
      common /shapeBC/shBCX,shBCY
      common /shapeK/shpk,shpkdx,shpkdy,shpn
      common /connect/Lnp,Lmnp,Lgp,Lmgp,lnods
      common /connectK/lgops,gpelem
      common /mesh/xm,dxm,dvm,numnp
      common /crack/Cyplane,Cxmin,Cxmax
c
      common /rkpm/ax,ay,afact1,afact2,nnc1,nnc2
      common /stepT/dt,time
      common /hyperelast/rho0,cc1,cc2,rambda      
      common /elastic/ el_E, el_V, el_G, el_K
      common /plastic/ pl_k0, pl_H, pl_EP
      common /temperature/Temp_0,Cp
      common /shear/r0,xc1,xc2,zimp,cr1,cr2
      common /algorithm/beta,gamma,theta
c
      common /step/nstep,nOutputFrq,nPlotFrq,iplot,iplot_gp
      common /ctrl/iMaterType,nnode,nintElem,iLumping
      common /meth/imethRKPM,imethFEM,imeth,iIntMeth,iInter
      common /adapt/iGL,iAdapt,numadp,iad_stage,Istage
      common /bound/nTract,nDispbcX,nDispbcY,nVelbcX,nVelbcY
      common /force/nTractf,LmTractF
      common /output/stsnpt,effEnpt,effSnpt,effTemp,effRate,effMax
      common /Type1/iCrack
      common /Type2/bweight
c
      !** iIntMeth = 1 : Predictor-corrector (Newmark-beta);
      !            = 2 : Central difference;
      !            = 3 : Improved Euler-theta;
      !** Parameters: 
      !      beta, gamma: Newmark-beta
      !      theta      : Improved Euler-theta
c
c
      i_ask_fhead=1
c
	! fhead is input from the command line
      if (i_ask_fhead .eq. 1) then
         write(*,*) 'please input the head of the files name :'
         read(*,'(a)') fhead
      endif

      call appext(fhead,'inp',frname)
      write(*,*) 'frname=',frname
      open(10,file =frname)
      
      ifwTimeStep=16
      
      call appext(fhead,'timestep',fwname)
      open(ifwTimeStep,file =fwname)
      
      ifwDisPT=17
      call appext(fhead,'dispt',fwname)
      open(ifwDisPT,file =fwname)
      
      ifwVelPT=18
      call appext(fhead,'velpt',fwname)
      open(ifwVelPT,file =fwname)
      
      ifwStsPT=19
      call appext(fhead,'stspt',fwname)
      open(ifwStsPT,file =fwname)
      
      ifwGmesh=24
      call appext(fhead,'gmesh',fwname)
      open(ifwGmesh,file =fwname)
      
      ifwMesh=25
      call appext(fhead,'mesh',fwname)
      open(ifwMesh,file =fwname)
      
      ifwEcho=26
      call appext(fhead,'echo',fwname)
      open(ifwEcho,file =fwname)
      
      ifwLoad=27
      call appext(fhead,'load',fwname)
      open(ifwLoad,file =fwname)
      
      ifwVelLoad=28
      call appext(fhead,'velload',fwname)
      open(ifwVelLoad,file =fwname)

      ifwAdapt =29
c      call appext(fhead,'adapt',fwname)
c      open(ifwAdapt,file =fwname)
c            
c-----------------------------------------------------
c     
      call prepro(xm,dxm,dvm,xmk,dvmk,nxtime,nytime,
     &            LmTract,vTract,
     &            LmDispbcX,LmDispbcY,vDispbcX,vDispbcY,
     &            LmVelbcX,LmVelbcY,vVelbcX,vVelbcY,
     &            lnods,igp_elemind,igp_locintInd,gp_loc2D,     
     &            nelem,numnp,mgk,nndxdt,nndydt,
     &            iextf_LoadingType,ivel_loadingType,ifwEcho)
c
      write(*,*) 'prepro passed.'
c
      temp   = 0.0d0
      Istage = 0
      teps   = 1.99999999999d0
      eps    = 1.0d-20
      suffix = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
c
      if(imeth .eq. imethRKPM) then
         call assign2d_dxm(dxm,xm,lnods,nelem,numnp)
      endif
c
c---------------------------------------------
c
 38   continue  ! the conjunction point
c
c--------------------------------------------
c
      do ipt = 1, numnp
        do j = 1, 2
          xs(j,ipt) = xm(j,ipt)
        enddo
      enddo
c      
      do ik = 1, mgk
        xsk(1,ik) = xmk(1,ik)
        xsk(2,ik) = xmk(2,ik)
      enddo
c
c*** below it will do :
c   For FEM or RKPM method:
c     (1) search nodes for each Gaussian point
c	  generate the lists: Lgp(ik) and Lmgp(ip,ik)
c                             where  ip=1..Lgp(ik)
c                                    ik=1..mgk
c
c     (2) search nodes for each node
c	  generate the lists: Lnp(ip) and Lmnp(ip,jp)
c                             where  jp=1..Lnp(ip)
c                                    ip=1..numnp
c
c     (3) cal the values of shape functions at each node and each
c         Gaussian point
c         generate the arrays: 
c          shpn(jp,ip) -- shape function N_I (x_J) at each node x_J
c                                where ip=1..Lnp(jp)
c                                      jp=1..numnp
c          shpk(ik,ip) -- shape function N_I (x_K) at each gauss pt x_K
c                                where ip=1..Lgp(ik)
c                                      ik=1..mgk
c
c
c-----------------------------------------------------
c     Searching nodes for each gauss point
c-----------------------------------------------------
c
      call search2d(xmk,mgk,igp_elemind,
     &              igp_LocIntInd,gp_loc2d)
c
c
c-----------------------------------------------------------------
c     mass matrix and lumped nodal volume calculation 
c     
      dvtemp = 0.0d0
      do ipt = 1, numnp
         amass(ipt) = 0.0d0
      enddo
c
      alpha  = 0.0d0
      alpha0 = 0.0d0
      alpha1 = 0.0d0
      do ik = 1, mgk
	 jLoop = Lgp(ik)
	 do jp = 1, jLoop
	    alpha0 = alpha0 + rho0*dvmk(ik)*shpk(jp,ik)**2
	    alpha1 = alpha1 + rho0*dvmk(ik)
	    dvtemp = dvtemp + dvmk(ik)*shpk(jp,ik)
	 enddo ! jp
      enddo    ! ik
c
      alpha = alpha1/alpha0
      print *, 'new mgk', mgk
      do ik = 1, mgk
         jLoop = Lgp(ik)
         do jp = 1, jLoop
            jpt = Lmgp(jp,ik)
c
	    if (iLumping .eq. 0) then
              amass(jpt) = amass(jpt) +
     &                     rho0*dvmk(ik)*shpk(jp,ik)
c
              if(dvmk(ik) .le. eps ) then
		 print *, 'wrong dvmk(ik)', dvmk(ik) 
	      endif
c
           elseif(iLumping .eq. 1) then
              amass(jpt) = amass(jpt) + alpha *
     &                  rho0*dvmk(ik)*shpk(jp,ik)**2
	   endif ! iLumping
         enddo   ! jp
      enddo      ! ik
c
      write(*,*) 'dvtemp=',dvtemp
      write(*,*) 'mass is done'
c      
      tom = 0.0d0
      do ipt = 1, numnp
         tom = tom + amass(ipt)
         if(amass(ipt) .le. eps) then
            print *,ipt,'  singularity in M'
            stop
         endif
         amass(ipt) = 1.0d0/amass(ipt)
      enddo
      write(*,*) 'tom=',tom
c
c
      if (imeth.eq.imethRKPM) then
          call bcinv(LmDispbcX,LmDispbcY,
     &           LmVelbcX,LmVelbcY,
     &           xm,dxm,dvm,
     &           numnp)
      endif
c
c
		!*** initial condition
c
		! Initial displacement fields
      do ipt = 1, numnp
        do idim = 1, 2
          fext(idim,ipt)      = 0.0d0
          fint(idim,ipt)      = 0.0d0
          disp(idim,ipt)      = 0.0d0
          disp_last(idim,ipt) = 0.0d0
          vel(idim,ipt)       = 0.0d0
          acc(idim,ipt)       = 0.0d0
        enddo
      enddo
c                
		! Initial stress states
      do ik=1,mgk
         do ists=1,4
            stsgp(ists,ik)= 0.0d0
         enddo
	 effstsgp(ik) = 0.0d0
	 effstegp(ik) = 0.0d0
      enddo
c
c                
c                !  Initial Traction Force
c
      istep=0
      call trabc(LmTract,vTract,xm,dxm,dvm,fext,nTract,
     &           numnp,
     &           istep,nstep,iextf_LoadingType,ifwLoad)
                !
		! Initial Accel
                ! fint=0 at time0	
      do ipt = 1, numnp
        acc(1,ipt) = amass(ipt)*fext(1,ipt)
        acc(2,ipt) = amass(ipt)*fext(2,ipt)
      enddo
c
		! init elast-plast parameters
      if (iMaterType.eq.2) then
          do ik=1,mgk
             effstegp(ik) = 0.0d0
             ys_k(ik)     = pl_k0
          end do
      elseif (iMaterType .eq. 4) then
	  do ik = 1,mgk
	     thermo(ik) = Temp_0
	  enddo
      elseif (iMaterType .eq. 5) then
	  do ik = 1,mgk
	     thermo(ik)       = Temp_0
	     Idamage_mark(ik) = 0
	     effRgp(ik)       = 0.0
             ys_k(ik)         = 0.0
	  enddo
      else
      endif
c
      Temp_max = Temp_0
c
      test_m1 = 0.0
      test_m2 = 0.0
      test_m3 = 0.0
c
c..........! output initial location
c
      if (1 .eq. 0) then          
         ifwout=31
         n_zero=0
         call Int2Str(n_zero,tmpstr1)
         call appstr('.newxy',tmpstr1,tmpstr2)
         call appstr(fhead,tmpstr2,fwname)
         write(*,*) 'fwname=',fwname
         open(ifwout,file=fwname)
c
         do ip = 1, numnp
            write(ifwout,'(1x,2e15.6)') xs(1,ip),xs(2,ip)
         enddo
         close(ifwout)
      endif
c
      ! output initial location in TECPLOT format
      if (Istage .eq.  0) then          
         ifwout=31
         n_zero=0
         call Int2Str(n_zero,tmpstr1)
         call appstr('.tecplt',tmpstr1,tmpstr2)
         call appstr(fhead,tmpstr2,fwname)
         write(*,*) 'fwname=',fwname
         open(ifwout,file=fwname)
c
         write(ifwout,'(1x,a)')  'VARIABLES="X","Y","S1","E1"'
         if (nnode.eq.4) then
            write(ifwout,'(1x,a,i6,a,i6,a)')
     &            'ZONE N =',numnp,', E =',nelem,
     &            'F = FEPOINT, ET = QUADRILATERAL'
         elseif (nnode.eq.3) then
            write(ifwout,'(1x,a,i6,a,i6,a)')
     &            'ZONE N =',numnp,', E =',nelem,
     &            'F = FEPOINT, ET = TRIANGLE' 
         endif
         do ip=1,numnp
            write(ifwout,'(1x,4e16.7)') xs(1,ip),xs(2,ip),temp,temp
         enddo
         do ie=1,nelem
            write(ifwout,'(1x,4i8)') (lnods(inode,ie),inode=1,nnode)
         enddo
               
         close(ifwout)
      endif  ! Istage = 0
c
      print *,'Istage =',Istage
c
c**************************************************************
c
c      time integration begin
c
c**************************************************************
c
      time    = 0.00
      iOutput = 0       ! the extersion of output file name
c
c
      istep_debug=1

      do 100 istep = 1, nstep
      
          time=time+dt
      
c         write(*,*) 'istep=',istep
c
         if (istep.ge.istep_debug) then
            itmp=itmp
         endif
c
c--------------------------------------------------------------
c     Flexible Time Step Calculation
c
c
c		! !!!!!!!!!!! ignore this part
c
      if (1.eq.0) then
         if (istep.eq.1) dt00 = dt
         elminx = 1.0e5
         elminy = 1.0e5

         do kp = 1, nndxdt-1
           kq = nxtime(kp)
           kr = nxtime(kp+1)
           elengx = dabs(xs(1,kr)-xs(1,kq))
           elminx = dmin1(elminx,elengx)
         enddo

         do kp = 1, nndydt-1
           kq = nytime(kp)
           kr = nytime(kp+1)
           elengy = dabs(xs(2,kr)-xs(2,kq))
           elminy = dmin1(elminy,elengy)
         enddo

         xle00 = dabs(xm(1,nxtime(1))-xm(1,nxtime(2)))

         dttem = elengx*dt00/xle00

         dt = dmin1(dttem,dt)
      endif		!endif ignore this part
      
c
c---------------------------------------------------------------
c     Predictor Phase
c

      if (iIntMeth.eq.1) then         
                      !** predictor-corrector
                      
         do ip = 1, numnp
            do idim = 1,2
               vel(idim,ip)  = vel(idim,ip)+(1.-gamma)*dt*acc(idim,ip)
               disp(idim,ip) = disp(idim,ip) + dt*vel(idim,ip)
     &                + dt**2*(0.5-beta)*acc(idim,ip)
            enddo
         enddo
      elseif (iIntMeth.eq.2)  then
c
c...............!** central difference
c         Find d_{theta}, vel_{theta} 
c         (The code is for forward tangent method)
c 
         if (istep.ne.1) then
            do ip=1,numnp
               do idim=1,2
                  vel(idim,ip)=vel(idim,ip)+dt*acc(idim,ip)
                  disp(idim,ip)=disp(idim,ip)+dt*vel(idim,ip)
               enddo
            enddo
         else      ! the first step, only dt/2 when eval veloc
            do ip=1,numnp
               do idim=1,2
                  vel(idim,ip) =vel(idim,ip) +0.5*dt*acc(idim,ip)
               enddo
            enddo
         endif
c      
         v_istep = real(istep) + 0.5d0
         call bcvel(LmDispbcX,LmDispbcY,
     &              vDispbcX,vDispbcY,
     &              LmVelbcX,LmVelbcY,
     &              vVelbcX,vVelbcY,
     &              vel,shpn,
     &              v_istep,nstep,ivel_LoadingType,
     &              ifwVelLoad)
c
            do ip=1,numnp
               do idim=1,2
                  disp(idim,ip)=disp(idim,ip)+0.5*dt*vel(idim,ip)
               enddo
            enddo
c
c         
      elseif (iIntMeth .eq. 3) then
c
c................! Improved Euler-theta method ! .........
c## Note: ##
c
c  The improved Euler-theta method is explained by B. Moran
c  in {Computers & Structures} Vol. 27, No. 2 pp. 241-247
c  (1987), which is suitable for the time integration for
c  elasto-viscoplasticity materials.
c
c  Note: here vel := v_{theta}
c
         do ip = 1, numnp
	    do idim = 1, 2
	       vel(idim,ip) = vel_last(idim,ip)
     &                      + dt*theta*acc_last(idim,ip) ! Eq.(*)
	    enddo
	 enddo 
c         
         v_istep = real(istep) + theta
         call bcvel(LmDispbcX,LmDispbcY,
     &              vDispbcX,vDispbcY,
     &              LmVelbcX,LmVelbcY,
     &              vVelbcX,vVelbcY,
     &              vel,shpn,
     &              v_istep,nstep,
     &              ivel_LoadingType,
     &              ifwVelLoad)
c         
c.........while the predict value of disp_(n+1) is calculated
c         as follows:
c
c               disp_{n+1} =disp_n + dt*vel_{theta}   (*)
c               substituting (*) inot 
c               disp_{theta} = (1-theta) d_n + theta d_{n+1} 
c     yields 
c               disp_{theta} = d_n + theta * dt * v_{theta}
c
c
         do ip=1,numnp
            do idim=1,2                    
               disp(idim,ip) = disp_last(idim,ip)
     &                       + theta*dt*vel(idim,ip)
            enddo
         enddo
c         
c..........the difference is hadapt-wd
c
c         v_istep=real(istep) + 0.5
c         call bcdisp(LmDispbcX,LmDispbcY,
c     &               vDispbcX,vDispbcY,
c     &               LmVelbcX,LmVelbcY,
c     &               vVelbcX,vVelbcY,
c     &               disp,vel,shpn,
c     &               v_istep,nstep,ivel_LoadingType)
c         
c         
c#Note# The predicted velocity array, vel, is: vel_theta , which is defined as
c       vel_theta := theta * vel_{n+1} + (1 - theta ) * vel_{n}  ! Eq. (**) 
c       if we use Euler forward scheme as predictor, then
c            vel_{p} = vel_n + dt * acc_n
c      substituting Eq.(**) into Eq.(*) yields
c            vel_theta = vel_n + theta * dt * acc_n
c
      else
         write(*,*) 'wrong number of iIntMeth=',iIntMeth
         stop         
      endif ! iIntMeth
c      
c
c      iUpdate_velBCnode=1
c
c		! d_disp is required to eval the increment 
c               ! of stresses for elast-plast material      
      do ip=1,numnp
         do idim=1,2
            d_disp(idim,ip) = disp(idim,ip) - disp_last(idim,ip)
         enddo
      enddo
c
c---------------------------------------------------------
c Calculate the stresses and evaluate the fint vector
c---------------------------------------------------------
c
      do ip=1,numnp
         do idim=1,2
            fint(idim,ip) = 0.0d0
         enddo
      enddo
c
      do ik=1,mgk
c       
         if ( iMaterType.eq.1 ) then
			!** hyperelastic material         
c         
            call calsts_hyp(xmk,disp,
     &                      sts_PK1,effsts,
     &                      mgk,ik)
c
            press(ik)=(sts_PK1(1,1)+sts_PK1(2,2)+sts_PK1(3,3))/3.
            effstsgp(ik) = effsts
            stsgp(1,ik)  = sts_PK1(1,1)
            stsgp(2,ik)  = sts_PK1(2,2)
            stsgp(3,ik)  = sts_PK1(3,3)
            stsgp(4,ik)  = sts_PK1(1,2) ! PK1 is not symmetric
c            
         elseif (iMaterType.eq.2) then
           !
           !** elastic-plastic material
           ! recover the stress values of last time step
           !
            do i=1,3
               do j=1,3
                  sts_CHY(i,j) = 0.0d0
               enddo
            enddo
            sts_CHY(1,1) = stsgp(1,ik)
            sts_CHY(2,2) = stsgp(2,ik)
            sts_CHY(3,3) = stsgp(3,ik)
            sts_CHY(1,2) = stsgp(4,ik)
            sts_CHY(2,1) = sts_CHY(1,2)
c         
            call calsts_ep(disp,d_disp,
     &                     sts_CHY,sts_PK1,
     &                     effstegp(ik),ys_k(ik),
     &                     effstsgp(ik),
     &                     mgk,ik)
c            
            press(ik)=(sts_CHY(1,1)+sts_CHY(2,2)+sts_CHY(3,3))/3.
c
            stsgp(1,ik)  = sts_CHY(1,1)
            stsgp(2,ik)  = sts_CHY(2,2)
            stsgp(3,ik)  = sts_CHY(3,3)
            stsgp(4,ik)  = sts_CHY(1,2)
c
         elseif (iMaterType .eq. 3) then
c
c............. Elasto-viscoplasticity ............... 
c
c  	! recover the stress values of last time step
c
            do i=1,3
               do j=1,3
                  sts_CHY(i,j) = 0.0d0
               enddo
            enddo
c
            sts_CHY(1,1) = stsgp(1,ik)
            sts_CHY(2,2) = stsgp(2,ik)
            sts_CHY(3,3) = stsgp(3,ik)
            sts_CHY(1,2) = stsgp(4,ik)
            sts_CHY(2,1) = sts_CHY(1,2)
c
c  XXXX eff_sts ?  eff_ste ?
c
            xkc1    = xsk(1,ik)
            xkc2    = xsk(2,ik)
            i_rate  = 0		! via Cauchy stress
            i_dtype = 0		! large deformation
            effsts  = effstsgp(ik)
            effep   = effstegp(ik)
c
	    call calsts_evp(disp,vel,
     &            sts_CHY,sts_PK1,
     &            effep,effsts,dt,theta,
     &            mgk,ik,
     &            xkc1,xkc2,i_rate,i_dtype)
c 
c
            effstsgp(ik) = effsts
            effstegp(ik) = effep
            stsgp(1,ik)  = sts_CHY(1,1)
            stsgp(2,ik)  = sts_CHY(2,2)
            stsgp(3,ik)  = sts_CHY(3,3)
            stsgp(4,ik)  = sts_CHY(1,2)
c            
c##
c
	    if (effsts .lt. 0.0 ) then
		print *, 'Effective Stress Error!!', effsts,effep
		stop
	    endif
c             
         elseif (iMaterType .eq. 4) then
c
c........Thermo-elasto-viscoplastic material (Adiabatic) 
c
            do i=1,3
	       do j=1,3
                  sts_KH(i,j) = 0.0
               enddo
            enddo
c
            xkc1    = xsk(1,ik)
            xkc2    = xsk(2,ik)
            sts_KH(1,1) = stsgp(1,ik)
	    sts_KH(2,2) = stsgp(2,ik) 
	    sts_KH(3,3) = stsgp(3,ik) 
	    sts_KH(1,2) = stsgp(4,ik)
	    sts_KH(2,1) = stsgp(4,ik)
c
	    effsts      = effstsgp(ik)
	    effep       = effstegp(ik)
	    Temp        = thermo(ik)
c
	    call sts_evpAT2d(disp,vel,
     &               sts_KH,sts_PK1,Temp,
     &               effep,effsts,dt,
     &               theta,ik,xkc1,xkc2)
c
	    effstsgp(ik) = effsts
	    effstegp(ik) = effep
	    thermo(ik)   = Temp
c
	    stsgp(1,ik)  = sts_KH(1,1)
	    stsgp(2,ik)  = sts_KH(2,2)
	    stsgp(3,ik)  = sts_KH(3,3)
	    stsgp(4,ik)  = sts_KH(1,2)
c
         elseif (iMaterType .eq. 5) then
c
c........Thermo-elasto-viscoplastic material (Adiabatic) 
c
            do i=1,3
	       do j=1,3
                  sts_KH(i,j) = 0.0
               enddo
            enddo
c
            sts_KH(1,1) = stsgp(1,ik)
	    sts_KH(2,2) = stsgp(2,ik) 
	    sts_KH(3,3) = stsgp(3,ik) 
	    sts_KH(1,2) = stsgp(4,ik)
	    sts_KH(2,1) = stsgp(4,ik)
c
	    effsts      = effstsgp(ik)
	    effep       = effstegp(ik)
	    Temp        = thermo(ik)
            idamage     = Idamage_mark(ik)
c
	    call sts_evpATD2d(disp,vel,
     &              sts_KH,sts_PK1,Temp,
     &                  effep,effsts,dt,
     &            theta,ik,idamage,dvmk(ik),
     &                  test_min,ifail)
c
	    effstsgp(ik)     = effsts
	    effstegp(ik)     = effep
	    thermo(ik)       = Temp
            effRgp(ik)       = test_min(2)
	    ys_k(ik)         = test_min(3)
	    Idamage_mark(ik) = idamage
c
	    if(idamage .eq. 2) then
	       ik_eID = gpelem(ik)
	       do i = 1, nintElem
		  ike = lgops(i,ik_eID)
		  Idamage_mark(ike) = 2
	       enddo
	    endif
c
	    stsgp(1,ik)  = sts_KH(1,1)
	    stsgp(2,ik)  = sts_KH(2,2)
	    stsgp(3,ik)  = sts_KH(3,3)
	    stsgp(4,ik)  = sts_KH(1,2)
c
            if (Temp .gt. Temp_max) then
		Temp_max = Temp
	    endif
c
c-----------------------------------------------------
c
            if (test_m1 .lt. test_min(1)) then
		test_m1 = test_min(1)
	    endif
c
            if (test_m2 .lt. test_min(2)) then
		test_m2 = test_min(2)
	    endif
c
            if (test_m3 .lt. test_min(3)) then
		test_m3 = test_min(3)
	    endif
c         
	 else
c         
            write(*,*) 'wrong value of iMaterType=',iMaterType            
            stop
         endif
c         
c         ! store the stress values of Gauss Point
c         !*** internal forces
c
         mLoop = Lgp(ik)
         do iz = 1, mLoop
            ip = Lmgp(iz,ik)
c            
c            		
c			! by XGXGXGXGXGXG
c
c	     fint(1,ip)=fint(1,ip)
c     &                 +dvmk(ik)*(stress(1,1)*shpkdx(iz,ik)
c     &                           +stress(1,2)*shpkdy(iz,ik) )
c
c	     fint(2,ip)=fint(2,ip)
c     &                 +dvmk(ik)*(stress(2,1)*shpkdx(iz,ik)
c     &                           +stress(2,2)*shpkdy(iz,ik) )
c
c
c     ! note: sts_PK1 is store in normal order, not transpose
c
	     fint(1,ip)=fint(1,ip)
     &                 +dvmk(ik)*(sts_PK1(1,1)*shpkdx(iz,ik)
     &                           +sts_PK1(2,1)*shpkdy(iz,ik) )
c     
	     fint(2,ip)=fint(2,ip)
     &                 +dvmk(ik)*(sts_PK1(1,2)*shpkdx(iz,ik)
     &                           +sts_PK1(2,2)*shpkdy(iz,ik) )
         enddo
      enddo	! enddo ik
c
c.......... call external force ......
c
      call trabc(LmTract,vTract,xm,dxm,dvm,fext,nTract,
     &           numnp,
     &           istep,nstep,iextf_LoadingType,ifwLoad)
c
c.........Calulate force-deflection curve.............
c
       TF = 0.0d0
       BW = 0.0d0
       do ip = 1, nTractf
        ipt  = LmTractF(ip)
        TF   = TF + fint(2,ipt)*bweight(ip)
        BW   = BW + bweight(ip)
       enddo
c        TF   = TF/BW
c
        v_istep = real(istep) + 0.5d0
	call esspar(v_istep,nstep,dt,
     &              ivel_LoadingType,
     &              dispar,velpar,accpar)
c
	DFL = dispar * time
c
c
c----------------------------------------------------------
c     Acceleration Calculation and Corrector Phase
c
      if (iIntMeth .eq. 1) then
            ! predictor-corrector
	    ! Acceleration Calculation and Corrector Phase
c                        
         do ip = 1, numnp
            do idim = 1, 2
               acc(idim,ip)  = amass(ip)*(fext(idim,ip)-fint(idim,ip))
               vel(idim,ip)  = vel(idim,ip)+gamma*dt*acc(idim,ip)
               disp(idim,ip) = disp(idim,ip)+beta*dt*dt*acc(idim,ip)
            enddo
         enddo  ! ip            
c
c           ! Essential Boundary Condition : Consistent BC Method
c
         v_istep = real(istep) + 1          
         call xgessbc(LmDispbcX,LmDispbcY,
     &                vDispbcX,vDispbcY,
     &                LmVelbcX,LmVelbcY,
     &                vVelbcX,vVelbcY,
     &                disp,vel,acc,shpn,
     &                istep,nstep,v_istep,ivel_LoadingType,     
     &                ifwVelLoad)
c
      elseif (iIntMeth.eq.2) then
		! central difference
         do ip = 1, numnp
            do idim = 1, 2
               acc(idim,ip) = amass(ip)*(fext(idim,ip)-fint(idim,ip))
	       vel(idim,ip) = vel_last(idim,ip)
     &                      + 0.5*dt*(acc_last(idim,ip) + acc(idim,ip))
               disp(idim,ip)= disp_last(idim,ip)
     &                      + dt*vel_last(idim,ip)
     &                      + 0.5*dt*dt*acc_last(idim,ip)
            enddo
         enddo
         v_istep = real(istep) + 1          
         call xgessbc(LmDispbcX,LmDispbcY,
     &                vDispbcX,vDispbcY,
     &                LmVelbcX,LmVelbcY,
     &                vVelbcX,vVelbcY,
     &                disp,vel,acc,shpn,
     &                istep,nstep,v_istep,ivel_LoadingType,     
     &                ifwVelLoad)
c
      elseif (iIntMeth .eq. 3) then
c
c.............! Improved Euler-theta method !..............
c
         do ip = 1, numnp
            do idim = 1, 2
               acc(idim,ip) = amass(ip)*(fext(idim,ip)-fint(idim,ip))
c
	       vel(idim,ip) = vel_last(idim,ip) 
     &                      + dt*( (1.0-theta)*acc_last(idim,ip) 
     &                      + theta*acc(idim,ip) )   
c
            enddo ! idim
         enddo    ! ip
c
c
c        ! Essential Boundary Condition : Consistent BC Method
c
         v_istep = real(istep) + 1.0d0
c
         call bcvel(LmDispbcX,LmDispbcY,
     &              vDispbcX,vDispbcY,
     &              LmVelbcX,LmVelbcY,
     &              vVelbcX,vVelbcY,
     &              vel,shpn,
     &              v_istep,nstep,ivel_LoadingType,
     &              ifwVelLoad)
c
c         do ip = 1, numnp
c	    do idim = 1,2
c	       disp(idim,ip) = disp_last(idim,ip) 
c     &                      + dt*( (1.0-theta)*vel_last(idim,ip) 
c     &                      + theta*vel(idim,ip) )   
c            enddo
c         enddo ! ip
c
c
c         call bcdisp(LmDispbcX,LmDispbcY,
c     &               vDispbcX,vDispbcY,
c     &               LmVelbcX,LmVelbcY,
c     &               vVelbcX,vVelbcY,
c     &               disp,vel,shpn,
c     &               v_istep,nstep,ivel_LoadingType)
cc
c
      else
         write(*,*) 'wrong number of iIntMeth=',iIntMeth
         stop         
      endif  !iIntMeth
c
c      	 ! store this time step values
      do ip=1,numnp
         do idim=1,2
            disp_last(idim,ip) = disp(idim,ip)
            vel_last(idim,ip)  = vel(idim,ip)
            acc_last(idim,ip)  = acc(idim,ip)
         enddo
      enddo
c      
c
c-----------------------------------------------------
c**     Update : xs,xsk
c
      do ipt = 1, numnp
         dixxx = 0.0d0
         diyyy = 0.0d0
         jLoop = Lnp(ipt)
         do jp = 1, jLoop
            jpt = Lmnp(jp,ipt)
            dixxx = dixxx + shpn(jp,ipt)*disp(1,jpt)
            diyyy = diyyy + shpn(jp,ipt)*disp(2,jpt)
         enddo
         xs(1,ipt) = xm(1,ipt) + dixxx
         xs(2,ipt) = xm(2,ipt) + diyyy
      enddo
c
      do ik = 1, mgk
         dixxx = 0.0d0
         diyyy = 0.0d0
         jLoop = Lgp(ik)
         do jp = 1, jLoop
            jpt = Lmgp(jp,ik)
            dixxx = dixxx + shpk(jp,ik)*disp(1,jpt)
            diyyy = diyyy + shpk(jp,ik)*disp(2,jpt)
         enddo
         xsk(1,ik) = xmk(1,ik) + dixxx
         xsk(2,ik) = xmk(2,ik) + diyyy
      enddo
c      
c..........! update the nodal disp, vel and acc
c
      do ipt = 1, numnp
         dixxx = 0.0d0
         diyyy = 0.0d0
         jLoop = Lnp(ipt)
         do jp = 1, jLoop
            jpt = Lmnp(jp,ipt)
            dixxx = dixxx + shpn(jp,ipt)*vel(1,jpt)
            diyyy = diyyy + shpn(jp,ipt)*vel(2,jpt)
         enddo
         disp_node(1,ipt) = dixxx
         disp_node(2,ipt) = diyyy
      enddo
c                
      do ipt = 1, numnp
         dixxx = 0.0d0
         diyyy = 0.0d0
         jLoop = Lnp(ipt)
         do jp = 1, jLoop
            jpt = Lmnp(jp,ipt)
            dixxx = dixxx + shpn(jp,ipt)*vel(1,jpt)
            diyyy = diyyy + shpn(jp,ipt)*vel(2,jpt)
         enddo
         vel_node(1,ipt) = dixxx
         vel_node(2,ipt) = diyyy
      enddo
c
      do ipt = 1, numnp
         dixxx = 0.0d0
         diyyy = 0.0d0
         jLoop = Lnp(ipt)
         do jp = 1, jLoop
            jpt = Lmnp(jp,ipt)
            dixxx = dixxx + shpn(jp,ipt)*acc(1,jpt)
            diyyy = diyyy + shpn(jp,ipt)*acc(2,jpt)
         enddo
         acc_node(1,ipt) = dixxx
         acc_node(2,ipt) = diyyy
      enddo
c
c......Interpolate effective variables
c      at the nodal points
c
c-----------------------------------------------------
c     output
c-----------------------------------------------------
c
             !***  output the result of point # iplot
c
      if (mod(istep,nPlotFrq).eq. 0) then
c      
         if ( (iplot.ne.0).and.(iplot_gp.ne.0) ) then 
            write(*,'(1x,i8,6e11.3)' ) 
     &       istep,disp(1,iplot),disp(2,iplot),
     &       xs(1,iplot),xs(2,iplot),
     &       effstsgp(iplot_gp),effstegp(iplot_gp)
         else
            if (iplot.ne.0) then
               write(*,'(1x,i8,6e11.3)' ) 
     &         istep,disp(1,iplot),disp(2,iplot),
     &         xs(1,iplot),xs(2,iplot)
            elseif ( iplot_gp.ne.0) then
               write(*,'(1x,i8,6e11.3)' ) 
     &         istep,effstsgp(iplot_gp),
     &         effstegp(iplot_gp)
            endif
         endif
c      
            write(ifwDisPT,'(1x,i8,3e14.7)')
     &         istep,time,temp_max
            write(ifwVelPT,'(1x,i8,3e13.5)')
     &         istep,time,vel(1,iplot),vel(2,iplot)
            write(ifwStsPT,'(1x,i8,4e13.5)')
     &         istep,time,test_m1,test_m2,test_m3
	    write(ifwTimeStep,3384)istep,time,temp_max
c         
         call flush(ifwDisPT)
         call flush(ifwVelPT)
         call flush(ifwStsPT)
	 call flush(ifwTimeStep)
                                
      endif ! nPlotFrq
c      
c-------------------------------------------------
c   (2)	 output the intermid result
c-------------------------------------------------
c
      if((mod(istep,nOutputFrq) .eq. 0)
     &   .or. (istep .eq. nstep)) then
c-----------------------------------------------
c.......Reproducing the stress at nodal points
c-----------------------------------------------
c
c......Interpolate effective variables
c      at the nodal points
c
       if (iMaterType .ne. 0) then
       do ipt = 1, numnp
	  do j = 1,4
	     stsnpt(j,ipt) = 0.0d0
	  enddo
	  effEnpt(ipt) = 0.0d0
	  effSnpt(ipt) = 0.0d0
	  effEopt(ipt) = 0.0d0
	  effTemp(ipt) = 0.0d0
	  effRate(ipt) = 0.0d0
	  effMax(ipt)  = 0.0d0
       enddo
c 
       do ik = 1, mgk
	  jLoop = Lgp(ik)
	  do jp = 1, jLoop
	     jpt = Lmgp(jp,ik)
c
	     effEopt(jpt) = effEopt(jpt)
     &                    + shpk(jp,ik)
          enddo
       enddo
c
       do ik = 1, mgk
	  jLoop = Lgp(ik)
	  do jp = 1, jLoop
	     jpt = Lmgp(jp,ik)
c
	      do k = 1, 4
		 stsnpt(k,jpt) = stsnpt(k,jpt)
     &          + stsgp(k,ik)*shpk(jp,ik)/effEopt(jpt)
              enddo
c
	      effEnpt(jpt) = effEnpt(jpt)
     &       + effstegp(ik)*shpk(jp,ik)/effEopt(jpt)
c
	      effSnpt(jpt) = effSnpt(jpt)
     &       + effstsgp(ik)*shpk(jp,ik)/effEopt(jpt)
c
	      effTemp(jpt) = effTemp(jpt)
     &       + thermo(ik)*shpk(jp,ik)/effEopt(jpt)
c
c
	     effRate(jpt) = effRate(jpt)
     &       + effRgp(ik)*shpk(jp,ik)/effEopt(jpt)
c
             effMax(jpt) = effMax(jpt)
     &       + ys_k(ik)*shpk(jp,ik)/effEopt(jpt)
c
          enddo ! jp
       enddo    ! ik
c 
	call  write2d(xm,xs,disp,vel,acc,
     &               disp_node,vel_node,acc_node,
     &               lnods,nelem,numnp,mgk,
     &               istep,iOutput,fhead)
c
        else ! iMaterType = 0 (special)
	    ifwout=31
c
            iOutput=iOutput+1
	    write(*,*) 'iOutput=',iOutput
	    write(*,3384)istep,dt,time
c
cc##
c
	    ist = Istage + 1
	    call Int2Str(iOutput,tmpstr1)
	    single = suffix(ist:ist)
	    call appstr(single,tmpstr1,tmpstr2)
	    call appstr('.gpsts',tmpstr2,tmpstr3)
	    call appstr(fhead,tmpstr3,fwname)
	    write(*,*) 'fwname=',fwname
	    open(ifwout,file=fwname)
c
	    write(ifwout,'(1x,a)')' VARIABLES = "X","Y",
     &"SS","EE","TT","S1","S12"'
c
	    do ik = 1, mgk
	       write(ifwout,'(1x,7e16.7)') xsk(1,ik),xsk(2,ik),
     &                      effstsgp(ik),effstegp(ik),thermo(ik),
     &                      stsgp(1,ik),stsgp(4,ik)
            enddo
	    close(ifwout)
	endif
      endif ! mod
      
 100  continue
c
c***********************************************************
c     time loop ends
c***********************************************************
c
c## H-Refinement begin:
c
       Istage = Istage + 1
c
c  Step (1): Weighted average algorithm
c
       do ipt = 1, numnp
          effEnpt(ipt)  = 0.0
	  effEopt(ipt)  = 0.0
	  Leff_adp(ipt) = 0
       enddo
c
       do ik = 1, mgk
	  jLoop = Lgp(ik)
	  do jp = 1, jLoop
	    jpt = Lmgp(jp,ik)
c
	    effEopt(jpt) = effEopt(jpt)
     &                   + shpk(jp,ik)
	  enddo ! jp
       enddo    ! ik
c
       do ik = 1, mgk
	  jLoop = Lgp(ik)
	  do jp = 1, jLoop
	     jpt= Lmgp(jp,ik)
c
	     effEnpt(jpt) = effEnpt(jpt)
     &       + effstegp(ik)*shpk(jp,ik)/effEopt(jpt)
	  enddo  ! jp
       enddo     ! ik
c
c
        if(iAdapt .eq. 1 .and. Istage .le. iad_stage) then
c
           print *, 'entery O. K.'
c
           call hadapt(xm,dxm,dvm,xmk,dvmk,
     &        LmTract,vTract,
     &        LmDispbcX,LmDispbcY,vDispbcX,vDispbcY,
     &        LmVelbcX,LmVelbcY,vVelbcX,vVelbcY,
     &        igp_elemind,igp_LocintInd,gp_loc2D, ! for FEM
     &        nelem,numnp,mgk,
     &        effEnpt,ifwAdapt,fhead)
c
           print *, 'exit O. K.'
c#
	   go to 38
        elseif(iAdapt .eq. 2 .and. Istage .le. iad_stage) then
c
c..........p (or wavelet adaptivity )
c
	   numnp      = numnp + 2 * numadp
	   dt         = 0.5 * dt
	   nstep      = 2 * nstep
	   nOutputFrq = 2 * nOutputFrq
c#
           go to 38
        endif ! iAdapt
c##
c
c##
c.....output the final deformation shape
c
c      write(ifwTimeStep,3384)istep,dt,time
c      
 97   format(2(1x,e12.5))
 98   format(3(1x,e12.5))
 3384 format(i7,2(1x,e12.5))
c
      close(10)
      close(ifwTimeStep)
      close(ifwDisPT)
      close(ifwVelPT)
      close(ifwStsPT)
      close(ifwLoad)
      close(ifwVelLoad)
c
      end
c
c
