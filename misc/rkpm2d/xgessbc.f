      subroutine xgessbc(LmDispbcX,LmDispbcY,
     &                   vDispbcX,vDispbcY,
     &                   LmVelbcX,LmVelbcY,
     &                   vVelbcX,vVelbcY,
     &                   disp,vel,acc,shpn,
     &                   istep,nstep,v_istep,ivel_LoadingType,
     &                   ifwVelLoad)
c
c------------------------------------------------------
c
c*** input
c
c------------------------------------------------------
c
      implicit none
      include 'parameter.h'
c      
      integer istep,nstep,ivel_LoadingType
      integer ifwVelLoad
      integer LmDispbcX(maxDispbcX),LmDispbcY(maxDispbcY),
     &        LmVelbcX(maxVelbcX),LmVelbcY(maxVelbcY)
      real*8  vDispbcX(maxDispbcX),vDispbcY(maxDispbcY),
     &        vVelbcX(maxVelbcX),vVelbcY(maxVelbcY) 
      real*8  shBCX(maxEssbcX,maxEssbcX),shBCY(maxEssbcY,maxEssbcY)
c      
      real*8  disp(2,maxNumnp),vel(2,maxNumnp),acc(2,maxNumnp)
      real*8  shpn(mnsch,maxNumnp)
      real*8  time,dt
c
c	! local vars
c
      integer LmEssbcX(maxEssbcX),LmEssbcY(maxEssbcY)
      real*8  vEssbcX_disp(maxEssbcX),vEssbcY_disp(maxEssbcY),  
     &        vEssbcX_vel(maxEssbcX),vEssbcY_vel(maxEssbcY),  
     &        vEssbcX_acc(maxEssbcX),vEssbcY_acc(maxEssbcY)  
c
      real*8  shbc
      real*8  fbc_disp(maxEssbc),
     &        fbc_vel(maxEssbc),
     &        fbc_acc(maxEssbc)
c           
      real*8 ftemp1,ftemp2,ftemp3      
      real*8 dispar,velpar,accpar
c      
      integer nEssbcX,nEssbcY,
     &        iEssbcX,jEssbcX,iEssbcY,jEssbcY,
     &        iDispbcX,iDispbcY,iVelbcX,iVelbcY
      integer ip,jp,iloop,mloop
c      
      integer imethRKPM,imethFEM,imeth,iIntMeth,iInter,
     &        nTract,nDispbcX,nDispbcY,nVelbcX,nVelbcY
c
      integer lnods(maxNode,maxElem),
     &        Lnp(maxNumnp),Lmnp(mnsch,maxNumnp),
     &        Lgp(maxGP),Lmgp(mnsch,maxGP)
c
      common /shapeBC/shBCX,shBCY
      common /connect/Lnp,Lmnp,Lgp,Lmgp,lnods
      common /stepT/dt,time
      common /meth/imethRKPM,imethFEM,imeth,iIntMeth,iInter
      common /bound/ nTract,nDispbcX,nDispbcY,nVelbcX,nVelbcY
c      
      real*8 v_istep
c      
c
c     Essential Boundary Condition : Consistent BC Method ( for RKPM )
c
c
      call esspar(v_istep,nstep,dt,ivel_LoadingType,
     &            dispar,velpar,accpar)
c
      if (1 .eq. 0) then
         write(ifwVelLoad,'(1x,2i8,2e15.6)') istep,nstep,
     &   real(istep)/real(nstep),velpar
      endif
c
      if ( imeth.eq.imethFEM ) then
c
         do iDispbcX=1,nDispbcX
            ip         = LmDispbcX(iDispbcX)
            disp(1,ip) = vDispbcX(iDispbcX)
	    vel(1,ip)  = 0.0            
            acc(1,ip)  = 0.0
         enddo
c
         do iDispbcY=1,nDispbcY
            ip         = LmDispbcY(iDispbcY)
            disp(2,ip) = vDispbcY(iDispbcY)
            vel(2,ip)  = 0.0            
            acc(2,ip)  = 0.0
         enddo

         do iVelbcX=1,nVelbcX
            ip         = LmVelbcX(iVelbcX)
            disp(1,ip) = vVelbcX(iVelbcX)*dispar
            vel(1,ip)  = vVelbcX(iVelbcX)*velpar
            acc(1,ip)  = vVelbcX(iVelbcX)*accpar
         enddo

         do iVelbcY=1,nVelbcY
            ip=LmVelbcY(iVelbcY)
            disp(2,ip) = vVelbcY(iVelbcY)*dispar
            vel(2,ip)  = vVelbcY(iVelbcY)*velpar
            acc(2,ip)  = vVelbcY(iVelbcY)*accpar
         enddo
         
         return
      endif  ! if imethFEM
c
c
c.........Below is for RKPM..............
c
c
      iEssbcX = 0
      iEssbcY = 0
c      
      do iDispbcX = 1,nDispbcX
         iEssbcX               = iEssbcX+1
         LmEssbcX(iEssbcX)     = LmDispbcX(iDispbcX)
         vEssbcX_disp(iEssbcX) = vDispbcX(iDispbcX)
         vEssbcX_vel(iEssbcX)  = 0.0
         vEssbcX_acc(iEssbcX)  = 0.0
      enddo
c
      do iDispbcY = 1,nDispbcY
         iEssbcY               = iEssbcY+1
         LmEssbcY(iEssbcY)     = LmDispbcY(iDispbcY)
         vEssbcY_disp(iEssbcY) = vDispbcY(iDispbcY)
         vEssbcY_vel(iEssbcY)  = 0.0
         vEssbcY_acc(iEssbcY)  = 0.0
      enddo
c      
      do iVelbcX=1,nVelbcX
         iEssbcX               = iEssbcX+1
         LmEssbcX(iEssbcX)     = LmVelbcX(iVelbcX)
         vEssbcX_disp(iEssbcX) = vVelbcX(iVelbcX)*dispar
         vEssbcX_vel(iEssbcX)  = vVelbcX(iVelbcX)*velpar
         vEssbcX_acc(iEssbcX)  = vVelbcX(iVelbcX)*accpar
      enddo
c
      do iVelbcY=1,nVelbcY
         iEssbcY               = iEssbcY+1
         LmEssbcY(iEssbcY)     = LmVelbcY(iVelbcY)
         vEssbcY_disp(iEssbcY) = vVelbcY(iVelbcY)*dispar
         vEssbcY_vel(iEssbcY)  = vVelbcY(iVelbcY)*velpar
         vEssbcY_acc(iEssbcY)  = vVelbcY(iVelbcY)*accpar
      enddo
c
      nEssbcX=nDispbcX+nVelbcX
      nEssbcY=nDispbcY+nVelbcY
c
c
c------------- (1) X-direction ----------
c  
c.......Important Trick...........
c
      do iEssbcX=1,nEssbcX
        ip         = LmEssbcX(iEssbcX) 
        disp(1,ip) = 0.0
        vel(1,ip)  = 0.0
	acc(1,ip)  = 0.0
      enddo

      do iEssbcX=1,nEssbcX
        ip = LmEssbcX(iEssbcX)
        ftemp1 = 0.0
        ftemp2 = 0.0
        ftemp3 = 0.0
        mLoop = Lnp(ip)
        do iloop = 1, mLoop
          jp = Lmnp(iloop,ip)
          shbc = shpn(iloop,ip)
          ftemp1 = ftemp1 + shbc*disp(1,jp)
          ftemp2 = ftemp2 + shbc*vel(1,jp)
          ftemp3 = ftemp3 + shbc*acc(1,jp)
        enddo
        fbc_disp(iEssbcX) = vEssbcX_disp(iEssbcX) - ftemp1
        fbc_vel(iEssbcX)  = vEssbcX_vel(iEssbcX)  - ftemp2
        fbc_acc(iEssbcX)  = vEssbcX_acc(iEssbcX)  - ftemp3
      enddo
c
      do iEssbcX = 1, nEssbcX
         ip         = LmEssbcX(iEssbcX)
	 do jEssbcX = 1, nEssbcX
            disp(1,ip) = disp(1,ip) 
     &                 + shBCX(iEssbcX,jEssbcX)*fbc_disp(jEssbcX)
            vel(1,ip)  = vel(1,ip)
     &                 + shBCX(iEssbcX,jEssbcX)*fbc_vel(jEssbcX)
            acc(1,ip)  = acc(1,ip)
     &                 + shBCX(iEssbcX,jEssbcX)*fbc_acc(jEssbcX)
	 enddo
      enddo
c      
c      
c----- (2) Y-Direction ---------
c
      do iEssbcY=1,nEssbcY
        ip = LmEssbcY(iEssbcY) 
        disp(2,ip) = 0.00
        vel(2,ip)  = 0.00
        acc(2,ip)  = 0.00
      enddo

      do iEssbcY=1,nEssbcY
        ip = LmEssbcY(iEssbcY)
        ftemp1 = 0.00
        ftemp2 = 0.00
        ftemp3 = 0.00
        mLoop = Lnp(ip)
        do iloop = 1, mLoop
          jp = Lmnp(iloop,ip)
          shbc = shpn(iloop,ip)
          ftemp1 = ftemp1 + shbc*disp(2,jp)
          ftemp2 = ftemp2 + shbc*vel(2,jp)
          ftemp3 = ftemp3 + shbc*acc(2,jp)
        enddo
        fbc_disp(iEssbcY) = vEssbcY_disp(iEssbcY) - ftemp1
        fbc_vel(iEssbcY)  = vEssbcY_vel(iEssbcY)  - ftemp2
        fbc_acc(iEssbcY)  = vEssbcY_acc(iEssbcY)  - ftemp3
      enddo
c
      do iEssbcY = 1, nEssbcY
         ip = LmEssbcY(iEssbcY)
	 do jEssbcY = 1, nEssbcY
            disp(2,ip) = disp(2,ip)
     &                 + shBCY(iEssbcY,jEssbcY)*fbc_disp(jEssbcY)
            vel(2,ip)  = vel(2,ip)
     &                 + shBCY(iEssbcY,jEssbcY)*fbc_vel(jEssbcY)
            acc(2,ip)  = acc(2,ip)
     &                 + shBCY(iEssbcY,jEssbcY)*fbc_acc(jEssbcY)
	 enddo
      enddo
c
      return
      end
c
