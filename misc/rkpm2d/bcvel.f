      subroutine bcvel(LmDispbcX,LmDispbcY,
     &                 vDispbcX,vDispbcY,
     &                 LmVelbcX,LmVelbcY,
     &                 vVelbcX,vVelbcY,
     &                 vel,shpn,
     &                 v_istep,nstep,ivel_LoadingType,
     &                 ifwVelLoad)
c
c--------------------------------------------------------
c
      implicit none
      include 'parameter.h'
c      
      real*8  v_istep
      integer nstep,ivel_LoadingType,ifwVelLoad
      integer LmVelbcX(maxVelbcX),LmVelbcY(maxVelbcY),
     &        LmDispbcX(maxDispbcX),LmDispbcY(maxDispbcY)
      real*8  vVelbcX(maxVelbcX),vVelbcY(maxVelbcY), 
     &        vDispbcX(maxDispbcX),vDispbcY(maxDispbcY) 
c      
      real*8  vel(2,maxNumnp),shpn(mnsch,maxNumnp)
      real*8  shBCX(maxEssbcX,maxEssbcX),shBCY(maxEssbcY,maxEssbcY)
      real*8  time,dt
c      
	! local vars
      integer LmEssbcX(maxEssbcX),LmEssbcY(maxEssbcY)
      real*8  vEssbcX(maxEssbcX),vEssbcY(maxEssbcY)  
c           
      real*8 ftemp,shbc,dispar,velpar,accpar
c      
      integer nEssbcX,nEssbcY,
     &        iEssbcX,jEssbcX,iEssbcY,jEssbcY,
     &        iVelbcX,iVelbcY,
     &        iDispbcX,iDispbcY
      integer ip,jp,iloop,mloop
c      
      integer lnods(maxNode,maxElem),
     &        Lnp(maxNumnp),Lmnp(mnsch,maxNumnp),
     &        Lgp(maxGP),Lmgp(mnsch,maxGP)
c
      integer imethRKPM,imethFEM,imeth,iIntMeth,iInter,
     &        nTract,nDispbcX,nDispbcY,
     &        nVelbcX,nVelbcY
c
      common /shapeBC/shBCX,shBCY
      common /connect/Lnp,Lmnp,Lgp,Lmgp,lnods
      common /stepT/dt,time
      common /meth/imethRKPM,imethFEM,imeth,iIntMeth,iInter
      common /bound/nTract,nDispbcX,nDispbcY,nVelbcX,nVelbcY
c      
c
c     Essential Boundary Condition : Consistent BC Method ( for RKPM )
c
c
      call esspar(v_istep,nstep,dt,ivel_LoadingType,
     &            dispar,velpar,accpar)
      
      if ( 1 .eq. 0) then
         write(ifwVelLoad,'(1x,f11.3,i8,2e15.6)') 
     &   v_istep,nstep,real(v_istep)/real(nstep),velpar
      endif
c
      if ( imeth.eq.imethFEM ) then
c      
         do iDispbcX = 1,nDispbcX
            ip        = LmDispbcX(iDispbcX)
            vel(1,ip) = 0.0
         enddo
c
         do iDispbcY=1,nDispbcY
            ip        = LmDispbcY(iDispbcY)
            vel(2,ip) = 0.0
         enddo
c
         do iVelbcX = 1,nVelbcX
            ip        = LmVelbcX(iVelbcX)
            vel(1,ip) = vVelbcX(iVelbcX)*velpar
         enddo
c
         do iVelbcY = 1,nVelbcY
            ip=LmVelbcY(iVelbcY)
            vel(2,ip)  = vVelbcY(iVelbcY)*velpar
         enddo
c         
         return
      endif  ! if imethFEM
c
c
c------------- below if for RKPM --------------	
c
c
      iEssbcX = 0
      iEssbcY = 0
c      
      do iDispbcX=1,nDispbcX
         iEssbcX           = iEssbcX+1
         LmEssbcX(iEssbcX) = LmDispbcX(iDispbcX)
         vEssbcX(iEssbcX)  = 0.0d0
      enddo
c
      do iDispbcY=1,nDispbcY
         iEssbcY           = iEssbcY+1
         LmEssbcY(iEssbcY) = LmDispbcY(iDispbcY)
         vEssbcY(iEssbcY)  = 0.0d0
      enddo
c
      do iVelbcX=1,nVelbcX
         iEssbcX           = iEssbcX+1
         LmEssbcX(iEssbcX) = LmVelbcX(iVelbcX)
         vEssbcX(iEssbcX)  = vVelbcX(iVelbcX)*velpar
      enddo
c
      do iVelbcY=1,nVelbcY
         iEssbcY           = iEssbcY+1
         LmEssbcY(iEssbcY) = LmVelbcY(iVelbcY)
         vEssbcY(iEssbcY)  = vVelbcY(iVelbcY)*velpar
      enddo
c
      nEssbcX = nVelbcX + nDispbcX 
      nEssbcY = nVelbcY + nDispbcY
c
c----- (1) X Direction ----------
c  
      if ( nEssbcX .gt. 0 ) then
c
         do iEssbcX=1,nEssbcX
           ip = LmEssbcX(iEssbcX)
           vel(1,ip) = 0.0d0
         enddo
c
         do iEssbcX=1,nEssbcX
           ip    = LmEssbcX(iEssbcX)
           ftemp = 0.0d0
           mLoop = Lnp(ip)
           do iloop = 1, mLoop
             jp    = Lmnp(iloop,ip)
             shbc  = shpn(iloop,ip)
             ftemp = ftemp + shbc*vel(1,jp)
           enddo  ! iloop
           vEssbcX(iEssbcX) = vEssbcX(iEssbcX) - ftemp
         enddo
c
         do iEssbcX = 1, nEssbcX
            ip         = LmEssbcX(iEssbcX)
            do jEssbcX = 1, nEssbcX
               vel(1,ip) = vel(1,ip)
     &       + shBCX(iEssbcX,jEssbcX)*vEssbcX(jEssbcX)
	    enddo
	 enddo
c
      endif 	! if nEssbcX > 0
c      
c----- (2) Y-Direction ---------
c
      if ( nEssbcY .gt. 0 ) then
c
         do iEssbcY = 1, nEssbcY
           ip        = LmEssbcY(iEssbcY)
           vel(2,ip) = 0.0d0
         enddo
c
         do iEssbcY=1,nEssbcY
           ip = LmEssbcY(iEssbcY)
           ftemp = 0.0d0
           mLoop = Lnp(ip)
           do iloop = 1, mLoop
              jp   = Lmnp(iloop,ip)
              shbc = shpn(iloop,ip)
              ftemp = ftemp + shbc*vel(2,jp)
           enddo  ! iloop
           vEssbcY(iEssbcY) = vEssbcY(iEssbcY) - ftemp
         enddo
c
         do iEssbcY = 1, nEssbcY
            ip = LmEssbcY(iEssbcY)
            do jEssbcY = 1, nEssbcY
               vel(2,ip) = vel(2,ip)
     &       + shBCY(iEssbcY,jEssbcY)*vEssbcY(jEssbcY)
	    enddo
	 enddo
c
      endif	! if nEssbcY > 0
c      
      return
      end
c
