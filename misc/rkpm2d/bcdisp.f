      subroutine bcdisp(LmDispbcX,LmDispbcY,
     &                  vDispbcX,vDispbcY,
     &                  LmVelbcX,LmVelbcY,
     &                  vVelbcX,vVelbcY,
     &                  disp,vel,shpn,
     &                  v_istep,nstep,ivel_LoadingType)
c
c--------------------------------------------------------------
c
c
c   Note:
c
c   In DYME2D (Dynamic Meshless (FEM/RKPM)) Solid 2-D code:
c
c   Three types essential b.c. are specified:
c
c   (1) the purely displacement b.c., by which, we mean:
c
c       both the velocity and acceleration on the boundary
c       are zero.
c
c   (2) the purely velocity b.c. , by which, we mean:
c
c       the acceleration on the b.c. is zero.
c       and, the displacement on the b.c. is prescribed
c       in exact fashion
c
c   (3) the general velocity b.c., 
c
c         the acceleration on the b.c. is prescribed;
c         the velocity on the b.c. is prescribed;
c    and, the displacement on the b.c. is prescribed;
c       
c   To simplify the input process, only velocity b.c card
c   and displacement b.c. card are read.
c   
c   Thus, in bcdisp.f, the displacement b.c. on both Vbc
c   and Dbc are enforced.
c
c
c
c--------------------------------------------------------------
c
      implicit none
      include 'parameter.h'
c      
      integer nstep,ivel_LoadingType
      integer LmDispbcX(maxDispbcX),LmDispbcY(maxDispbcY),
     &        LmVelbcX(maxVelbcX),LmVelbcY(maxVelbcY)
      real*8  vDispbcX(maxDispbcX),vDispbcY(maxDispbcY),
     &        vVelbcX(maxVelbcX),vVelbcY(maxVelbcY)
c      
      real*8  disp(2,maxNumnp)
      real*8  vel(2,maxNumnp)
      real*8  shpn(maxNumnp,mnsch)
      real*8  shBCX(maxEssbcX,maxEssbcX),shBCY(maxEssbcY,maxEssbcY)
      real*8  time,dt,v_istep
c      
c	! local vars
c
      integer LmEssbcX(maxEssbcX),LmEssbcY(maxEssbcY)
      real*8  vEssbcX_disp(maxEssbcX),vEssbcY_disp(maxEssbcY)  
c           
      real*8 ftemp,shbc
      real*8 dispar,velpar,accpar
c      
      integer lnods(maxNode,maxElem),
     &        Lnp(maxNumnp),Lmnp(mnsch,maxNumnp),
     &        Lgp(maxGP),Lmgp(mnsch,maxGP)
c
      integer nEssbcX,nEssbcY,
     &        iEssbcX,jEssbcX,iEssbcY,jEssbcY,
     &        iDispbcX,iDispbcY,
     &        iVelbcX,iVelbcY,
     &        ip,jp,iloop,mloop
c      
      integer imethRKPM,imethFEM,imeth,iIntMeth,iInter,
     &        nTract,nDispbcX,nDispbcY,
     &        nVelbcX,nVelbcY
c
      common /shapeBC/ shBCX,shBCY
      common /connect/Lnp,Lmnp,Lgp,Lmgp,lnods
      common /stepT/ dt,time
      common /meth/imethRKPM,imethFEM,imeth,iIntMeth,iInter
      common /bound/ nTract,nDispbcX,nDispbcY,nVelbcX,nVelbcY
c
c
c Essential Boundary Condition : Consistent BC Method ( for RKPM )
c
c
      call esspar(v_istep,nstep,dt,ivel_LoadingType,
     &            dispar,velpar,accpar)
c
      if ( imeth.eq.imethFEM ) then
         do iDispbcX = 1,nDispbcX
            ip         = LmDispbcX(iDispbcX)
            disp(1,ip) = vDispbcX(iDispbcX)
         enddo
c
         do iDispbcY = 1,nDispbcY
            ip         = LmDispbcY(iDispbcY)
            disp(2,ip) = vDispbcY(iDispbcY)
         enddo
c
         do iVelbcX = 1,nVelbcX
            ip         = LmVelbcX(iVelbcX)
            disp(1,ip) = vVelbcX(iVelbcX)*dispar
         enddo
c
         do iVelbcY = 1,nVelbcY
            ip         = LmVelbcY(iVelbcY)
            disp(2,ip) = vVelbcY(iVelbcY)*dispar
         enddo
c
         return
      endif  ! if imethFEM
c
c
c---------- below is for RKPM
c 
      iEssbcX = 0
      iEssbcY = 0
c      
      do iDispbcX = 1,nDispbcX
         iEssbcX               = iEssbcX + 1
         LmEssbcX(iEssbcX)     = LmDispbcX(iDispbcX)
         vEssbcX_disp(iEssbcX) = vDispbcX(iDispbcX)
      enddo
c
      do iDispbcY = 1,nDispbcY
         iEssbcY               = iEssbcY + 1
         LmEssbcY(iEssbcY)     = LmDispbcY(iDispbcY)
         vEssbcY_disp(iEssbcY) = vDispbcY(iDispbcY)
      enddo
c      
      do iVelbcX = 1,nVelbcX
         iEssbcX               = iEssbcX + 1
         LmEssbcX(iEssbcX)     = LmVelbcX(iVelbcX)
         vEssbcX_disp(iEssbcX) = vVelbcX(iVelbcX)*dispar
      enddo
c
      do iVelbcY = 1,nVelbcY
         iEssbcY               = iEssbcY + 1
         LmEssbcY(iEssbcY)     = LmVelbcY(iVelbcY)
         vEssbcY_disp(iEssbcY) = vVelbcY(iVelbcY)*dispar
      enddo
c
      nEssbcX = nDispbcX + nVelbcX
      nEssbcY = nDispbcY + nVelbcY
c
c----- (1)  X-Direction ----------
c
      if ( nEssbcX .gt. 0 ) then
c
         do iEssbcX=1, nEssbcX
           ip         = LmEssbcX(iEssbcX)
           disp(1,ip) = 0.00
         enddo

         do iEssbcX=1,nEssbcX
           ip = LmEssbcX(iEssbcX)
           ftemp = 0.0
           mLoop = Lnp(ip)
           do iloop = 1, mLoop
             jp = Lmnp(iloop,ip)
             shbc = shpn(iloop,ip)
             ftemp = ftemp + shbc*disp(1,jp)
           enddo
           vEssbcX_disp(iEssbcX) = vEssbcX_disp(iEssbcX) - ftemp
         enddo
c
         do iEssbcX = 1, nEssbcX
            ip         = LmEssbcX(iEssbcX)
            do jEssbcX = 1, nEssbcX
               disp(1,ip) = disp(1,ip)
     &       + shBCX(iEssbcX,jEssbcX)*vEssbcX_disp(jEssbcX)
            enddo
         enddo
c
      endif    ! nEssbcX > 0            
c      
c----- (2) Y- Direction ---------
c

      if ( nEssbcY .gt. 0 ) then
c
         do iEssbcY=1,nEssbcY
           ip = LmEssbcY(iEssbcY)
           disp(2,ip) = 0.00d0
         enddo
c
         do iEssbcY=1,nEssbcY
           ip = LmEssbcY(iEssbcY)
           ftemp = 0.0d0
           mLoop = Lnp(ip)
           do iloop = 1, mLoop
             jp   = Lmnp(iloop,ip)
             shbc = shpn(iloop,ip)
             ftemp = ftemp + shbc*disp(2,jp)
           enddo
           vEssbcY_disp(iEssbcY) = vEssbcY_disp(iEssbcY) - ftemp
         enddo
c
         do iEssbcY = 1, nEssbcY
            ip         = LmEssbcY(iEssbcY)
            do jEssbcY = 1, nEssbcY
	       disp(2,ip) = disp(2,ip) 
     &       + shBCY(iEssbcY,jEssbcY)*vEssbcY_disp(jEssbcY)
            enddo
         enddo
c
      endif    ! nEssbcY >0
c
      return
      end
c
