      subroutine search2d(xmk,mgk,igp_elemind,
     &                    igp_LocIntInd,gp_loc2d)
c
c****************************************************************     
c 
c            
c*** below it will do :
c   For FEM or RKPM method:
c     (1) search nodes for each Gaussian point
c	  generate the lists: Lgp(ik) and Lmgp(ik,ip)
c                             where  ip=1..Lgp(ik)
c                                    ik=1..mgk
c
c     (2) search nodes for each node
c	  generate the lists: Lnp(ip) and Lmnp(ip,jp)
c                             where  jp=1..Lnp(ip)
c                                    ip=1..numnp
c
c     (2) cal the values of shape functions at each node and each
c         Gaussian point
c         generate the arrays: 
c             shpn(jp,ipt) -- shape function N (x ) at each node x
c                                            i  j                j
c                                where ip=1..Lnp(jp)
c                                      jp=1..numnp
c             shpk(jp,ik) -- shape function N (x ) at each gauss pt x
c                                            j  k                    k
c                                where ip=1..Lgp(ik)
c                                      ik=1..mgk
c
c
c-------------------------------------------------------------
c
      implicit none
      include 'parameter.h'
c
      real*8 xm(2,maxNumnp),dxm(2,maxNumnp),dvm(maxNumnp)
      real*8 xmk(2,maxGP)
      real*8 gp_loc2D(2,maxIntElem)
c  
      integer numnp,mgk
      integer Lnp(maxNumnp),Lmnp(mnsch,maxNumnp),
     &        Lgp(maxGP),Lmgp(mnsch,maxGP)
      integer igp_elemind(maxGP),igp_LocIntInd(maxGP),
     &        lnods(maxNode,maxElem)
c
      real*8 shpk(mnsch,maxGP),shpkdx(mnsch,maxGP),
     &       shpkdy(mnsch,maxGP),
     &       shpn(mnsch,maxNumnp)
c
c.....Local Variables............
c
      real*8  b(6),bd(2,6),bdd(3,6),
     &        cpt(2),cjt(2),dcjt(2),
     &        shp,shpd(2),shpdd(3),
     &        shapeFEM(maxNode),shpdevFEM(2,maxNode) 
      real*8  xxk,yyk,xxn,yyn,xr,yr,wjt,
     &        xsi,eta,teps
      integer ie,ik,ip,ipt,inode,jp,jpt,jLoop,nite,
     &        iLocint
c
      integer iMaterType,nnode,nintElem,iLumping,
     &        imethRKPM,imethFEM,imeth,iIntMeth,iInter,
     &        iGL,iAdapt,numadp,iad_stage,Istage
c
      integer i,n_pt,lnsp(mnsch)
      real*8  xbar(2)
      integer iCrack
      real*8  Cyplane,Cxmin,Cxmax
c
      common /connect/Lnp,Lmnp,Lgp,Lmgp,lnods
      common /shapeK/shpk,shpkdx,shpkdy,shpn
      common /mesh/xm,dxm,dvm,numnp
      common /ctrl/iMaterType,nnode,nintElem,iLumping
      common /meth/imethRKPM,imethFEM,imeth,iIntMeth,iInter
      common /adapt/iGL,iAdapt,numadp,iad_stage,Istage
      common /crack/Cyplane,Cxmin,Cxmax
      common /Type1/iCrack
c            
c-------------------------------------------------------------
c
c     Searching nodes for each gauss point
c
     	teps  = 1.99999999d0
        if (imeth .eq. imethRKPM) then
         do ik = 1, mgk
           xxk = xmk(1,ik)
           yyk = xmk(2,ik)
           nite = 0
           do jpt = 1, numnp
             xxn = xm(1,jpt)
             yyn = xm(2,jpt)
             xr = dabs((xxk-xxn)/dxm(1,jpt))
             yr = dabs((yyk-yyn)/dxm(2,jpt))
             if (xr .lt. teps .and. yr .le. teps) then
               nite = nite + 1
               if (nite .gt. mnsch) then
                  write(*,*)
     & 'nite > mnsch when searching nodes for each gauss point'
                  write(*,*) 'nite=',nite,'  jpt=',jpt
                  stop
               endif
               Lmgp(nite,ik) = jpt
             end if
           enddo
           Lgp(ik) = nite
         enddo
      elseif (imeth .eq. imethFEM) then
         do ik = 1,mgk
            ie = igp_elemind(ik)
            Lgp(ik) = nnode
            do inode = 1,nnode
               Lmgp(inode,ik)=lnods(inode,ie)
            enddo
         enddo
      endif   
c
c.....Modify connectivity map via visibility condition
c
      if ((iCrack .eq. 1).and.(imeth .eq. imethRKPM)) then
         do ik = 1, mgk
	      do i = 1,2
		 xbar(i) = xmk(i,ik)
              enddo
c
	      n_pt = Lgp(ik)
	      do ip = 1, n_pt
		 Lnsp(ip) = Lmgp(ip,ik)
              enddo
c
	  call visibleY2d(xbar,n_pt,Lnsp)
c
c.......Updated
c
	 Lgp(ik) = n_pt
	 do ip = 1, n_pt
	    Lmgp(ip,ik) = Lnsp(ip) 
	 enddo
c
	 enddo ! ik
       endif     ! iCrack = 1
c         
c
c     Searching nodes for each node
c
c 
      if (imeth .eq. imethRKPM) then      
         do ipt = 1, numnp
           xxk = xm(1,ipt)
           yyk = xm(2,ipt)
           nite = 0
           do jpt = 1, numnp
             xxn = xm(1,jpt)
             yyn = xm(2,jpt)
             xr = dabs((xxk-xxn)/dxm(1,jpt))
             yr = dabs((yyk-yyn)/dxm(2,jpt))
             if (xr .lt. teps .and. yr. lt. teps) then
               nite = nite + 1
               if (nite .gt. mnsch) then
                  write(*,*)
     &               'nite > mnsch when searching nodes for each node'
                  write(*,*) 'nite=',nite,'  jpt=',jpt
                  stop
               endif
c
               Lmnp(nite,ipt) = jpt
             end if
           enddo
           Lnp(ipt) = nite
         enddo
      elseif ( imeth.eq.imethFEM ) then
               ! 
               ! notice the Kronecker delta property of FEM shape function,
               ! the only node in the support of any node is itself
               !
         do ipt=1,numnp
            Lnp(ipt)   = 1
            Lmnp(1,ipt)= ipt
         enddo
c         
      endif
c            
      print *,'Searching is done'
c
c.....Modify connectivity map via visibility condition
c
      if ((iCrack.eq.1) .and. (imeth.eq.imethRKPM))then
         do ipt = 1, numnp
	      do i = 1,2
                 xbar(i) = xm(i,ipt)
              enddo
c
              n_pt  = Lnp(ipt)
	      do jp = 1, n_pt
		 Lnsp(jp) = Lmnp(jp,ipt)
              enddo
c
	 call visibleY2d(xbar,n_pt,Lnsp)
c
c.......Updated
c
	       Lnp(ipt) = n_pt
	       do jp = 1, n_pt
		  Lmnp(jp,ipt) = Lnsp(jp)
               enddo
c
	  enddo ! ipt
      endif     ! iCrack = 1
c
      print *,'Searching is done'
c
c--------------------------------------------------------------------
c
c  Eval and store Shape function NJ(xk): 
c  the shape functions at gaussion points in the support of any nodes
c
      do ik = 1, mgk
      
         if (imeth .eq. imethRKPM) then      
c      
            cpt(1) = xmk(1,ik)
            cpt(2) = xmk(2,ik)
c
	    jLoop = Lgp(ik)
	    do jp =1, jLoop
	       Lnsp(jp) = Lmgp(jp,ik)
	    enddo
c            
c ** to consider the different base function ------------------
c
            call correct(b,bd,bdd,cpt,xm,dxm,dvm,numnp,
     &                   iInter,jLoop,Lnsp)
c
c -------------------------------------------------------------
c            
            jLoop = Lgp(ik)
            do jp = 1, jLoop
               jpt     = Lmgp(jp,ik)
               cjt(1)  = xm(1,jpt)
               cjt(2)  = xm(2,jpt)
               dcjt(1) = dxm(1,jpt)
               dcjt(2) = dxm(2,jpt)
               wjt     = dvm(jpt)
c
c ** to consider the different base function ------------------
c
              call RKPMshape(shp,shpd,shpdd,b,bd,bdd,
     &                       cpt,cjt,dcjt,wjt,iInter)
c
c--------------------------------------------------------------
c
               shpk(jp,ik)   = shp
               shpkdx(jp,ik) = shpd(1)
               shpkdy(jp,ik) = shpd(2)
c
            enddo
c 
         elseif (imeth.eq.imethFEM) then
         
            ie     =igp_elemind(ik)
            iLocint=igp_LocintInd(ik)
c            
            xsi = gp_loc2D(1,iLocint)
            eta = gp_loc2D(2,iLocint)
c
            if (nnode.eq.4) then
               call evl_FEM_shpdev4(lnods,xm,xsi,eta,ie,
     &                              shapeFEM,shpdevFEM)
            elseif (nnode.eq.3) then
               call evl_FEM_shpdev3(lnods,xm,xsi,eta,ie,
     &                              shapeFEM,shpdevFEM)
            else
               write(*,*) 'In main: Unknown number nnode=',nnode
               stop
            endif
            do inode=1,nnode            
               shpk(inode,ik)   = shapeFEM(inode)
               shpkdx(inode,ik) = shpdevFEM(1,inode)
               shpkdy(inode,ik) = shpdevFEM(2,inode)   
            enddo
         endif 
      enddo   ! ik
c      
      write(*,*) 'cal shape function for Gauss points passed.'
c
c
c     Store Shape function NJ(xJ)
c     --------------------------
c
c       print *,'shape for x'

      if (imeth .eq. imethRKPM) then 
         do ipt = 1, numnp
            cpt(1) = xm(1,ipt)
            cpt(2) = xm(2,ipt)
c
	    jLoop = Lnp(ipt)
	    do jp = 1, jLoop
	       Lnsp(jp) = Lmnp(jp,ipt)
	    enddo
c
c-----------------------------------------------------------
c
c** consider different base function ------------------------
c
            call correct(b,bd,bdd,cpt,xm,dxm,dvm,numnp,
     &                   iInter,jLoop,Lnsp)
c
c------------------------------------------------------------
c
            do jp = 1, jLoop
               jpt     = Lmnp(jp,ipt)
               cjt(1)  = xm(1,jpt)
               cjt(2)  = xm(2,jpt)
               dcjt(1) = dxm(1,jpt)
               dcjt(2) = dxm(2,jpt)
               wjt     = dvm(jpt)
c              
c** consider different base function ------------------------              
c
              call RKPMshape(shp,shpd,shpdd,b,bd,bdd,
     &                       cpt,cjt,dcjt,wjt,iInter)
c
c------------------------------------------------------------
c
              shpn(jp,ipt) = shp
            enddo
         enddo
      elseif (imeth .eq. imethFEM) then
              !
              ! notice the Kronecker delta property of FEM shape function,
              ! the only node in the support of any node is itself
              ! and the shape function value of N_I (X_I ) = 1
              !                                    
         do ipt = 1,numnp
            shpn(1,ipt) = 1.0
         enddo
      endif 
         
      write(*,*) 'cal shape function for nodal points passed.'
c
c
c-----------------------------------------------------------------
 5987 format(1x,i5,1x,i5)
c
c      close(15)
c
      return
      end
c
