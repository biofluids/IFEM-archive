      subroutine hadapt(xm,dxm,dvm,xmk,dvmk,
     &           LmTract,vTract,
     &           LmDispbcX,LmDispbcY,vDispbcX,vDispbcY,
     &           LmVelbcX,LmVelbcY,vVelbcX,vVelbcY,
     &           igp_elemind,igp_LocintInd,gp_loc2D, ! for FEM
     &           nelem,numnp,mgk,
     &           effEnpt,ifwAdapt,fhead)
c
c***************************************************************
c
c     The currect h-adaptivity subroutine is only valid for
c     quadraliteral element, or quadraliteral particle distribution.
c
c     Documentation:
c
c     Output:
c
c     gp_loc2D:       the natural coordinate for FEM;
c     igp_LocintInd:  the global address of the element
c                     that the Gauss pt. contains;
c
c     Definitions of various connectivities:
c
c     Leff_ceL:   [ic:nelem]   |-->   [ic:icn_adp]; 
c
c     (icn_adp:   total number of cells that are needed refined;)
c    
c     Leff_ceG:   [ic:icn_adp] |-->   [ic:nelem]
c
c     IGbc:       [id:numnp]   |-->   [iDispbc:nDispbcX]
c                                     [iDispbc:nDispbcY]
c                                     [iVelbc :nVelbcX]
c                                     [iVelbc :nVelbcY]
c     
c     Leff_adp:  [ipt:numnp]   |-->   [iadp:numadp]
c
c     Note:
c
c     This is not a one-one maping; one cell can maps up to two
c     segments under current restriction, (in general, to all 
c     four segements, we exclude that possibility).
c     By assuming the monotonic connection, those two segments
c     are adjacent.
c
c     JTract:     [it:nTract]  |--> [iad:it_mark]
c
c
c      August, 1998
c
c
c****************************************************************
c
      implicit none
      include 'parameter.h'
c
      real*8  xm(2,maxNumnp),dxm(2,maxNumnp),dvm(maxNumnp),
     &        xmk(2,maxGP),dvmk(maxGP)
c
      integer LmDispbcX(maxDispbcX),LmDispbcY(maxDispbcY),
     &        LmVelbcX(maxVelbcX),LmVelbcY(maxVelbcY),
     &        LmTract(2,maxTract)
c
      real*8  vDispbcX(maxDispbcX),vDispbcY(maxDispbcY),
     &        vVelbcX(maxVelbcX),vVelbcY(maxVelbcY),
     &        vTract(2,2,maxTract)
c
      integer nelem,numnp,mgk,ifwAdapt
      real*8  gp_loc2D(2,maxIntElem)
      real*8  gp_weight2D(maxIntElem)
      integer igp_elemind(maxGP),igp_locintInd(maxGP)
c
c.....local array & variables...............
c
      real*8  eff_cut,effi,effj,cr_eff,eps0,eps1,rr,
     &        xn1,yn1,xn2,yn2,xn3,yn3,xn4,yn4,
     &        xn5,yn5,xn6,yn6,xn7,yn7,xn8,yn8,xn9,yn9,
     &        xi,yi,xj,yj,xsi,eta,ajj,dr,
     &        dvtt1,dvtt2
c
      real*8  shape(maxNode),
     &        effEopt(maxNumnp),
     &        effEnpt(maxNumnp)
c
      integer i,id5,id6,id7,id8,id9,jcm,
     &        iadp,ncount,ncount1,ipt,jpt,jp,jc,je,
     &        ic,ie,ik,ind,indn,jnd,icm,
     &        ncheck,ncheck1,ncheck2,
     &        kg,iLocint,n1,n2,n3,n4,
     &        m1,m2,m3,m4,jLoop,
     &        nelem_new,numnp_new,
     &        iDispbc,iVelbc,iT,it_mark,
     &        icn_adp,ic_mark,id_mark,ied
c
      integer Leff_adp(maxNumnp),Leff_adpG(maxNumnp),
     &        Leff_ceL(maxElem),Leff_ceG(maxElem), ! maxAdapt !
     &        IGbc(maxNumnp),
     &        JTract(maxTract)
c
      character*60 fhead,fwname,tmpstr1,tmpstr2
c
      real*8 ax,ay,dt,time,afact,afact1,afact2,
     &       r0,xc1,xc2,zimp,cr,cr1,cr2,dxn,dyn,
     &       det1,det2,det3,det4
c
      integer nstep,nOutputFrq,nPlotFrq,iplot,iplot_gp,
     &        imethRKPM,imethFEM,imeth,iIntMeth,iGL,
     &        iMaterType,nintElem,nnode,idim,inode,nnc,
     &        nnc1,nnc2,
     &        nTract,nDispbcX,nDispbcY,nVelbcX,nVelbcY,
     &        iInter,iAdapt,numadp,iad_stage,Istage
c
      integer lnods(maxNode,maxElem),
     &        Lnp(maxNumnp),Lmnp(mnsch,maxNumnp),
     &        Lgp(maxGP),Lmgp(mnsch,maxGP)
c
      common /rkpm/ax,ay,afact1,afact2,nnc1,nnc2
      common /stepT/dt,time
      common /step/ nstep,nOutputFrq,nPlotFrq,iplot,iplot_gp
      common /shear/r0,xc1,xc2,zimp,cr1,cr2
      common /ctrl/iMaterType,nnode,nintElem
      common /meth/imethRKPM,imethFEM,imeth,iIntMeth,iInter
      common /adapt/iGL,iAdapt,numadp,iad_stage,Istage
      common /bound/nTract,nDispbcX,nDispbcY,nVelbcX,nVelbcY
      common /connect/Lnp,Lmnp,Lgp,Lmgp,lnods
c
      real*8 xyloc_of_elem4(2,4)
      data ((xyloc_of_elem4(idim,inode),idim=1,2),inode=1,4)
     ./
     .  -1.,-1.,
     .   1.,-1.,
     .   1., 1.,
     .  -1., 1.
     ./
c
c  The h-adaptive scheme: 
c
c                            y
c                            ^
c                            |
c                            |
c                            |
c                            |
c                   4--------7--------3
c                   |        |        |
c                   |   4    |   3    |
c                   |        |        |
c                   8--------9--------6 -------------------> x
c                   |        |        |
c                   |   1    |   2    |
c                   |        |        |
c                   1--------5--------2
c
c
c--------------------------------------------------------
c
c  Goal: Leff_adp: the connectivity array for adaptivity;
c
c--------------------------------------------------------
c
c.......Initialization the connectivity array
c
      do ic = 1, maxElem
	 Leff_ceL(ic) = 0
	 Leff_ceG(ic) = 0
      enddo
c
      do ipt = 1, maxNumnp
	 IGbc(ipt)      = 0
	 Leff_adp(ipt)  = 0
	 Leff_adpG(ipt) = 0
	 dvm(ipt)       = 0.0
      enddo
c
      if (Istage .eq. 1) then
	 cr    = cr1
	 afact = afact1
	 nnc   = nnc1
      elseif(Istage .eq. 2) then
	 cr    = cr2
	 afact = afact2
	 nnc   = nnc2
      endif
c
      eps0  = 1.00d-14
      eps1  = 1.00d-8
c
c.....Evaluate the refinement criteria......
c     --------------------------------
c
      if (iGL .eq. 0) then ! Global J_2
c
           do ipt = 1, numnp
              effEopt(ipt) = effEnpt(ipt)
           enddo
c
	   call permute(effEopt,numnp,maxNumnp)
c
           eff_cut = effEopt(numadp)
	   iadp    = 0
c 
	   do 1111 ipt = 1, numnp
	      xi = xm(1,ipt)
	      yi = xm(2,ipt)
c
	      if (effEnpt(ipt) .ge. eff_cut) then 
		 iadp = iadp + 1
		 Leff_adp(ipt)   = iadp
		 Leff_adpG(iadp) = ipt
              endif
 1111      continue
c
	   dr     = numadp * Istage
           numadp = numadp + int(dr)  
c 
c........the local criteria
c 
       elseif(iGL  .eq. 1) then ! Local J_2
	   iadp = 0           ! assume that ax, ay > 1.0 !
           do 2222 ipt = 1, numnp
	      effi = effEnpt(ipt)
	      jLoop= Lnp(ipt)
	      do jp   = 1, jLoop
		 jpt    = Lmnp(jp,ipt)
                 effj   = effEnpt(jpt)
                 cr_eff = (effi - effj)/effi
c
                 if (cr_eff .ge. cr ) then
		     iadp            = iadp + 1
		     Leff_adp(ipt)   = iadp
		     Leff_adpG(iadp) = ipt
		     go to 2222
		 endif
	      enddo  !jpt
 2222       continue !ipt
c
       else
                   !nothing happen  
       endif       !iGL
c
c.........Output the adaptive mesh
c
      call appext(fhead,'adaptH',fwname)
      open(ifwAdapt,file=fwname)
      do ipt = 1, numnp
	  if(Leff_adp(ipt) .ne. 0) then
	     write(ifwAdapt,97)xm(1,ipt), xm(2,ipt)
	  endif
      enddo
      call flush(ifwAdapt)
      close(ifwAdapt)
c
       print *,'adapted points', iadp
c
c  -------------------------------------------
c  Find the cell that is needed to be refined
c       ( The following adaptive algorithm 
c         is only for quadrilateral element )
c  -------------------------------------------
c
c      icn_adp = 0
c      do  ic = 1, nelem
c	  xn1 = xm(1,lnods(1,ic))
c	  yn1 = xm(2,lnods(1,ic))
c	  xn2 = xm(1,lnods(2,ic))
c         yn2 = xm(2,lnods(2,ic))
c	  xn3 = xm(1,lnods(3,ic))
c	  yn3 = xm(2,lnods(3,ic))
c	  xn4 = xm(1,lnods(4,ic))
c	  yn4 = xm(2,lnods(4,ic))
c
c          ncount = 0
c
c          do jp  = 1, iadp
c	     jpt = Leff_adpG(jp)
c	     xi  = xm(1,jpt)
c	     yi  = xm(2,jpt)
c
c	     det1 = abs((xn2-xn1)*(yi-yn1) - (yn2-yn1)*(xi-xn1))
c            det2 = abs((xn3-xn2)*(yi-yn2) - (yn3-yn2)*(xi-xn2))
c            det3 = abs((xn4-xn3)*(yi-yn3) - (yn4-yn3)*(xi-xn3))
c	     det4 = abs((xn1-xn4)*(yi-yn4) - (yn1-yn4)*(xi-xn4))
c
c             if (det1 .le. eps1) then
c	        if (abs(xn2 - xn1) .gt. abs(yn2-yn1)) then
c                   if (xn2 .gt. xn1) then
c		      if ((xn1 .le. xi) .and. (xi .le. xn2)) then
c		         ncount = ncount + 1  
c		      endif
c		   elseif( xn1 .gt. xn2) then
c		      if ((xn2 .le. xi) .and. (xi .le. xn1)) then
c		         ncount = ncount + 1  
c		      endif
c                   endif
c		elseif(abs(yn2-yn1) .gt. abs(xn2 - xn1)) then
c		   if (yn2 .gt. yn1) then
c		      if ((yn1 .le. yi) .and. (yi .le. yn2)) then
c		           ncount = ncount + 1  
c		      endif
c		   elseif( yn1 .gt. yn2) then
c		      if ((yn2 .le. yi) .and. (yi .le. yn1)) then
c		           ncount = ncount + 1  
c		      endif
c		   endif
c		endif
c             elseif (det2 .le. eps1) then
c	        if (abs(xn3 - xn2) .gt. abs(yn3-yn2)) then
c                   if (xn3 .gt. xn2) then
c		      if ((xn2 .le. xi) .and. (xi .le. xn3)) then
c		         ncount = ncount + 1  
c		      endif
c		   elseif(xn2 .gt. xn3) then
c		      if ((xn3 .le. xi) .and. (xi .le. xn2)) then
c		         ncount = ncount + 1  
c		      endif
c                   endif
c                elseif(abs(yn3-yn2) .gt. abs(xn3 - xn2)) then
c		   if (yn3 .gt. yn2) then
c		      if ((yn2 .le. yi) .and. (yi .le. yn3)) then
c		           ncount = ncount + 1  
c		      endif
c		   elseif(yn2 .gt. yn3) then
c		      if ((yn3 .le. yi) .and. (yi .le. yn2)) then
c		          ncount = ncount + 1  
c		      endif
c		   endif
c		endif
c             elseif (det3 .le. eps1) then
c	        if (abs(xn4-xn3) .gt. abs(yn4-yn3)) then
c	           if (xn4 .gt. xn3) then
c		      if ((xn3 .le. xi) .and. (xi .le. xn4)) then
c		         ncount = ncount + 1  
c		      endif
c		   elseif(xn3 .gt. xn4) then
c		      if ((xn4 .le. xi) .and. (xi .le. xn3)) then
c		         ncount = ncount + 1  
c		      endif
c                   endif
c	        elseif(abs(yn4-yn3) .gt. abs(xn4-xn3)) then
c		   if (yn4 .gt. yn3) then
c		      if ((yn3 .le. yi) .and. (yi .le. yn4)) then
c		         ncount = ncount + 1  
c		      endif
c		   elseif(yn3 .gt. yn4) then
c		      if ((yn4 .le. yi) .and. (yi .le. yn3)) then
c		          ncount = ncount + 1  
c		      endif
c		   endif
c		endif
c             elseif (det4 .le. eps1) then
c	        if (abs(xn1-xn4) .gt. abs(yn1-yn4)) then
c	           if (xn1 .gt. xn4) then
c		      if ((xn4 .le. xi) .and. (xi .le. xn1)) then
c		         ncount = ncount + 1  
c		      endif
c		   elseif(xn4 .gt. xn1) then
c		      if ((xn1 .le. xi) .and. (xi .le. xn4)) then
c		         ncount = ncount + 1  
c		      endif
c                   endif
c	        elseif(abs(yn1-yn4) .gt. abs(xn1-xn4)) then
c	           if(yn1 .gt. yn4) then
c	              if ((yn4 .le. yi) .and. (yi .le. yn1)) then
c		         ncount = ncount + 1  
c		      endif
c		   elseif(yn4 .gt. yn1) then
c		      if ((yn1 .le. yi) .and. (yi .le. yn4)) then
c		         ncount = ncount + 1  
c		      endif
c		   endif
c		endif
c             endif
cc
c	    enddo  ! jp
c
c	    if(ncount .ge. nnc) then
c	       icn_adp           = icn_adp + 1
c	       Leff_ceL(ic)      = icn_adp   !! The Local address
c	       Leff_ceG(icn_adp) = ic        !! The Global address
c            endif
c	 enddo     !ic
cc
c
c
	 icn_adp = 0
         do ic  = 1, nelem
	    ncount = 0                          !!Warning misuse symbol !!
c
	    do jnd = 1, 4
	       jpt = lnods(jnd,ic)
	       if (Leff_adp(jpt) .ne. 0) then   !! count how many pt. in 
                  ncount = ncount + 1           !! a cell that are marked 
	       endif
            enddo
c
	    if(ncount .ge. nnc) then
	       icn_adp           = icn_adp + 1
	       Leff_ceL(ic)      = icn_adp      !! The Local address
	       Leff_ceG(icn_adp) = ic           !! The Global address
            endif
         enddo 
c
	 print *, 'adapted cell', icn_adp
	 nelem_new = nelem + 3 * icn_adp
c 
c        old element number : nelem - icn_adp
c 
c.......unfinished ??
c
         id_mark = numnp
	 ic_mark = nelem
	 do 3333 ic = 1, icn_adp
            icm = ic - 1
	    ie  = Leff_ceG(ic)
c
            n1 = lnods(1,ie)
            n2 = lnods(2,ie)
            n3 = lnods(3,ie)
            n4 = lnods(4,ie)
c
            xn1 = xm(1,n1)
	    yn1 = xm(2,n1)
	    xn2 = xm(1,n2)
	    yn2 = xm(2,n2)
	    xn3 = xm(1,n3)
            yn3 = xm(2,n3)
	    xn4 = xm(1,n4)
	    yn4 = xm(2,n4)
c
            xn5 = (xn1 + xn2)/2.0d0
            yn5 = (yn1 + yn2)/2.0d0
            xn6 = (xn2 + xn3)/2.0d0
            yn6 = (yn2 + yn3)/2.0d0
            xn7 = (xn3 + xn4)/2.0d0
            yn7 = (yn3 + yn4)/2.0d0 
            xn8 = (xn4 + xn1)/2.0d0
            yn8 = (yn4 + yn1)/2.0d0
c
            xn9 = (xn1 + xn2 + xn3 + xn4)/4.0d0
            yn9 = (yn1 + yn2 + yn3 + yn4)/4.0d0
c
c..........Redefine the dilation parameter at nodal pt.
c
            dxn = dabs(xn2 - xn1)
            dyn = dabs(yn2 - yn1)
c
	    if (dxn .lt. dabs(xn3-xn2)) then
               dxn = dabs(xn3 - xn2)
	    endif
c
	    if(dxn .lt. dabs(xn4-xn3)) then
               dxn = dabs(xn4 - xn3)
            endif
c
	    if(dxn .lt. dabs(xn1-xn4)) then
               dxn = dabs(xn1 - xn4)
            endif
c
c................................................
c
	    if (dyn .lt. dabs(yn3-yn2)) then
               dyn = dabs(yn3 - yn2)
	    endif
c
	    if(dyn .lt. dabs(yn4-yn3)) then
               dyn = dabs(yn4 - yn3)
	    endif
c
	    if(dyn .lt. dabs(yn1-yn4)) then
               dyn = dabs(yn1 - yn4)
            endif
c
	    dxn = afact*ax*dxn
	    dyn = afact*ay*dyn
c
c
c..........previous strategy ......................
c
c           dxn = afact*dxm(1,1)
c           dyn = afact*dxm(2,1)
c
c            do i = 2,4
c	       ind = lnods(i,ie)
c	       if (dxn .lt. afact*dxm(1,ind)) then
c		  dxn = afact*dxm(1,ind)
c	       endif
cc
c	       if (dyn .lt. afact*dxm(2,ind)) then
c		  dyn = afact*dxm(2,ind)
c	       endif
c	    enddo
c
c....................................................
c
            if(dxm(1,n1) .gt. dxn) then
               dxm(1,n1) = dxn
            endif               ! the judgement is important
c                               ! for the mixed size refinement
            if(dxm(2,n1) .gt. dyn) then
	       dxm(2,n1) = dyn
            endif
c
            if(dxm(1,n2) .gt. dxn) then
	       dxm(1,n2) = dxn
            endif
c
	    if(dxm(2,n2) .gt. dyn) then
	       dxm(2,n2) = dyn
            endif
c
            if(dxm(1,n3) .gt. dxn) then
	       dxm(1,n3) = dxn
            endif
c
	    if(dxm(2,n3) .gt. dyn) then
	       dxm(2,n3) = dyn
            endif
c
            if(dxm(1,n4) .gt. dxn) then
	       dxm(1,n4) = dxn
            endif
c
            if(dxm(2,n4) .gt. dyn) then
	       dxm(2,n4) = dyn
            endif
c
c
c............Add the new particles:
c----------------------------------------------------
c
c...........I...............
c
	    if (Istage .ge. 2) then
               jLoop = Lnp(n1)
	       do jp = 1, jLoop
	          jpt = Lmnp(jp,n1)
	          xj  = xm(1,jpt)
	          yj  = xm(2,jpt)
	          rr  = (xn5 - xj)**2 + (yn5 - yj)**2
	          rr  = dsqrt(rr)
	          if (rr .le. eps1) then
		     id5  = jpt
		     go to 3341          ! the point exists
	          endif                  ! at last refinement
	       enddo
	    endif
c
c
	    if (icm .lt. 1) go to 3340
            do jc = 1, icm
	       je = Leff_ceG(jc)
	       do ind = 1,4
		  indn  = ind + 1  - int(ind/4)*4 
	          m1    = lnods(ind,je)
		  m2    = lnods(indn,je)
                  xj    = (xm(1,m1) + xm(1,m2))/2.0d0
		  yj    = (xm(2,m1) + xm(2,m2))/2.0d0
                  rr    = (xn5-xj)**2 + (yn5-yj)**2
		  rr    = dsqrt(rr)
c
                  if(rr .le. eps1) then
	            jcm= nelem + 4*(jc - 1) + ind
		    id5 = lnods(indn,jcm)
		    go to 3341            ! the point exists
                  endif                   ! at current refinement
	       enddo
            enddo ! jc
c
c
 3340       continue
	    id_mark    = id_mark + 1
	    id5        = id_mark
	    xm(1,id5)  = xn5
	    xm(2,id5)  = yn5
c
            dxm(1,id5) = dxn
            dxm(2,id5) = dyn
c
 3341       continue
c
c...........II.............
c
	    if(Istage .ge. 2) then
               jLoop = Lnp(n2)
	       do jp = 1, jLoop
	          jpt = Lmnp(jp,n2)
	          xj  = xm(1,jpt)
	          yj  = xm(2,jpt)
	          rr  = (xn6 - xj)**2 + (yn6 - yj)**2
	          rr  = dsqrt(rr)
	          if (rr .le. eps1) then
		      id6 = jpt
		      go to 3351 ! the point exists
	          endif          ! at last refinement
	       enddo             
            endif
c
	    if (icm .lt. 1) go to 3350
            do jc = 1, icm
	       je = Leff_ceG(jc)
	       do ind = 1,4
		  indn  = ind + 1  - int(ind/4)*4 
	          m1  = lnods(ind,je)
		  m2  = lnods(indn,je) 
                  xj  = (xm(1,m1) + xm(1,m2))/2.0d0
		  yj  = (xm(2,m1) + xm(2,m2))/2.0d0
                  rr  = (xn6-xj)**2 + (yn6-yj)**2
		  rr  = dsqrt(rr)
c
                  if(rr .le. eps1) then
	            jcm = nelem + 4*(jc - 1) + ind
		    id6 = lnods(indn,jcm)
		    go to 3351            ! the point exists
                  endif                   ! at current refinement
	       enddo
            enddo ! jc
c
c
 3350       continue
	    id_mark    = id_mark + 1
	    id6        = id_mark
	    xm(1,id6)  = xn6
	    xm(2,id6)  = yn6
c
            dxm(1,id6) = dxn
            dxm(2,id6) = dyn
c
c
 3351       continue
c
c...................III......................
c
	    if (Istage .eq. 2) then
               jLoop = Lnp(n3)
	       do jp = 1, jLoop
	          jpt = Lmnp(jp,n3)
	          xj  = xm(1,jpt)
	          yj  = xm(2,jpt)
	          rr  = (xn7-xj)**2 + (yn7-yj)**2
	          rr  = dsqrt(rr)
	          if (rr .le. eps1) then
		     id7  = jpt
		     go to 3361      ! the point exists
	          endif              ! at the last refinement
	       enddo
            endif
c
            if (icm .lt. 1) go to 3360
            do jc = 1, icm
	       je = Leff_ceG(jc)
	       do ind = 1,4
		  indn  = ind + 1 - int(ind/4)*4 
	          m1    = lnods(ind,je)
		  m2    = lnods(indn,je) 
                  xj    = (xm(1,m1) + xm(1,m2))/2.0d0
		  yj    = (xm(2,m1) + xm(2,m2))/2.0d0
                  rr    = (xn7-xj)**2 + (yn7-yj)**2
		  rr    = dsqrt(rr)
c
                  if(rr .le. eps1) then
	            jcm= nelem + 4*(jc - 1) + ind
		    id7 = lnods(indn,jcm)
		    go to 3361            ! the point exists
                  endif                   ! at current refinement
	       enddo
            enddo ! jc
c
c
 3360       continue
	    id_mark     = id_mark + 1
	    id7         = id_mark
	    xm(1,id7)   = xn7
	    xm(2,id7)   = yn7
c
            dxm(1,id7)  = dxn
            dxm(2,id7)  = dyn
c
c.....................IV.....................
c
 3361       continue
	    if(Istage .eq. 2) then
              jLoop = Lnp(n4)
	      do jp = 1, jLoop
	         jpt = Lmnp(jp,n4)
	         xj  = xm(1,jpt)
	         yj  = xm(2,jpt)
	         rr  = (xn8-xj)**2 + (yn8-yj)**2
	         rr  = dsqrt(rr)
	         if (rr .le. eps1) then
		    id8 = jpt
		    go to 3371         ! the point exists
	         endif                 ! at the last refinement
	      enddo                    ! jp
            endif                      ! Istage
c
	    if (icm .lt. 1) go to 3370
            do jc = 1, icm
	       je = Leff_ceG(jc)
	       do ind = 1,4
		  indn  = ind + 1  - int(ind/4)*4 
	          m1    = lnods(ind,je)
		  m2    = lnods(indn,je) 
                  xj    = (xm(1,m1) + xm(1,m2))/2.0d0
		  yj    = (xm(2,m1) + xm(2,m2))/2.0d0
                  rr    = (xn8-xj)**2 + (yn8-yj)**2
		  rr    = dsqrt(rr)
c
                  if(rr .le. eps1) then
	            jcm= nelem + 4*(jc - 1) + ind
		    id8 = lnods(indn,jcm)
		    go to 3371            ! the point exists
                  endif                   ! at current refinement
	       enddo
	     enddo !jc
c
 3370       continue
	    id_mark     = id_mark + 1
	    id8         = id_mark
	    xm(1,id8)   = xn8
	    xm(2,id8)   = yn8
c
            dxm(1,id8) = dxn
            dxm(2,id8) = dyn
c
c
c.....................V......................
c
 3371       continue
	    id_mark = id_mark + 1
	    id9     = id_mark
	    xm(1,id9)  = xn9
	    xm(2,id9)  = yn9
	    dxm(1,id9) = dxn
	    dxm(2,id9) = dyn
c
c..........the old nodal point still use the previous
c
c
c...........add the connectivity in a dump way
c          (which is redundant at the moment)
c
	       icm = ic_mark + 1
	       lnods(1,icm) = lnods(1,ie)  !(pt. 1)
	       lnods(2,icm) = id5          !(pt. 5)  !!Problem ???
	       lnods(3,icm) = id9          !(pt. 9)
	       lnods(4,icm) = id8          !(pt. 8)
c   
	       icm = ic_mark + 2
	       lnods(1,icm) = id5          !(pt. 5)
	       lnods(2,icm) = lnods(2,ie)  !(pt. 2)
	       lnods(3,icm) = id6          !(pt. 6)
	       lnods(4,icm) = id9          !(pt. 9)
c
	       icm = ic_mark + 3
	       lnods(1,icm) = id9          !(pt. 9)
	       lnods(2,icm) = id6          !(pt. 6)
	       lnods(3,icm) = lnods(3,ie)  !(pt. 3)
	       lnods(4,icm) = id7          !(pt. 7)
c
	       icm = ic_mark + 4           
               lnods(1,icm) = id8          !(pt. 8)
	       lnods(2,icm) = id9          !(pt. 9)
               lnods(3,icm) = id7          !(pt. 7)
               lnods(4,icm) = lnods(4,ie)  !(pt. 4)
c
               ic_mark = icm
c
 3333     continue
c
c
	  if (ic_mark .eq. nelem_new + icn_adp) then
              print *, 'nodal adaptivity is ready'
          else
              print *, 'something is wrong'
	      stop
	  endif
c
	  numnp_new = id_mark
c
c---------------------------------------------
c.........Update the boundary connectivity    
c---------------------------------------------
c
       if (nDispbcX .eq. 0) go to 3390
       do iDispbc = 1, nDispbcX
          ncheck       = LmDispbcX(iDispbc) ! release global address
	  IGbc(ncheck) = iDispbc            ! store local address
       enddo
c
       ncount = 0 
       do ic = 1, icn_adp      ! icn_adp: total number 
          ie = Leff_ceG(ic) ! of cell are refined;
c
          n1 = lnods(1,ie)
          n2 = lnods(2,ie)
          n3 = lnods(3,ie)
          n4 = lnods(4,ie)
c
c------------------------------------------------
c
	  m1 = IGbc(n1)
	  m2 = IGbc(n2)
	  m3 = IGbc(n3)
	  m4 = IGbc(n4)
c
	  icm = ic - 1
	  if(m1 .ne. 0 .and. m2 .ne. 0) then
              ncount = ncount + 1
	      ind    = nelem + 4*icm + 1
	      ied    = lnods(2,ind)
              LmDispbcX(nDispbcX+ncount) = ied
              vDispbcX(nDispbcX + ncount)  = (vDispbcX(m1)
     &                                     +  vDispbcX(m2))/2.0d0
          endif
c
          if(m2 .ne. 0 .and. m3 .ne. 0) then
              ncount = ncount + 1
	      ind    = nelem + 4*icm + 2
	      ied    = lnods(3,ind)
              LmDispbcX(nDispbcX+ncount)  = ied
              vDispbcX(nDispbcX + ncount) = (vDispbcX(m2)
     &                                    +  vDispbcX(m3))/2.0d0
          endif
c
         if(m3 .ne. 0 .and. m4 .ne. 0) then
              ncount = ncount + 1
	      ind    = nelem + 4*icm + 3
	      ied    = lnods(4,ind)
              LmDispbcX(nDispbcX+ncount)  = ied
              vDispbcX(nDispbcX + ncount) = (vDispbcX(m3)
     &                                    +  vDispbcX(m4))/2.0d0
         endif
c
         if(m4 .ne. 0 .and. m1 .ne. 0) then
              ncount = ncount + 1
	      ind    = nelem + 4*icm + 4
	      ied    = lnods(1,ind)
              LmDispbcX(nDispbcX+ncount)  = ied
              vDispbcX(nDispbcX + ncount) = (vDispbcX(m4)
     &                                    +  vDispbcX(m1))/2.0d0
	 endif
        enddo
c
	nDispbcX = nDispbcX + ncount
c
 3390   continue
c...............................................................
c
       do ipt = 1, maxNumnp
	  IGbc(ipt) = 0
       enddo
c
       if(nDispbcY .eq. 0) go to 3392
       ncount = 0 
       do iDispbc = 1, nDispbcY
	  ncheck       = LmDispbcY(iDispbc)
	  IGbc(ncheck) = iDispbc
       enddo
c
       ncount = 0 
       do ic = 1, icn_adp
          ie = Leff_ceG(ic)
c
          n1 = lnods(1,ie)
          n2 = lnods(2,ie)
          n3 = lnods(3,ie)
          n4 = lnods(4,ie)
c
c-------------------------------------
c
	  m1 = IGbc(n1)
	  m2 = IGbc(n2)
	  m3 = IGbc(n3)
	  m4 = IGbc(n4)
c
	  icm = ic - 1
	  if(m1 .ne. 0 .and. m2 .ne. 0) then
              ncount = ncount + 1
	      ind    = nelem  + 4*icm + 1
	      ied    = lnods(2,ind)
              LmDispbcY(nDispbcY + ncount) = ied
              vDispbcY(nDispbcY  + ncount) = (vDispbcY(m1)
     &                                     +  vDispbcY(m2))/2.0d0
          endif
c
          if(m2 .ne. 0 .and. m3 .ne. 0) then
              ncount = ncount + 1
	      ind    = nelem  + 4*icm + 2
	      ied    = lnods(3,ind)
              LmDispbcY(nDispbcY + ncount) = ied
              vDispbcY(nDispbcY + ncount)  = (vDispbcY(m2)
     &                                     +  vDispbcY(m3))/2.0d0
          endif
c
          if(m3 .ne. 0 .and. m4 .ne. 0) then
              ncount = ncount + 1
	      ind    = nelem  + 4*icm + 3
	      ied    = lnods(4,ind)
              LmDispbcY(nDispbcY + ncount) = ied 
              vDispbcY(nDispbcY + ncount)  = (vDispbcY(m3)
     &                                     +  vDispbcY(m4))/2.0d0
          endif
c
          if(m4 .ne. 0 .and. m1 .ne. 0) then
             ncount = ncount + 1
	     ind    = nelem  + 4*icm + 4
	     ied    = lnods(1,ind)
             LmDispbcY(nDispbcY + ncount) = ied
             vDispbcY(nDispbcY  + ncount) = (vDispbcY(m4)
     &                                    +  vDispbcY(m1))/2.0d0
	  endif
        enddo ! ic 
c
        nDispbcY = nDispbcY + ncount
c
 3392   continue
c...............................................................
c
       do ipt = 1, maxNumnp
	  IGbc(ipt) = 0
       enddo
c
       if(nVelbcX .eq. 0) go to 3394
       do iVelbc = 1, nVelbcX
	  ncheck       = LmVelbcX(iVelbc)
	  IGbc(ncheck) = iVelbc
       enddo
c
       ncount = 0 
       do ic = 1, icn_adp
          ie = Leff_ceG(ic)
c
          n1 = lnods(1,ie)
          n2 = lnods(2,ie)
          n3 = lnods(3,ie)
          n4 = lnods(4,ie)
c
c-----------------------------------------
c
	  m1 = IGbc(n1)
	  m2 = IGbc(n2)
	  m3 = IGbc(n3)
	  m4 = IGbc(n4)
c
	  icm = ic - 1
	  if(m1 .ne. 0 .and. m2 .ne. 0) then
              ncount = ncount + 1
	      ind    = nelem + 4*icm + 1
	      ied    = lnods(2,ind)
              LmVelbcX(nVelbcX + ncount) = ied
              vVelbcX(nVelbcX + ncount)  = (vVelbcX(m1)
     &                                   +  vVelbcX(m2))/2.0d0
          endif
c
          if(m2 .ne. 0 .and. m3 .ne. 0) then
              ncount = ncount + 1
	      ind    = nelem + 4*icm + 2
	      ied    = lnods(3,ind)
              LmVelbcX(nVelbcX+ncount)  = ied
              vVelbcX(nVelbcX + ncount) = (vVelbcX(m2)
     &                                  +  vVelbcX(m3))/2.0d0
          endif
c
          if(m3 .ne. 0 .and. m4 .ne. 0) then
              ncount = ncount + 1
	      ind    = nelem + 4*icm + 3
	      ied    = lnods(4,ind)
              LmVelbcX(nVelbcX+ncount)  = ied
              vVelbcX(nVelbcX + ncount) = (vVelbcX(m3)
     &                                  +  vVelbcX(m4))/2.0d0
          endif
c
	  if(m4 .ne. 0 .and. m1 .ne. 0) then
              ncount = ncount + 1
	      ind    = nelem + 4*icm + 4
	      ied    = lnods(1,ind)
              LmVelbcX(nVelbcX+ncount)  = ied
              vVelbcX(nVelbcX + ncount) = (vVelbcX(m1)
     &                                  +  vVelbcX(m4))/2.0d0
	  endif
        enddo ! ic
c
	nVelbcX = nVelbcX + ncount
c
 3394   continue
c
c...............................................................
c
       do ipt = 1, maxNumnp
	  IGbc(ipt) = 0
       enddo
c
       if(nVelbcY .eq. 0) go to 3396
       do iVelbc = 1, nVelbcY
          ncheck       = LmVelbcY(iVelbc)
	  IGbc(ncheck) = iVelbc
       enddo
c
       ncount = 0 
       do ic = 1, icn_adp
          ie = Leff_ceG(ic)
c
          n1 = lnods(1,ie)
          n2 = lnods(2,ie)
          n3 = lnods(3,ie)
          n4 = lnods(4,ie)
c
c--------------------------------------------------------
c
	  m1 = IGbc(n1)
	  m2 = IGbc(n2)
	  m3 = IGbc(n3)
	  m4 = IGbc(n4)
c
	  icm = ic - 1
	  if( m1 .ne. 0 .and. m2 .ne. 0) then
              ncount = ncount + 1
	      ind    = nelem  + 4*icm + 1
	      ied    = lnods(2,ind)
              LmVelbcY(nVelbcY + ncount) = ied
              vVelbcY(nVelbcY + ncount)  = (vVelbcY(m1)
     &                                   +  vVelbcY(m2))/2.0d0
c
          endif
c
          if(m2 .ne. 0 .and. m3 .ne. 0) then
               ncount = ncount + 1
	       ind    = nelem  + 4*icm + 2
	       ied    = lnods(3,ind)
               LmVelbcY(nVelbcY + ncount) = ied
               vVelbcY(nVelbcY + ncount)  = (vVelbcY(m2)
     &                                    +  vVelbcY(m3))/2.0d0
          endif
c
          if(m3 .ne. 0 .and. m4 .ne. 0) then
               ncount = ncount + 1
	       ind    = nelem  + 4*icm + 3
	       ied    = lnods(4,ind)
               LmVelbcY(nVelbcY+ncount) = ied
               vVelbcY(nVelbcY + ncount)= (vVelbcY(m3)
     &                                  +  vVelbcY(m4))/2.0d0
          endif
c
	  if(m4 .ne. 0 .and. m1 .ne. 0) then
               ncount = ncount + 1
	       ind    = nelem  + 4*icm + 4
	       ied    = lnods(1,ind)
               LmVelbcY(nVelbcY+ncount) = ied
               vVelbcY(nVelbcY + ncount)= (vVelbcY(m1)
     &                                  +  vVelbcY(m4))/2.0d0
	  endif
       enddo
c
	nVelbcY = nVelbcY + ncount
c
 3396  continue
c
c      ! End of the Essential B.C. update
c
c
c--------------------------------------------------------
c
c............Update the traction B.C. ...............
c
c   Note:
c
c   The following algorithm is only suitable for
c   quadrilateral elements in a single connected
c   region, in which, both height and width dimension
c   contains more than one cell.
c   The traction boundary segments are numbered
c   in a counterclock-wise fashion.
c
c--------------------------------------------------------
c
      if (nTract .eq. 0) go to 6667
c
	 ncount = 0
	 ncount1= 0
         do iT = 1, nTract
	    ncheck1 = LmTract(1,iT)
	    ncheck2 = LmTract(2,iT)
c
	    do ic = 1, icn_adp
	       ie = Leff_ceG(ic)
	       icm= ic - 1 
               n1 = lnods(1,ie)
               n2 = lnods(2,ie)
               n3 = lnods(3,ie)
               n4 = lnods(4,ie)
c
	       if(n1 .eq. ncheck1 .and. n2 .eq. ncheck2) then
                  m1 = nTract + ncount + 1 
                  m2 = nTract + ncount + 2 
c
		  ind= nelem  + 4*icm + 1
                  LmTract(1,m1) = lnods(1,ind) 
                  LmTract(2,m1) = lnods(2,ind)
c
                  LmTract(1,m2) = lnods(1,ind + 1)
                  LmTract(2,m2) = lnods(2,ind + 1)
c
	          vTract(1,1,m1) = vTract(1,1,iT)
	          vTract(2,1,m1) = vTract(2,1,iT)
c
	          vTract(1,2,m1) = 0.5*(vTract(1,1,iT)
     &                           +      vTract(1,2,iT))
	          vTract(2,2,m1) = 0.5*(vTract(2,1,iT)
     &                           +      vTract(2,2,iT))
c
	          vTract(1,1,m2) = vTract(1,2,m1)
	          vTract(2,1,m2) = vTract(2,2,m1)
c
	          vTract(1,2,m2) = vTract(1,2,iT)
	          vTract(2,2,m2) = vTract(2,2,iT)
c
		  ncount  = ncount  + 2
		  ncount1 = ncount1 + 1
	       elseif(n2 .eq. ncheck1 .and. n3 .eq. ncheck2) then
                  m1 = nTract + ncount + 1 
                  m2 = nTract + ncount + 2 
c
		  ind= nelem + 4*icm + 2
                  LmTract(1,m1) = lnods(1,ind) 
                  LmTract(2,m1) = lnods(2,ind)
c
                  LmTract(1,m2) = lnods(1,ind + 1)
                  LmTract(2,m2) = lnods(2,ind + 1)
c
	          vTract(1,1,m1) = vTract(1,1,iT)
	          vTract(2,1,m1) = vTract(2,1,iT)
c
	          vTract(1,2,m1) = 0.5*(vTract(1,1,iT)
     &                           +      vTract(1,2,iT))
	          vTract(2,2,m1) = 0.5*(vTract(2,1,iT)
     &                           +      vTract(2,2,iT))
c
	          vTract(1,1,m2) = vTract(1,2,m1)
	          vTract(2,1,m2) = vTract(2,2,m1)
c
	          vTract(1,2,m2) = vTract(1,2,iT)
	          vTract(2,2,m2) = vTract(2,2,iT)
c
		  ncount  = ncount  + 2
		  ncount1 = ncount1 + 1
	       elseif(n3 .eq. ncheck1 .and. n4 .eq. ncheck2) then
                  m1 = nTract + ncount + 1 
                  m2 = nTract + ncount + 2 
c
		  ind= nelem + 4*icm + 3
                  LmTract(1,m1) = lnods(1,ind) 
                  LmTract(2,m1) = lnods(2,ind)
c
                  LmTract(1,m2) = lnods(1,ind + 1)
                  LmTract(2,m2) = lnods(2,ind + 1)
c
	          vTract(1,1,m1) = vTract(1,1,iT)
	          vTract(2,1,m1) = vTract(2,1,iT)
c
	          vTract(1,2,m1) = 0.5*(vTract(1,1,iT)
     &                           +      vTract(1,2,iT))
	          vTract(2,2,m1) = 0.5*(vTract(2,1,iT)
     &                           +      vTract(2,2,iT))
c
	          vTract(1,1,m2) = vTract(1,2,m1)
	          vTract(2,1,m2) = vTract(2,2,m1)
c
	          vTract(1,2,m2) = vTract(1,2,iT)
	          vTract(2,2,m2) = vTract(2,2,iT)
c
		  ncount  = ncount  + 2
		  ncount1 = ncount1 + 1
	       elseif(n4 .eq. ncheck1 .and. n1 .eq. ncheck2) then
                  m1  = nTract + ncount + 1 
                  m2  = nTract + ncount + 2 
c
		  ind = nelem + 4*icm + 4
                  LmTract(1,m1) = lnods(1,ind) 
                  LmTract(2,m1) = lnods(2,ind)
c
                  LmTract(1,m2) = lnods(1,ind - 3)
                  LmTract(2,m2) = lnods(2,ind - 3)
c
	          vTract(1,1,m1) = vTract(1,1,iT)
	          vTract(2,1,m1) = vTract(2,1,iT)
c
	          vTract(1,2,m1) = 0.5*(vTract(1,1,iT)
     &                           +      vTract(1,2,iT))
	          vTract(2,2,m1) = 0.5*(vTract(2,1,iT)
     &                           +      vTract(2,2,iT))
c
	          vTract(1,1,m2) = vTract(1,2,m1)
	          vTract(2,1,m2) = vTract(2,2,m1)
c
	          vTract(1,2,m2) = vTract(1,2,iT)
	          vTract(2,2,m2) = vTract(2,2,iT)
c
		  ncount  = ncount  + 2
		  ncount1 = ncount1 + 1
               endif
	    enddo    ! ic
	    JTract(iT) = ncount1
	 enddo       ! iT
c
c
c####
 6666   continue
c
c
 6667   continue
c
        it_mark = ncount1
c
c.......Shift the connectivity array lnodes
c
	ncount = 0
        do 7777 ic = 1, ic_mark
	   if(Leff_ceL(ic) .ne. 0) then
              ncount = Leff_ceL(ic) 
	      go to 7777    ! Do not shift, or, discard !!
	   endif 
	   ie = ic - ncount ! Shift a position by ``ncount''
	   do jnd = 1,4
	      lnods(jnd,ie) = lnods(jnd,ic)
           enddo
 7777   continue
c
	ncount = 0
        do 8888 iT = 1, it_mark
	   if(JTract(iT) .ne. 0) then
              ncount = JTract(iT) 
	      go to 8888    ! Do not shift, or, discard !!
	   endif 
	   ie = iT - ncount ! Shift a position by ``ncount''
	   do ind = 1,2
	      LmTract(ind,ie) = LmTract(ind,iT)
	      do jpt = 1,2
		 vTract(ind,jpt,ie)  = vTract(ind,jpt,iT)
	      enddo
           enddo
 8888  continue
c 
c 
c------------------------------------------------------
c     Update xmk: the Gauss quadrature pt. coordinate
c------------------------------------------------------
c
        kg    = 0
	dvtt1 = 0.0d0
	dvtt2 = 0.0d0
c
        call gauss2D(gp_loc2D,gp_weight2D,nintElem,nnode)
c
        do ie = 1, nelem_new
	     xn1 = xm(1,lnods(1,ie))
	     yn1 = xm(2,lnods(1,ie))
	     xn2 = xm(1,lnods(2,ie))
	     yn2 = xm(2,lnods(2,ie))
	     xn3 = xm(1,lnods(3,ie))
	     yn3 = xm(2,lnods(3,ie))
	     xn4 = xm(1,lnods(4,ie))
	     yn4 = xm(2,lnods(4,ie))
c
	     do iLocint=1,nintElem
		kg  = kg + 1
		xsi = gp_loc2D(1,iLocint)
		eta = gp_loc2D(2,iLocint)
c
                if (imeth.eq.imethFEM) then
                   igp_elemind(kg)   = ie
		   igp_LocintInd(kg) = iLocint
                endif !FEM
c
		call evl_FEM_shape4(shape,xsi,eta)
		xmk(1,kg)=shape(1)*xn1+shape(2)*xn2
     &                   +shape(3)*xn3+shape(4)*xn4
                xmk(2,kg)=shape(1)*yn1+shape(2)*yn2
     &                 +shape(3)*yn3+shape(4)*yn4
		call ajacob4(ajj,xn1,xn2,xn3,xn4,
     &                       yn1,yn2,yn3,yn4,
     &                       xsi,eta)
c
                dvmk(kg) = ajj*gp_weight2D(iLocint)
		dvtt2    = dvtt2 + dvmk(kg)
c
		if (dvmk(kg) .le. eps0 ) then
		   print *,'Guass weight is wrong!'
		   stop
		endif
             enddo ! iLocint
c
	  do jnd = 1, nnode
	     jpt = lnods(jnd,ie)
	     call ajacob4(ajj,xn1,xn2,xn3,xn4,yn1,yn2,
     &                    yn3,yn4,xyloc_of_elem4(1,jnd),
     &                    xyloc_of_elem4(2,jnd))
	     dvtt1    = dvtt1 + ajj
	     dvm(jpt) = dvm(jpt) + ajj
	  enddo
      enddo ! endLoop for ie       
c
c.....Update the total number of Gauss point
c
      mgk = kg
      ind = id_mark - numnp
      print *, 'mgk =', mgk
      print *, 'mass-first', dvtt1,dvtt2
c
c------------------------------------------
c..........Output the adaptivity mesh......
c------------------------------------------
c
      ifwAdapt = ifwAdapt
      if (1 .eq. 1) then
         call Int2Str(Istage,tmpstr1)
	 call appstr('.Amesh',tmpstr1,tmpstr2)
	 call appstr(fhead,tmpstr2,fwname)
         write(*,*) 'fwname=',fwname
	 open(ifwAdapt,file=fwname)
c
	 write(ifwAdapt,'(1x,a)') 'VARIABLES="X","Y"'
	 write(ifwAdapt,'(1x,a,i6,a,i6,a)')
     &        'ZONE N = ',numnp_new,', E=',nelem_new,
     &        ', F = FEPOINT, ET = QUADRILATERAL'
c
         do ipt = 1, numnp_new
	     write(ifwAdapt,9200) xm(1,ipt),xm(2,ipt) !! unit wrong !!
         enddo
c
	 do ie = 1,nelem_new
	    write(ifwAdapt,'(1x,4i8)')
     &      (lnods(inode,ie), inode = 1, nnode)
	 enddo
	 close (ifwAdapt)
      endif
      if(1 .eq. 1) then
         call appext(fhead,'Afmesh',fwname)
	 open (ifwAdapt,file=fwname)
c##
c	      nelem_old = nelem - icn_adp
c	      nelem_new = nelem + 3 * icn_adp
c##
         do ipt = numnp +1, id_mark
	     write(ifwAdapt,9201) ipt, dxm(1,ipt),dxm(2,ipt) !! unit wrong !!
         enddo
	 close (ifwAdapt)
      endif
      if(1 .eq. 0) then 
         call appext(fhead,'Agmesh',fwname)
         open (ifwAdapt,file=fwname)
	 do ik = 1, mgk 
	    write(ifwAdapt,9200)xmk(1,ik),xmk(2,ik)
         enddo
	 call flush(ifwAdapt)
	 close (ifwAdapt)
      endif
      if(1 .eq. 1) then
         call appext(fhead,'DispX',fwname)
         open (ifwAdapt,file=fwname)
         do iDispbc = 1, nDispbcX
	    jpt = LmDispbcX(iDispbc)
	    write(ifwAdapt,9202) 
     &      jpt, xm(1,jpt),xm(2,jpt),
     &      dxm(1,jpt),dxm(2,jpt),vDispbcX(iDispbc)
	 enddo
	 call flush(ifwAdapt)
	 close (ifwAdapt)
      endif
      if(1 .eq. 1) then
         call appext(fhead,'DispY',fwname)
         open (ifwAdapt,file=fwname)
         do iDispbc = 1, nDispbcY
	    jpt = LmDispbcY(iDispbc)
	    write(ifwAdapt,9202)
     &      jpt, xm(1,jpt),xm(2,jpt),
     &      dxm(1,jpt),dxm(2,jpt),vDispbcY(iDispbc)
	 enddo
	 call flush(ifwAdapt)
	 close (ifwAdapt)
      endif
      if(1 .eq. 1) then
         call appext(fhead,'VelX',fwname)
         open (ifwAdapt,file=fwname)
         do iVelbc = 1, nVelbcX
	    jpt = LmVelbcX(iVelbc)
	    write(ifwAdapt,9202)
     &      jpt, xm(1,jpt),xm(2,jpt),
     &      dxm(1,jpt),dxm(2,jpt),vVelbcX(iVelbc)
	 enddo
	 close (ifwAdapt)
      endif
      if(1 .eq. 1) then
         call appext(fhead,'VelY',fwname)
         open (ifwAdapt,file=fwname)
         do iVelbc = 1, nVelbcY
	    jpt = LmVelbcY(iVelbc)
	    write(ifwAdapt,9202)
     &      jpt, xm(1,jpt),xm(2,jpt),
     &      dxm(1,jpt),dxm(2,jpt),vVelbcY(iVelbc)
	 enddo
	 call flush(ifwAdapt)
	 close (ifwAdapt)
      endif
c
c.......Update nelem.......
c
	nelem = nelem_new
	numnp = numnp_new
	print *, 'new nelem', nelem
	print *, 'new numnp', numnp
c
c.......Update dt, nstep, nOutputFrq
c
	if (iGL .eq. 0) then
	   dt         = 0.5 * dt
	   nstep      = nstep 
	   nstep      = 2 * nstep 
           nOutputFrq = int(nstep/10)
	elseif(iGL .eq. 1) then
	   dt         = 0.5 * dt
	   nstep      = nstep 
           nstep      = 2 * nstep
           nOutputFrq = int(nstep/10)
	endif
c
c
 97      format(2(1x,e14.6))
 9200    format(2(1x,e12.5))
 9201    format(i5,2(1x,e12.5))
 9202    format(i5,5(1x,e12.5))
c
         return
	 end
c
c
