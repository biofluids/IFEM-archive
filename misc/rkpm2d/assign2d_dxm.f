       subroutine assign2d_dxm(dxm,xm,lnods,nelem,numnp)
c
c-----------------------------------------------------
c
c     Assign the support size for non-uniform
c     particle distribution
c
c
c
c     Shaofan Li
c
c     February, 1999
c
c-----------------------------------------------------
c
       implicit none
       include 'parameter.h'
c
       real*8 xm(2,maxNumnp),dxm(2,maxNumnp)
       integer lnods(maxNode,maxElem),
     &         nelem,numnp
c
c
       integer ie,ipt,ng(4),
     &         inode
       real*8 xn(4),yn(4),dxe(4),dye(4),
     &        xx1,xx2,yy1,yy2,
     &        eps
c
       real*8 ax,ay,afact1,afact2
       integer nnc1,nnc2
c
       common /rkpm/ax,ay,afact1,afact2,nnc1,nnc2
c
      eps = 1.0d-10
c
      do ipt = 1, numnp
	 dxm(1,ipt) = 0.0
	 dxm(2,ipt) = 0.0
      enddo
c
      do ie = 1, nelem
c
	 do inode = 1,4
	    ng(inode) = lnods(inode,ie) 
	    xn(inode) = xm(1,ng(inode))
	    yn(inode) = xm(2,ng(inode))
	 enddo
c
c......(1)................
c
            xx1 = dabs(xn(2) - xn(1))
	    xx2 = dabs(xn(4) - xn(1))
            yy1 = dabs(yn(2) - yn(1))
	    yy2 = dabs(yn(4) - yn(1))
c
            dxe(1) = dmax1(xx1,xx2)
            dye(1) = dmax1(yy1,yy2)
c
c......(2)....................
c
            xx1 = dabs(xn(3) - xn(2))
	    xx2 = dabs(xn(1) - xn(2))
            yy1 = dabs(yn(3) - yn(2))
	    yy2 = dabs(yn(1) - yn(2))
c
            dxe(2) = dmax1(xx1,xx2)
            dye(2) = dmax1(yy1,yy2)
c
c......(3)....................
c
            xx1 = dabs(xn(4) - xn(3))
	    xx2 = dabs(xn(2) - xn(3))
            yy1 = dabs(yn(4) - yn(3))
	    yy2 = dabs(yn(2) - yn(3))
c
            dxe(3) = dmax1(xx1,xx2)
            dye(3) = dmax1(yy1,yy2)
c
c......(4)....................
c
            xx1 = dabs(xn(1) - xn(4))
	    xx2 = dabs(xn(3) - xn(4))
            yy1 = dabs(yn(1) - yn(4))
	    yy2 = dabs(yn(3) - yn(4))
c
            dxe(4) = dmax1(xx1,xx2)
            dye(4) = dmax1(yy1,yy2)
c
c
c...........dxe,dye are ready
c
            do inode = 1,4
c
c........(X)
c
               if (dxm(1,ng(inode)) .lt. eps) then
		   dxm(1,ng(inode)) = ax * dxe(inode)
	       elseif(dxm(1,ng(inode)) .lt. ax*dxe(inode)) then
		   dxm(1,ng(inode)) = 0.5*(ax * dxe(inode)
     &                              + dxm(1,ng(inode)))
	       endif
c
c.........(Y)
c
               if (dxm(2,ng(inode)) .lt. eps) then
		   dxm(2,ng(inode)) = ay * dye(inode)
	       elseif(dxm(2,ng(inode)) .lt. ay*dye(inode)) then
		   dxm(2,ng(inode)) = 0.5*(ay * dye(inode)
     &                              + dxm(2,ng(inode)))
	       endif
c
            enddo
c
c   
       enddo ! ie
c
c
       return
       end
c
