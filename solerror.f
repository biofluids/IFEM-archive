
	subroutine solerror (xn,er,res,pr,ien,eg,fg)

	implicit none
	include "global.h"
	real*8 xn(nsd,nn_loc),er(nsd,nn_loc),err(nsd,nnc),res(nsd,nnc)
	real*8 pr(nsd,nn_loc),eg(nsd,nnc),fg(nsd,nnc)
	integer ien(nen,ne)
	real*8 norma,normb, hn(nnc),hm(nn_loc)
	integer inn,isd

	call fclear(er,nn_loc*nsd)
	call fclear(err,nnc*nsd)
	call fclear(pr,nn_loc*nsd)
c	call berror1(xn,fi,pr,ien)
	call scatter(pr,res,nsd,.true.,hn,hm)
c	call addit(pr,res,ien,nsd)     !! dr is recycled for residual


 100  continue
	call berror2(xn,er,pr,ien)
	call scatter(pr,fg,nsd,.true.,hn,hm)
c	call addit(pr,fg,ien,nsd)
	call getnorm(er,er,nsd*nn_loc,normb)
c	stop
	call gather(er,err,nsd,.true.,hn,hm)
	do inn=1,nnc
	   do isd=1,nsd
	      err(isd,inn)=err(isd,inn)+(res(isd,inn)-fg(isd,inn))/eg(isd,inn)
	   enddo
	   write(*,*) err(1,inn),err(2,inn),err(3,inn)
	enddo
c	stop
	call getnorm(err,err,nsd*nnc,norma)
	norma = norma + epsilon
	stop
	write(6,*) abs(normb/norma-1.0)

	if(abs(normb/norma-1.0).gt.0.25) goto 100


	normb = 150*150*0.5

	do inn=1,nnc
	   err(1,inn) = -abs(err(1,inn))/normb *9.0       
	   err(2,inn) = -abs(err(2,inn))/normb *9.0                    
	   err(3,inn) = -abs(err(3,inn))/normb *9.0                    
	   write(*,*) err(1,inn),err(2,inn),err(3,inn)
	enddo

	call scatter(er,err,nsd,.true.,hn,hm)
	return
	end

c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine berror2(xn, dr, p, ien)

	implicit none
	include "global.h"

	integer ien(nen,ne)
	real* 8 xn(nsd,nn_loc),dr(nsd,nn_loc),p(nsd,nn_loc)
	real* 8 x(3,8),d(3,8),grs(3)

	real* 8 eft0,det,sh(0:3,8),xr(3,3),cf(3,3),sx(3,3)

	integer inl, isd, ie, iq, iflag,node

	iflag = 0

      do ie=1,nec
	 do inl=1,nen
	    do isd=1,nsd
	       x(isd,inl) = xn(isd,ien(inl,ie))
	    enddo
	    d(1,inl) = dr(1,ien(inl,ie))
	    d(2,inl) = dr(2,ien(inl,ie))
	    d(3,inl) = dr(3,ien(inl,ie))
	   enddo

	   do iq=1,nquad
	      if (nen.eq.4) then
		 include "sh3d4n.h"
	      else if (nen.eq.8) then
		 include "sh3d8n.h"
	      end if

	      eft0 = abs(det) * wq(iq)

	      grs(1) = 0.0
	      grs(2) = 0.0
	      grs(3) = 0.0
	      do inl=1,nen
		 grs(1)=grs(1)+sh(0,inl)*d(1,inl)        
		 grs(2)=grs(2)+sh(0,inl)*d(2,inl)        
		 grs(3)=grs(3)+sh(0,inl)*d(3,inl)        
	      end do

	      do inl=1,nen
		 sh(0,inl) = sh(0,inl)*eft0
	      enddo

	      do inl=1,nen
		 node=ien(inl,ie)
		 p(1,node)=p(1,node)+grs(1)*sh(0,inl)
		 p(2,node)=p(2,node)+grs(2)*sh(0,inl)
		 p(3,node)=p(3,node)+grs(3)*sh(0,inl)
	      enddo
	      
	   enddo
	enddo

	return
	end
