c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	S. Aliabadi                                                          c
c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	  subroutine blockfnode(xloc,shrk,dn,don,ien,cnn,ncnn,hn,hm,fn)
	  
	  implicit none
	  include "global.h"
	  include "malloc.h"

	  integer ien(nen,nec),cnn(maxconn,nqdc),ncnn(nqdc)
	  real* 8 xloc(nsd,nn_loc),fn(nnc),dn(nnc)
	  real* 8 shrk(0:nsd,maxconn,nquad*nec)
	  real* 8 don(nn_on)
	  
	  
	  
	  real* 8 frk(nquad*nec),floc(nn_loc),p(nn_loc),q(nn_loc),hm(nn_loc)
	  real* 8 pn(nnc),qn(nnc),finc(nnc),hn(nnc)
	  pointer (frkptr,frk),(flocptr,floc),(pptr,p),(qptr,q)
	  pointer (pnptr,pn),(qnptr,qn),(fincptr,finc)
	  
	  real* 8 x(nsdpad,nenpad),d(maxconn),f(nenpad)
	  
	  real* 8 eft0,det
	  real* 8 sh(0:nsdpad,nenpad)
	  real* 8 xr(nsdpad,nsdpad),cf(nsdpad,nsdpad),sx(nsdpad,nsdpad)
	  
	  real* 8 frg,res_f,xnorm
	  
	  integer inl,isd,ie,iq,node,qp,itr,ifp,maxrecl,ibase
	  
	  character*9 dout
	  integer i,j,i1,i2,i3,i4
	  integer lk,ir,status(MPI_STATUS_SIZE)
	  logical assemble
	  
	  flocptr = malloc(nn_loc*fsize) 
	  pptr = malloc(nn_loc*fsize) 
	  qptr = malloc(nn_loc*fsize) 
	  pnptr = malloc(nnc*fsize) 
	  qnptr = malloc(nnc*fsize) 
	  fincptr = malloc(nnc*fsize)
	  frkptr = malloc(nec*nquad*fsize)
	  
	  fn(:) = dn(:)
	  call gather (floc,fn,1,hn,hm)
	  
	  qp = 0
	  q(:) = 0
	  do ie = 1,nec

		do inl = 1,nen
		  do isd = 1,nsd
			x(isd,inl) = xloc(isd,ien(inl,ie))
		  end do
		end do

		do iq = 1,nquad

		  qp = qp + 1
		  do inl=1,ncnn(qp)
			d(inl) =  don(cnn(inl,qp))
		  enddo

		  frk(qp) = 0.0
		  do inl=1,ncnn(qp)
			frk(qp)=frk(qp)+shrk(0,inl,qp)*d(inl)
		  end do

		  if (nen.eq.4) then
			include "sh3d4n.h"
		  else if (nen.eq.8) then
			include "sh3d8n.h"
		  end if
		  
		  eft0 = abs(det) * wq(iq)

		  do inl=1,nen
			node = ien(inl,ie)
			q(node)=q(node)+sh(0,inl)*eft0 
		  enddo

		end do
	  end do
	  call scatter (q,qn,1,.true.,hn,hm)



      do itr = 1,50
		p(:) = 0.0
		qp = 0
		do ie=1,nec 
		  
		  do inl = 1,nen
			do isd = 1,nsd
			  x(isd,inl) = xloc(isd,ien(inl,ie))
			end do
		  end do
		  
		  
		  do iq=1,nquad
			qp = qp + 1
			do inl=1,nen
			  f(inl) =  floc(ien(inl,ie))
			enddo
			
			
			
			if (nen.eq.4) then
			  include "sh3d4n.h"
			else if (nen.eq.8) then
			  include "sh3d8n.h"
			end if
			
			
			eft0 = abs(det) * wq(iq)
			
c  Now use RKPM shape functions
			
			
			frg = 0.0
			do inl = 1,nen
			  frg=frg+sh(0,inl)*f(inl)
			enddo
			
			
			res_f = (frg - frk(qp))*eft0

			do inl=1,nen
			  node = ien(inl,ie)
			  p(node)=p(node)-res_f*sh(0,inl)
			enddo
		  enddo
		  
		enddo
		call scatter (p,pn,1,.true.,hn,hm)
		
		do i=1,nnc
		  finc(i) =  pn(i)/qn(i)
		  fn(i) = fn(i) + 1.8*finc(i)
		enddo
		
		call gather (floc,fn,1,hn,hm)
		call getnorm(finc,finc,nnc,xnorm)
	  enddo                
	  
	  call free(flocptr) 
	  call free(pptr) 
	  call free(qptr) 
	  call free(pnptr) 
	  call free(qnptr) 
	  call free(fincptr)
	  call free(frkptr)
	  
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	  
      ibase = ichar('0')		!! integer value for char '0'
	  
	  dout = "dddd.0000"
	  
      i4 = idisk/1000
      i3 = (idisk-i4*1000)/100
      i2 = (idisk-i4*1000-i3*100)/10
      i1 = (idisk-i4*1000-i3*100-i2*10)/1
      
      i4 = i4 + ibase
      i3 = i3 + ibase
      i2 = i2 + ibase
      i1 = i1 + ibase

      dout(6:6) = char(i4)
      dout(7:7) = char(i3)
      dout(8:8) = char(i2)
      dout(9:9) = char(i1)
ccccccccc
	  maxrecl = maxnnc * 8
      lk = 1
      if (myid.gt.0) then
		call MPI_RECV(lk,1,MPI_INTEGER,myid-1,5,MPI_COMM_WORLD,status,ir)
	  end if
	  call ewd_open(dout, ifp)
	  call ewd_lseek(ifp, myid*maxrecl, 0)
	  call ewd_write(ifp,fn,nnc*8)
	  call ewd_close(ifp)
      if (myid.lt.numproc-1) then
		call MPI_SEND(lk,1,MPI_INTEGER,myid+1,5,MPI_COMM_WORLD,ir)
	  end if
ccccccccc

	  return
	  end




