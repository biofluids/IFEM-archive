c  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  S. Aliabadi                                                          c
c  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	  subroutine locate(xloc, ien, floc) 

      implicit none
	  include "global.h"

      integer ien(nen,nec)
	  real* 8 xloc(nsd,nn_loc),floc(nn_loc)
      real* 8 x(nsdpad,nenpad),f(nenpad),fi,fia,fiam

      real* 8 eft0,det
      real* 8 sh(0:nsdpad,nenpad)
      real* 8 xr(nsdpad,nsdpad),cf(nsdpad,nsdpad),sx(nsdpad,nsdpad)

      real* 8 a,b,c,d,e,m,r,rhs,fpr,y,dy,eps
	  real* 8 bb,cc,dd,ee
	  integer i,inl,ie,isd,iq

      integer ir,status(MPI_STATUS_SIZE)

      eps = 0.001
	  a = liq
	  y = 0.5
	  m = 1.0-delta(5)
	  delta(6) = y
	  if(.not.conserve) return

	  do i=1,3

		b = 0.0
		c = 0.0
		d = 0.0
		e = 0.0

		do ie=1,nec
		  do inl=1,nen
			do isd=1,nsd
			  x(isd,inl) = xloc(isd,ien(inl,ie))
			enddo
			f(inl) = floc(ien(inl,ie))
		  enddo

		  do iq=1,nquad
			if (nen.eq.4) then
			  include "sh3d4n.h"
			else if (nen.eq.8) then
			  include "sh3d8n.h"
			end if
			eft0 = abs(det) * wq(iq)  

			fi = 0.0
			do inl=1,nen
			  fi = fi+ sh(0,inl)*f(inl)
			enddo

			fia = fi**delta(5)
			fiam = (1-fi)**delta(5)
			if(fi.gt.1.0-eps) b = b + fi*eft0
			if((fi.gt.0.0).and.(fi.le.y)) c = c + fia*eft0
			if((fi.gt.y).and.(fi.lt.1.0-eps)) d = d + eft0
			if((fi.gt.y).and.(fi.lt.1.0-eps)) e = e + fiam*eft0
		  enddo
		enddo


		call MPI_BARRIER(MPI_COMM_WORLD,ir)
		call MPI_REDUCE (b,bb,1,MPI_DOUBLE_PRECISION,
	1		 MPI_SUM,0,MPI_COMM_WORLD,ir)
		call MPI_REDUCE (c,cc,1,MPI_DOUBLE_PRECISION,
	1		 MPI_SUM,0,MPI_COMM_WORLD,ir)
		call MPI_REDUCE (d,dd,1,MPI_DOUBLE_PRECISION,
	1		 MPI_SUM,0,MPI_COMM_WORLD,ir)
		call MPI_REDUCE (e,ee,1,MPI_DOUBLE_PRECISION,
	1		 MPI_SUM,0,MPI_COMM_WORLD,ir)
		call MPI_BCAST(bb,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ir)
		call MPI_BCAST(cc,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ir)
		call MPI_BCAST(dd,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ir)
		call MPI_BCAST(ee,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ir)
		call MPI_BARRIER(MPI_COMM_WORLD,ir)
		
		b = bb
		c = cc
		d = dd
		e = ee
		r = a-d-b
		rhs = r-c*y**m+e*(1-y)**m
		fpr = c*m*y**(m-1) + e*m*(1.0-y)**(m-1)
		dy = rhs/fpr
		y = y + dy
		if(y.lt.0.1) y = 0.1
		if(y.gt.0.9) y = 0.9

	  enddo

	  delta(6) = y

	  if(((y.lt.0.25).or.(y.gt.0.75)).and.(myid.eq.0)) write(6,100) y
	  if(((y.lt.0.25).or.(y.gt.0.75)).and.(myid.eq.0)) write(7,100) y
 100  format('INTERFSCE IS STRUGGLING TO ADJUST AT ', f4.4 
	1	   ,/'INCREASE GMRES ITERATIONS FOR FAST RECOVERY.')

      return
      end
