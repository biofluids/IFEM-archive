c  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  S. Aliabadi                                                          c
c  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	  subroutine blockgfrk(xloc,shrk,shrkb,qon,p,hk,
	1	   ien,rng,cnn,ncnn)

	  implicit none
	  include "global.h"

      integer ien(nen,nec),cnn(maxconn,nqdc),ncnn(nqdc)
      integer rng(neface,nec)
	  real* 8 xloc(nsd,nn_loc),shrk(0:nsd,maxconn,nquad*nec)
	  real* 8 shrkb(maxconn,neface*nec)
      real* 8 qon(nn_on)
	  real* 8 p(nn_on),hk(nec)

	  real* 8 x(nsdpad,nenpad),q(maxconn)
      real* 8 du(nsdpad,maxconn)

      real* 8 eft0,det
      real* 8 sh(0:nsdpad,nenpad)
      real* 8 xr(nsdpad,nsdpad),cf(nsdpad,nsdpad),sx(nsdpad,nsdpad)

	  real* 8 qrt,qrx,qry,qrz
	  real* 8 res_fi,vi,area,gg,coef 
	  real* 8 a1,a2,a3,b1,b2,b3,c1,c2,c3
	  real* 8 u,v,w
	  real* 8 hg,ka,vel
	  real* 8 res_f,temp,dtinv,beta
	  integer inl,ie,isd,idf,iq,node,qp,qb,node1,node2,node3,node4
	  integer ieface,irng

	  dtinv = 1.0/dt/alpha
	  if(steady) dtinv = 0.0
	  coef = sqrt(hmin*hmin + hmax*hmax)

	  qp = 0
	  qb = 0
      do ie=1,nec 

		do inl = 1,nen
		  do isd = 1,nsd
			x(isd,inl) = xloc(isd,ien(inl,ie))
		  end do
		end do


		do iq=1,nquad
		  qp = qp + 1

		  do inl=1,ncnn(qp)
			q(inl) =  qon(cnn(inl,qp))
	      enddo

		  hg = hk(ie)


		  if (nen.eq.4) then
			include "sh3d4n.h"
		  else if (nen.eq.8) then
			include "sh3d8n.h"
		  end if

		  

		  eft0 = abs(det) * wq(iq) * alpha

c  Now use RKPM shape functions

		  u = 1.0
		  v = 0.0
		  w = 0.0
		  qrx = 0.0
		  qry = 0.0
		  qrz = 0.0
		  qrt = 0.0
		  do inl=1,ncnn(qp)
			temp = q(inl)
			qrx=qrx+shrk(1,inl,qp)*temp
			qry=qry+shrk(2,inl,qp)*temp            
			qrz=qrz+shrk(3,inl,qp)*temp            
			qrt=qrt+shrk(0,inl,qp)*temp*dtinv
		  end do

		  res_f = qrt+u*qrx+v*qry+w*qrz
		  vel = sqrt(u*u+v*v+w*w)
		  ka = delta(3)*0.5*hg*vel

		  do inl=1,ncnn(qp)
			node = cnn(inl,qp)
			p(node)=p(node)+res_f*shrk(0,inl,qp)*eft0
			p(node)=p(node)+
	1			 ka*(qrx*shrk(1,inl,qp)+
     &   		 qry*shrk(2,inl,qp)+qrz*shrk(3,inl,qp))*eft0
		  enddo

		enddo


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccInside the element loop
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        
        do ieface=1,neface
		  irng = rng(ieface,ie)
		  if(irng.gt.0) then
			qb = qb + 1
			if(bcf(irng).eq.1) then
			  node1 = ien(map(ieface,1,etype),ie)
			  node2 = ien(map(ieface,2,etype),ie)
			  node3 = ien(map(ieface,3,etype),ie)
			  a1 = xloc(1,node2)-xloc(1,node1)
			  a2 = xloc(2,node2)-xloc(2,node1)
			  a3 = xloc(3,node2)-xloc(3,node1)
			  b1 = xloc(1,node3)-xloc(1,node1)
			  b2 = xloc(2,node3)-xloc(2,node1)
			  b3 = xloc(3,node3)-xloc(3,node1)

			  c1 =   a2*b3-a3*b2
			  c2 = -(a1*b3-a3*b1)
			  c3 =   a1*b2-a2*b1
			  area = 0.5*sqrt(c1*c1+c2*c2+c3*c3)
			  if(nen.eq.8) then
				node4 = ien(map(ieface,4,etype),ie)
				a1 = xloc(1,node4)-xloc(1,node1)
				a2 = xloc(2,node4)-xloc(2,node1)
				a3 = xloc(3,node4)-xloc(3,node1)

				c1 =   a2*b3-a3*b2
				c2 = -(a1*b3-a3*b1)
				c3 =   a1*b2-a2*b1
				area = area + 0.5*sqrt(c1*c1+c2*c2+c3*c3)
			  endif
              

			  do inl=1,ncnn(nec*nquad+qb)
				q(inl) =  qon(cnn(inl,nec*nquad+qb))
			  enddo

			  vi = 0.0
			  do inl=1,ncnn(nec*nquad+qb)
				vi=vi+shrkb(inl,qb)*q(inl)
			  end do

			  res_fi = coef*vi

			  do inl=1,ncnn(nec*nquad+qb)
				node = cnn(inl,nec*nquad+qb)
				p(node)=p(node)+res_fi*shrkb(inl,qb)*area
			  enddo

			endif
		  endif
		enddo

	  enddo

	  return
	  end
