c  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  S. Aliabadi                                                          c
c  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	  subroutine blockfrk(xloc,shrk,shrkb,don,doon,p,q,hk,
	1	   ien,rng,cnn,ncnn)

	  implicit none
	  include "global.h"

      integer ien(nen,nec),cnn(maxconn,nqdc),ncnn(nqdc)
      integer rng(neface,nec)
	  real* 8 xloc(nsd,nn_loc),shrk(0:nsd,maxconn,nquad*nec)
	  real* 8 shrkb(maxconn,neface*nec)
      real* 8 don(nn_on),doon(nn_on)
	  real* 8 p(nn_on),q(nn_on),hk(nec)

	  real* 8 x(nsdpad,nenpad),d(maxconn),do(maxconn)
      real* 8 du(nsdpad,maxconn)

      real* 8 eft0,det
      real* 8 sh(0:nsdpad,nenpad)
      real* 8 xr(nsdpad,nsdpad),cf(nsdpad,nsdpad),sx(nsdpad,nsdpad)

	  real* 8 drt,drx,dry,drz
	  real* 8 res_fi,fi,area,gg,coef 
	  real* 8 a1,a2,a3,b1,b2,b3,c1,c2,c3
	  real* 8 u,v,w
	  real* 8 hg,ka,vel
	  real* 8 res_f,temp,dtinv,beta
	  integer inl,ie,isd,idf,iq,node,qp,qb,node1,node2,node3,node4
	  integer ieface,irng

	  dtinv = 1.0/dt
	  if(steady) dtinv = 0.0
	  beta = 1.0 -alpha
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
			d(inl) =  don(cnn(inl,qp))
			do(inl) = doon(cnn(inl,qp))
		  enddo

		  hg = hk(ie)


		  if (nen.eq.4) then
			include "sh3d4n.h"
		  else if (nen.eq.8) then
			include "sh3d8n.h"
		  end if

		  

		  eft0 = abs(det) * wq(iq)

c  Now use RKPM shape functions

		  u = 1.0
		  v = 0.0
		  w = 0.0
		  drx = 0.0
		  dry = 0.0
		  drz = 0.0
		  drt = 0.0
		  do inl=1,ncnn(qp)
			temp = alpha*d(inl)+beta*do(inl)
			drx=drx+shrk(1,inl,qp)*temp
			dry=dry+shrk(2,inl,qp)*temp            
			drz=drz+shrk(3,inl,qp)*temp            
			drt=drt+shrk(0,inl,qp)*(d(inl)-do(inl))*dtinv
		  end do

		  res_f = drt+u*drx+v*dry+w*drz
		  vel = sqrt(u*u+v*v+w*w)
		  ka = delta(3)*0.5*hg*vel

		  do inl=1,ncnn(qp)
			node = cnn(inl,qp)
			p(node)=p(node)-res_f*shrk(0,inl,qp)*eft0
			p(node)=p(node)-
	1			 ka*(drx*shrk(1,inl,qp)+
	2			 dry*shrk(2,inl,qp)+drz*shrk(3,inl,qp))*eft0
		  enddo

		  do inl=1,ncnn(qp)
			node = cnn(inl,qp)
			temp = u*shrk(1,inl,qp)+v*shrk(2,inl,qp)+w*shrk(3,inl,qp)
			temp = shrk(0,inl,qp)*dtinv + alpha*temp
			q(node) = q(node)+temp*shrk(0,inl,qp)*eft0
			temp = shrk(1,inl,qp)**2+shrk(2,inl,qp)**2+shrk(3,inl,qp)**2
			q(node) = q(node)+ka*alpha*temp*eft0
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
			  gg = bvf(irng)
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
				d(inl) =  don(cnn(inl,nec*nquad+qb))
			  enddo

			  fi = 0.0
			  do inl=1,ncnn(nec*nquad+qb)
				fi=fi+shrkb(inl,qb)*d(inl)
			  end do

			  res_fi = coef*(fi-gg)

			  do inl=1,ncnn(nec*nquad+qb)
				node = cnn(inl,nec*nquad+qb)
				p(node)=p(node)-res_fi*shrkb(inl,qb)*area
			  enddo

			  do inl=1,ncnn(nec*nquad+qb)
				node = cnn(inl,nec*nquad+qb)
				q(node) = q(node)+shrkb(inl,qb)*shrkb(inl,qb)*area
			  enddo
			endif
		  endif
		enddo

	  enddo

	  
	  return
	  end
