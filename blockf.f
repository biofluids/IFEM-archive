c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	S. Aliabadi                                                          c
c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine blockf(xloc,dloc,doloc,uloc,p,q,hk,ien)

	implicit none
	include "global.h"

      integer ien(nen,nec)
      real* 8 xloc(nsd,nn_loc),dloc(nn_loc),doloc(nn_loc)
      real* 8 uloc(nsd,nn_loc)
	real* 8 p(nn_loc),q(nn_loc),hk(nec)

	real* 8 x(nsdpad,nenpad),d(nenpad),do(nenpad)
      real* 8 du(nsdpad,nenpad)

      real* 8 eft0,det
      real* 8 sh(0:nsdpad,nenpad)
      real* 8 xr(nsdpad,nsdpad),cf(nsdpad,nsdpad),sx(nsdpad,nsdpad)

	real* 8 drt,drx,dry,drz
	real* 8 u,v,w
	real* 8 hg,ka,vel
	real* 8 res_f,temp,dtinv,beta
	integer inl, ie, isd, idf, iq, node

	dtinv = 1.0/dt
	if(steady) dtinv = 0.0
	beta = 1.0 -alpha

      do ie=1,nec 

      do inl=1,nen
      do isd=1,nsd
      x(isd,inl)  = xloc(isd,ien(inl,ie))
      du(isd,inl) = uloc(isd,ien(inl,ie))
      enddo
      d(inl) =  dloc(ien(inl,ie))
      do(inl) = doloc(ien(inl,ie))
      enddo

	hg = hk(ie)

        do iq=1,nquad

        if (nen.eq.4) then
                include "sh3d4n.h"
        else if (nen.eq.8) then
                include "sh3d8n.h"
        end if

        eft0 = abs(det) * wq(iq)

        u = 0.0
	  v = 0.0
	  w = 0.0
        drx = 0.0
	  dry = 0.0
	  drz = 0.0
        drt = 0.0
	  do inl=1,nen
        temp = alpha*d(inl)+beta*do(inl)
	  drx=drx+sh(1,inl)*temp
	  dry=dry+sh(2,inl)*temp            
	  drz=drz+sh(3,inl)*temp            
	  drt=drt+sh(0,inl)*(d(inl)-do(inl))*dtinv
	  u=u+sh(0,inl)*du(udf,inl)
	  v=v+sh(0,inl)*du(vdf,inl)
	  w=w+sh(0,inl)*du(wdf,inl)
        end do

        res_f = drt+u*drx+v*dry+w*drz
        vel = sqrt(u*u+v*v+w*w)
	  ka = delta(3)*0.25*0.5*hg*vel

       do inl=1,nen
       node = ien(inl,ie)
       p(node)=p(node)-res_f*sh(0,inl)*eft0
       p(node)=p(node)-
     &         ka*(drx*sh(1,inl)+dry*sh(2,inl)+drz*sh(3,inl))*eft0
       enddo

	 do inl=1,nen
	 node = ien(inl,ie)
	 temp = u*sh(1,inl)+v*sh(2,inl)+w*sh(3,inl)
	 temp = sh(0,inl)*dtinv + alpha*temp
	 q(node) = q(node)+temp*sh(0,inl)*eft0
	 temp = sh(1,inl)**2+sh(2,inl)**2+sh(3,inl)**2
	 q(node) = q(node)+ka*alpha*temp*eft0
	 enddo

      enddo
      enddo

      return
      end
