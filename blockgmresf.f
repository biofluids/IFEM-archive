c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	S. Aliabadi                                                          c
c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine blockgmresf(xloc,uloc,qloc,p,hk,ien)

	implicit none
	include "global.h"

      integer ien(nen,nec)
      real* 8 xloc(nsd,nn_loc),uloc(nsd,nn_loc)
      real* 8 p(nn_loc),qloc(nn_loc),hk(nec)

      real* 8 x(nsdpad,nenpad),q(nenpad)
      real* 8 du(nsdpad,nenpad)

      real* 8 eft0,det
      real* 8 sh(0:nsdpad,nenpad),ph(0:nsdpad,nenpad)
      real* 8 xr(nsdpad,nsdpad),cf(nsdpad,nsdpad),sx(nsdpad,nsdpad)

	real* 8 qrt,qrx,qry,qrz
	real* 8 u,v,w
	real* 8 hg,ka,vel
	real* 8 res_f,temp,dtinv
	integer inl, ie, isd, idf, iq, node

	dtinv = 1.0/dt/alpha
	if(steady) dtinv = 0.0

       do ie=1,nec 

       do inl=1,nen
       do isd=1,nsd
       x(isd,inl) =  xloc(isd,ien(inl,ie))
       du(isd,inl) =  uloc(isd,ien(inl,ie))
       enddo
       q(inl) =  qloc(ien(inl,ie))
       enddo

        hg = hk(ie)

        do iq=1,nquad

        if (nen.eq.4) then
                include "sh3d4n.h"
        else if (nen.eq.8) then
                include "sh3d8n.h"
        end if

            eft0 = abs(det) * wq(iq) * alpha


        u = 0.0
	  v = 0.0
	  w = 0.0
        qrx = 0.0
	  qry = 0.0
	  qrz = 0.0
        qrt = 0.0
	  do inl=1,nen
	      qrx=qrx+sh(1,inl)*q(inl)          
	      qry=qry+sh(2,inl)*q(inl)          
	      qrz=qrz+sh(3,inl)*q(inl)          
	      qrt=qrt+sh(0,inl)*q(inl)*dtinv
	      u=u+sh(0,inl)*du(udf,inl)
	      v=v+sh(0,inl)*du(vdf,inl)
	      w=w+sh(0,inl)*du(wdf,inl)
	  end do

        res_f = qrt+u*qrx+v*qry+w*qrz
        vel  = sqrt(u*u+v*v+w*w)
	  ka = delta(3)*0.25*0.5*hg*vel

       do inl=1,nen
       node = ien(inl,ie)
       p(node)=p(node)+res_f*sh(0,inl)*eft0
       p(node)=p(node)+
     &         ka*(qrx*sh(1,inl)+qry*sh(2,inl)+qrz*sh(3,inl))*eft0
       enddo

      enddo
      enddo

      return
      end
