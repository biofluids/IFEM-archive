ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   Lucy Zhang
c   Calculate mesh velocities for ALE.  7/1/99
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine velocity(xn,xnold,xloc,meshveln,meshvel,refvel,
     +     finv,dloc,ien,hn,hm)

      implicit none
      include "global.h"
      
      integer ien(nen,nec)
      real* 8 xn(nsd,nnc), xnold(nsd,nnc)
      real* 8 meshveln(nsd,nnc),dx,refvel(nsd,nquad,nec)
      real* 8 meshvel(nsd,nn_loc),convel(nsd)
      real* 8 f(nsd,nsd,nquad,nec),finv(nsd,nsd,nquad,nec)
      real* 8 d(ndf,nen),mv(nsd,nen),dloc(ndf,nn_loc),diffv(nsd,nen)
      integer inn, isd,node,i,j,ie,inl,iq
      real* 8 hn(nnc),hm(nn_loc),eft0,sh(0:nsdpad,nenpad)
      real* 8 xr(nsdpad,nsdpad),cf(nsdpad,nsdpad),sx(nsdpad,nsdpad)
      real* 8 x(nsd,nen),xloc(nsd,nn_loc),det
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c... calculate mesh velocities, vhat.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call fclear(refvel,nsd*nquad*nec)
      do inn = 1,nnc
         do isd = 1,nsd
            dx = xn(isd,inn) - xnold(isd,inn)
            meshveln(isd,inn) = dx/dt
         enddo
      enddo
      call gather(meshvel,meshveln,nsd,hn,hm)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c... calculate referential velocities, w
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do ie = 1,nec
         do inl=1,nen
            do isd=1,nsd
               d(isd,inl)=dloc(isd,ien(inl,ie))
               x(isd,inl)=xloc(isd,ien(inl,ie))
               mv(isd,inl)=meshvel(isd,ien(inl,ie))
               diffv(isd,inl)=d(isd,inl)-mv(isd,inl)
            enddo
         enddo
         do iq=1,nquad
            if (nen.eq.4) then
               include "sh3d4n.h"
            elseif (nen.eq.8) then
               include "sh3d8n.h"
            endif
            eft0=abs(det)*wq(iq)

            do isd=1,nsd
               convel(isd)=0.0
            enddo

            do inl=1,nen
               convel(xsd)=convel(xsd)+sh(0,inl)*diffv(xsd,inl)
               convel(ysd)=convel(ysd)+sh(0,inl)*diffv(ysd,inl)
               convel(zsd)=convel(zsd)+sh(0,inl)*diffv(zsd,inl)
            enddo

            refvel(isd,iq,ie)=0
        
            do i = 1,nsd
               do j = 1,nsd
                  refvel(i,iq,ie)=refvel(i,iq,ie)+finv(i,j,iq,ie)*convel(j)
               enddo
            enddo
         enddo
      enddo

      return
      end





